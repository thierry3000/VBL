/*
 *  CellsSystem-D-3.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 09/01/17.
 *  
 *  Thierry Fredrich 01.02.2018
 *  added TBB support --> use prallel_do()
 *
 */

#include <cassert>
#include <cmath>
#include <tbb/parallel_do.h>
#include <stdexcept>

#if VBL_USE_TUMORCODE
  #include <mwlib/helpers-vec.h>
  #include <common/hdfio.h>
  #include <mwlib/lattice-data.h>
#endif

#if VBL_USE_TUMORCODE
void CellsSystem::Set_Tumorcode_Continuous_lattice(LatticeDataQuad3d &field_ld)
{
  field_ld = field_ld;
}
#endif

Triangulation_3 *DelTri;
#ifndef useSerialApproach

void ApplyGeometricCalculation::operator()(const Triangulation_3::Vertex &item) const
{
  apply_geometry(item);
}


void apply_geometry(const Triangulation_3::Vertex &iter)
{
  Triangulation_3::Vertex_handle a = DelTri->insert(iter.point());
  //if(!DelTri->is_infinite(a))
  {
    unsigned long k_mt = a->info();
    std::vector<Triangulation_3::Vertex_handle> vn_per_thread;
    DelTri->finite_adjacent_vertices(a, std::back_inserter(vn_per_thread));
    //at this stage, the no. of neighbor includes the infinite parts
    p_to_current_CellsSystem->Set_neigh(k_mt, vn_per_thread.size());
    
#ifndef NDEBUG
    printf("I apply geometry calculations on vertex %i\n", a->info());
#endif
    //default of isonCH is true!!!
    p_to_current_CellsSystem->Set_isonCH(k_mt, false);      // bool true if ON convex hull
    //p_to_current_CellsSystem->Set_env_surf(k_mt, 0.);       // contact area with environment
    p_to_current_CellsSystem->Set_g_env(k_mt, 0.);          // geom factor with environment
    p_to_current_CellsSystem->Set_contact_surf(k_mt, 0.0);  // init total contact area
    // weighted radius of current cell
    double rk = (p_to_current_CellsSystem->Get_type(k_mt)->Get_extension_coeff())*p_to_current_CellsSystem->Get_r(k_mt);
#ifndef NDEBUG
    printf("neigh[k_mt]: %i, ncells: %i\n", p_to_current_CellsSystem->Get_neigh(k_mt), p_to_current_CellsSystem->Get_ncells());
#endif
    // in this loop we prepare the list of neighbors and we compute contact areas
    
    for(int kk=0; kk < p_to_current_CellsSystem->Get_neigh(k_mt) ; kk++)					
    {
      auto neighbor_id = vn_per_thread[kk]->info(); 	// id of neighbor
      
#ifndef NDEBUG
      std::cout << "starting second for loop" << std::endl;
      printf(" kk: %i, neigh[k]: %i\n", kk, p_to_current_CellsSystem->Get_neigh(k_mt));
      assert(neighbor_id < p_to_current_CellsSystem->Get_ncells());
      //printf("k_mt: %i, neighbor: %i, kk: %i, vneigh.size(): %i \n", k_mt, neighbor_id, kk, vneigh.size());
      std::cout << "neighbor: " << neighbor_id << std::endl;
#endif

      p_to_current_CellsSystem->Set_vneigh_quick(k_mt, kk, neighbor_id); // store name of neighbor
      
      //***
      // weighted radius of neighboring cell
      double rkk = (p_to_current_CellsSystem->Get_type(neighbor_id)->Get_extension_coeff())*p_to_current_CellsSystem->Get_r(neighbor_id);
      double dd = p_to_current_CellsSystem->Distance( k_mt,neighbor_id ); 			// distance between current cell and neigbor
      // *********** check for debugging
#ifndef NDEBUG // this check is only performed in the Debug build, other wise it is not present which saves time
      if(dd != dd || dd <= 0)
      {
        std::cout << k_mt << "-th cell, neighbor " << neighbor_id << ", undefined distance " << dd << std::endl;
      }
#endif
      // *********** end check for debugging
      p_to_current_CellsSystem->Set_vdist_quick(k_mt, kk, dd);
      double sqr_dd = SQR(dd);

      if( dd < (rk+rkk) ) // here we compute the contact area
      {
        double this_vc_surf= -PI*(sqr_dd-SQR(rk-rkk))*(sqr_dd-SQR(rk+rkk))/(4*sqr_dd);
        p_to_current_CellsSystem->Set_vcsurf_quick(k_mt, kk, this_vc_surf );
      
        if( this_vc_surf < 0 )
        {
          p_to_current_CellsSystem->Set_vcsurf_quick(k_mt, kk, 0.0 );
          p_to_current_CellsSystem->Set_gnk_quick(k_mt, kk, 0.0 );
        }
        else
        {
          p_to_current_CellsSystem->Set_gnk_quick(k_mt, kk, (this_vc_surf/dd) );
        }
        //contact_surf[k_mt] += vcsurf[k_mt][nFV];
        double total_contact_surface = p_to_current_CellsSystem->Get_contact_surf(k_mt);
        p_to_current_CellsSystem->Set_contact_surf(k_mt, (total_contact_surface + this_vc_surf));
      }
      else
      {
        p_to_current_CellsSystem->Set_vcsurf_quick(k_mt, kk, 0.0 );
      }
    }
  
    double this_env_surf = p_to_current_CellsSystem->Get_surface(k_mt)/(p_to_current_CellsSystem->Get_neigh(k_mt)+1);
    p_to_current_CellsSystem->Set_env_surf(k_mt, this_env_surf);
    double this_g_env = this_env_surf/p_to_current_CellsSystem->Get_r(k_mt);
    p_to_current_CellsSystem->Set_g_env(k_mt, this_g_env);
    //vn_per_thread.clear();
#ifndef NDEBUG
    printf("finished apply geometry\n");
#endif
  }
}

void ParallelApply()
{
  tbb::parallel_do(DelTri->finite_vertices_begin(), DelTri->finite_vertices_end(), ApplyGeometricCalculation());
}

#endif

// ************ this section contains the interface to CGAL ************ //
// 
// REMARK: the CGAL methods are used ONLY in this part of the program, and are isolated 
// from the rest of the code to avoid interferences.
//
void CellsSystem::Geometry()
{
#ifdef useSerialApproach
#pragma message(" useSerialApproach defined")
#endif
  unsigned long k;
  
  /* I use a pointer now
   * this allows to either use the parallel code (tbb) or the serial variant.
   * I is important to have both. parallel version works only with a minimum of 5 cells
   * (see tools/tests/info_insert_with_pair_iterator_parallel.cpp )
   */
  //Triangulation_3 *DelTri; // I use a pointer now
  //
  
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::min();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::min();
  double min_z = std::numeric_limits<double>::max();
  double max_z = std::numeric_limits<double>::min();
  // cout << "Triangulation OK" << endl;

  // insertion of cells' centres into the vector of Point structures
  // the vector of Points is defined in CellsSystem.h, and the initial reserve is set in 
  // CellsSystem.cpp

  v.clear(); // vector of all points

//Delaunay triangulation
  
#ifndef useSerialApproach
  if(ncells > 10)
  {
    for(k=0; k<ncells; k++)
    {
#ifndef NDEBUG
      if(std::isnan( x[k]))
        printf("x is nan at %i\n", k);
      if(std::isnan( y[k]))
        printf("y is nan at %i\n", k);
      if(std::isnan( z[k]))
        printf("z is nan at %i\n", k);
      std::cout.flush();
#endif
      v.push_back( std::make_pair(Point(x[k],y[k],z[k]),k) );
    
      // find min, max to estimate lock data structure
      if(x[k]<min_x)
      {
        min_x=x[k];
      }
      if(x[k]>max_x)
      {
        max_x=x[k];
      }
      if(y[k]<min_y)
      {
        min_y=y[k];
      }
      if(y[k]>max_y)
      {
        max_y=y[k];
      }
      if(z[k]<min_z)
      {
        min_z=z[k];
      }
      if(z[k]>max_z)
      {
        max_z=z[k];
      }
    }
#ifndef NDEBUG
    min_x=min_y=min_z=-10.0;
    max_x=max_y=max_z=10.0;
    if(std::isnan(min_x) or std::isnan(min_y) or std::isnan(min_z) or std::isnan(max_x) or std::isnan(max_y) or std::isnan(max_z))
    {
      printf("found nan in bounding box\n");
      std::cout.flush();
    }
#endif
    Triangulation_3::Lock_data_structure locking_ds(CGAL::Bbox_3(min_x, min_y, min_z, max_x, max_y, max_z), 20);
    DelTri = new Triangulation_3( v.begin(),v.end(), &locking_ds );
  }
  else
  {
    for(k=0; k<ncells; k++)
    {
      v.push_back( std::make_pair(Point(x[k],y[k],z[k]),k) );
    }
    DelTri = new Triangulation_3( v.begin(),v.end() );
  }
  
#else  // in the serial case, we do no need a Lock_data_structure

  //std::cout << "found " << ncells << " in the system" << std::endl;
  for(k=0; k<ncells; k++)
  {
#ifndef NDEBUG
    //printf("push back point: %f, %f, %f, %i\n", x[k], y[k], z[k], k);
#endif
    v.push_back( std::make_pair(Point(x[k],y[k],z[k]),k) );
  }
  
  DelTri = new Triangulation_3( v.begin(),v.end() );
  CGAL_assertion( DelTri->number_of_vertices() == ncells );
  assert(DelTri->is_valid());
#endif
  
  //prepare neighborhood and initialize
  for(int i=0; i< ncells; i++)
  {
    //int no_current_neighbors = 25;
    //neigh[i] = no_current_neighbors;          // number of all neighbors finit+ infinite
    //vneigh[i].resize(no_current_neighbors);		// init names of neighbors
    //vcsurf[i].resize(no_current_neighbors);		// init contact areas
    //vdist[i].resize(no_current_neighbors);		// allocate distance vect
    //gnk[i].resize(no_current_neighbors);			// allocate geom factors
    isonCH[i]=true;       // by default we assume the cell is on the convex hull, we change that if we detect that it is not!
    isonAS[i] = false;    // by default we are not on the alph shape
    g_env[k]=0.;
  }

#ifdef W_timing
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //begin neighborhood timing
#endif

#ifdef useSerialApproach
  /**
   * this was the OLD way!
   * serial access to the triangulation -> no memory problems
   * this piece of code stores all Vertex_handles in a container 
   * and marches throught the container. 
   * Getting the neighbor for every entry of the container
   */
  std::vector<Triangulation_3::Vertex_handle> vn_per_thread;
  for (Finite_vertices_iterator vit = DelTri->finite_vertices_begin(); vit != DelTri->finite_vertices_end(); ++vit)
  {
    int k_mt = vit->info();
    DelTri->finite_adjacent_vertices(vit, std::back_inserter(vn_per_thread));
    
    //reset to finite value!
    //k_mt could have nasty values in case of infinite
    int current_nneigh = vn_per_thread.size();
    if(current_nneigh > MAX_N_NEIGH-1)
    {
      std::cout << "current_nneigh: " << current_nneigh << std::endl;
      std::cout << "Error: only " << MAX_N_NEIGH << " neigbors are supported by the array" << std::endl;
      throw std::runtime_error("more neighbours than supported!!!");
    }
    neigh[k_mt] = current_nneigh;
  
    isonCH[k_mt]=false;				// bool true if ON convex hull

    // weighted radius of current cell
    double rk = (type[k_mt]->Get_extension_coeff())*r[k_mt];

    contact_surf[k_mt] = 0.;			// init total contact area
    
    // in this loop we prepare the list of neighbors and we compute contact areas
    for(int kk=0; kk < current_nneigh ; kk++)					
    {
#ifdef myDebugComments
      std::cout << "starting second for loop" << std::endl;
#endif

#ifdef myDebugComments
      printf(" kk: %i, neigh[k]: %i, vn_per_thread.size(): %i\n", kk, neigh[k_mt],vn_per_thread.size());
#endif
      Triangulation_3::Vertex_handle a = vn_per_thread[kk];
      unsigned long neighbor_id = a->info();
      
#ifndef NDEBUG
assert(neighbor_id < ncells);
#endif
#ifdef myDebugComments
printf("k: %i, neighbor: %i, kk: %i, vneigh.size(): %i \n", k, neighbor_id, kk, vneigh.size());
      std::cout << "neighbor: " << neighbor_id << std::endl;
#endif
      vneigh[k_mt][kk] = neighbor_id; 				// store name of neighbor
      //***
      // weighted radius of neighboring cell
      double rkk = (type[neighbor_id]->Get_extension_coeff())*r[neighbor_id];	
      double dd = Distance( k_mt,neighbor_id ); 			// distance between current cell and neigbor
      // *********** check for debugging
      
#ifdef ENABLE_RELEASE_DEBUG_OUTPUT
      if(dd != dd || dd <= 0 || k_mt == neighbor_id)
      {
        std::cout << k_mt << "-th cell, neighbor " << neighbor_id << ", undefined distance " << dd << std::endl;
      }
#endif
      // *********** end check for debugging
      
      vdist[k_mt][kk] = dd;
      double sqr_dd = SQR(dd);

      if( dd < (rk+rkk) ) // here we compute the contact area
      {
        vcsurf[k_mt][kk] = -PI*(sqr_dd-SQR(rk-rkk))*(sqr_dd-SQR(rk+rkk))/(4*sqr_dd);
        if( vcsurf[k_mt][kk] < 0 )
        {
          vcsurf[k_mt][kk] = 0;
          gnk[k_mt][kk] = 0;
        }
        else
        {
          gnk[k_mt][kk] = vcsurf[k_mt][kk]/dd;
        }
        contact_surf[k_mt] += vcsurf[k_mt][kk];
      }
      else
      {
        vcsurf[k_mt][kk] = 0.;
      }
    }
    env_surf[k_mt] = surface[k_mt]/(neigh[k_mt]+1);			// contact area with environment
    g_env[k_mt] = env_surf[k_mt]/r[k_mt];					// geom factor with environment
    

    // *** fine del calcolo del fattore geometrico con l'ambiente
    vn_per_thread.clear();
  }//end for( k = 0;k<ncells; ++k) 

#else

  ParallelApply();
  
#endif
  
#ifdef W_timing
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    myTiming.time_geometry_neighborhood = end-begin;
    myTiming.geometry_neighborhood = myTiming.time_geometry_neighborhood.count();
#endif

  // compute fixed alpha shape with ALPHAVALUE defined in sim.h
  
  // we know apriori, that there are more than 5 cells, otherwise Geometry is
  // not called, see CellsSystem.cpp
//   if( ncells < 5 )	// if there are less than 5 cells they are all on alphashape
//   {
//     for(k=0; k<ncells; k++)
//     {
//       isonAS[k] = true;
//     }
//   }
//   else 
  
  //Fixed_alpha_shape_3 as(v.begin(),v.end(),ALPHAVALUE); // this is not so fast, since we traverse v again!
  Fixed_alpha_shape_3 as((*DelTri),ALPHAVALUE);
  std::vector<Vertex_handle> vertices_on_alpha_shape;
  as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::REGULAR);
  as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::SINGULAR);

  // set boolean true for vertices on alpha shape
  for(k=0; k<vertices_on_alpha_shape.size(); k++)
  {
    isonAS[vertices_on_alpha_shape[k]->info()] = true;
  }
  
  delete DelTri;
  
  // questo loop azzera i g_env di tutte le cellule che non stanno sull'alpha shape
//   for(k=0; k<ncells; k++)
//   {
//     if( !isonAS[k] ) 
//       g_env[k]=0.;
//   }
    
  // the following loop identifies those cells that are in contact with blood vessels
#pragma omp parallel
{
  std::array<double,3> cellpos; // store the cell coordinates in a 3-vector
  unsigned long k_mt;
#pragma omp for
  for(k_mt=0; k_mt<ncells; k_mt++)
  {    
    cellpos[0]=x[k_mt];
    cellpos[1]=y[k_mt];
    cellpos[2]=z[k_mt];

    isonBV[k_mt] = 0; // by default, cells are not close to blood vessels 
    g_bv[k_mt] = 0; // by default, there is no contact term with blood vessels
    
    for( int nvessel=0; nvessel<nbv; nvessel++ ) // loop over all blood vessels
    {
      double x0[3]; // position on blood vessel axis closest to cell
      double dbv = BloodVesselVector[nvessel].DistanceFromVessel( cellpos, x0 ); // distance between cell and blood vessel
      if( dbv < BloodVesselVector[nvessel].GetBloodVesselR() + r[k_mt] ) // if the cell's center and the blood vessel axis are closer than the sum of the radii, then there is contact
      {
        isonBV[k] = nvessel+1;		// note the shift, which is done to use the value 0 as false and otherwise use the true value to store the blood vessel position in the blood vessel vector
        bv_surf[k_mt] = surface[k_mt];
        g_bv[k_mt] = bv_surf[k_mt]/r[k_mt]; // this are bold statements; here we assume that the cell behaves like a disk (not a sphere), and that half of the surface area faces the blood vessel
        // the last approximation should probably be mitigated with the inclusion of a surface-modulating parameter or even better by an improved geometrical modelling
        break;	// blood vessel found, we jump out of the blood vessel loop	
      }
    }
    // g_bv[k] = 0; // **** TEMPORARY, used only to eliminate blood vessels from calculations !!!   
  }
}// #pragma omp parallel
}

// calcoli minimi nel caso di cellule disperse
void CellsSystem::NoGeometry()
{
  for(unsigned long k=0; k<ncells; k++)
  {
    env_surf[k] = surface[k];				// qui si calcola la superficie esposta all'ambiente
    g_env[k] = env_surf[k]/r[k];			// fattore geometrico verso l'ambiente
    contact_surf[k] = 0;					// la sup. di contatto con le altre cellule e' nulla nel caso di cellule disperse
    g_bv[k] = 0;							// no contact with blood vessels
  }
}

int CellsSystem::checkNeighbourhood_consistency(std::string atPlace)
{
  for(int i = 0; i<ncells; ++i)
  {
    for(int k =0; k<neigh[i]; ++k)
    {
      if( neigh[i] > ncells )
      {
	std::cout << "checking consistency at: " << atPlace << std::endl;
	std::cout << "non plausible neigh: " << neigh[i] << " ncells: " << ncells << std::endl;
	return 1;
      }
      if( i == vneigh[i][k] )
      {
	std::cout << "checking consistency at: " << atPlace << std::endl;
	std::cout<< "checkNeighbourhood_consistency faild for cell: " << i << " and neigh #"<< k << " which is: " << vneigh[i][k] << std::endl;
	return 1;
      }
    }
  }
  return 0;
}
