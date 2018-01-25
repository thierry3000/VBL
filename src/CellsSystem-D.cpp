/*
 *  CellsSystem-D-3.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 09/01/17.
 *
 */

#include "CellsSystem.h"

using namespace vbl;

// ************ this section contains the interface to CGAL ************ //
// 
// REMARK: the CGAL methods are used ONLY in this part of the program, and are isolated 
// from the rest of the code to avoid interferences.
//
void CellsSystem::Geometry()
{
  unsigned long k;
  /* I use a pointer now
   * this allows to either use the parallel code (tbb) or the serial variant.
   * I is important to have both. parallel version works only with a minimum of 5 cells
   * (see tools/tests/info_insert_with_pair_iterator_parallel.cpp )
   */
  //Triangulation_3 *DelTri; // I use a pointer now
  std::shared_ptr<Triangulation_3> DelTri;
  
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
  //arr_of_vn_pointers.clear();
  //arr_of_vn_pointers.resize(ncells);

//Delaunay triangulation
  if(ncells > 10)
  {
    for(k=0; k<ncells; k++)
    {
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
    Triangulation_3::Lock_data_structure locking_ds(CGAL::Bbox_3(min_x, min_y, min_z, max_x, max_y, max_z), 20);
    //DelTri = new Triangulation_3( v.begin(),v.end(), &locking_ds );
    DelTri = std::make_shared<Triangulation_3>( v.begin(),v.end(), &locking_ds );
  }
  else
  {
    for(k=0; k<ncells; k++)
    {
      v.push_back( std::make_pair(Point(x[k],y[k],z[k]),k) );
    }
    //DelTri = new Triangulation_3( v.begin(),v.end() );
    DelTri = std::make_shared<Triangulation_3>( v.begin(),v.end() );
  }
  /* Use pointer here!
   * this makes sure that the memory space is safely freed after this function has finished
   */
  CGAL::Compact_container<Vertex_handle> *myContainer = new CGAL::Compact_container<Vertex_handle>();
  
// next we find the list of connections and we estimate contact areas
//#define myDebugComments
// loop over finite vertices (k is the name of current cell)
// note: we do not now how many vertices are finite at this stage
// therefore we need to do this in a single thread
  const int MAX_NEIGHBOR = 15;
  int MAX_NEIGH = 20;
  //for (Finite_vertices_iterator vit = DelTri->finite_vertices_begin(); vit != DelTri->finite_vertices_end(); ++vit)
  for (All_vertices_iterator vit = DelTri->all_vertices_begin(); vit != DelTri->all_vertices_end(); ++vit)
  {
    k = vit->info();
//     arr_of_vn_pointers[k] = std::shared_ptr<std::vector<Vertex_handle>>(new std::vector<Vertex_handle>);
//     DelTri->incident_vertices(vit, std::back_inserter(*arr_of_vn_pointers[k])); // list of neighbors
//     
//     //these are the not thread safe operations!
//     neigh[k] = (*arr_of_vn_pointers[k]).size();// number of all neighbors
//     vneigh[k].resize(neigh[k]);							// init names of neighbors
//     vcsurf[k].resize(neigh[k]);							// init contact areas
//     vdist[k].resize(neigh[k]);							// allocate distance vect
//     gnk[k].resize(neigh[k]);							// allocate geom factors
    myContainer->insert(vit);
    //k = vit->info();
    //arr_of_vn_pointers[k] = std::shared_ptr<std::vector<Vertex_handle>>(new std::vector<Vertex_handle>);
    //DelTri->incident_vertices(vit, std::back_inserter(*arr_of_vn_pointers[k])); // list of neighbors
    
    //these are the not thread safe operations!
    //neigh[k] = (*arr_of_vn_pointers[k]).size();// number of all neighbors
    vneigh[k].resize(ncells);							// init names of neighbors
    vcsurf[k].resize(ncells);							// init contact areas
    vdist[k].resize(ncells);							// allocate distance vect
    gnk[k].resize(ncells);							// allocate geom factors
    //isonCH[k] = false;
    //env_surf[k] = 0.;				// contact area with environment
    //g_env[k] = 0.;				// geom factor with environment
  }
//   for(int t =0;t<ncells;++t)
//   {
//     vneigh[t].resize(MAX_NEIGHBOR);							// init names of neighbors
//     vcsurf[t].resize(MAX_NEIGHBOR);							// init contact areas
//     vdist[t].resize(MAX_NEIGHBOR);							// allocate distance vect
//     gnk[t].resize(MAX_NEIGHBOR);							// allocate geom factors
//   }
  //printf("ncells: %i, myContainer.size(): %i\n", ncells, myContainer.size());
  //assert(ncells == myContainer.size());
#ifdef myDebugComments
  std::cout << "starting geometry in parallel" << std::endl;
  std::cout.flush();
#endif
//#pragma omp parallel
{
  std::vector<Vertex_handle> vn_per_thread;
//#pragma omp for
  for( unsigned long cnter = 0; cnter <myContainer->size(); ++cnter)
  {
    unsigned long k_mt = (*myContainer)[cnter]->info();
    //DelTri->incident_vertices((*myContainer)[cnter], std::back_inserter(vn_per_thread));
    DelTri->finite_incident_vertices((*myContainer)[cnter], std::back_inserter(vn_per_thread));
    neigh.at(k_mt) = vn_per_thread.size();
    printf("neigh[k_mt]: %i, ncells: %i\n", neigh[k_mt], ncells);
    assert(neigh[k_mt]<=ncells);// there should be at most ncells neighbors!
    
    if( !DelTri->is_infinite((*myContainer)[cnter]) )
    {
  #if 0
      vneigh[k].resize(neigh[k]);							// init names of neighbors
      vcsurf[k].resize(neigh[k]);							// init contact areas
      vdist[k].resize(neigh[k]);							// allocate distance vect
      gnk[k].resize(neigh[k]);							// allocate geom factors
  #endif
      isonCH[k_mt]=false;				// bool true if ON convex hull
      env_surf[k_mt] = 0.;				// contact area with environment
      g_env[k_mt] = 0.;				// geom factor with environment

      // weighted radius of current cell
      double rk = (type[k_mt]->Get_extension_coeff())*r[k_mt];

      contact_surf[k_mt] = 0.;			// init total contact area
      
      unsigned long nFV = 0;				// init number of finite vetices
      // in this loop we prepare the list of neighbors and we compute contact areas 
      for(int kk=0; kk < neigh[k_mt] ; kk++)					
      {
  #ifdef myDebugComments
        std::cout << "starting second for loop" << std::endl;
  #endif
        if(!DelTri->is_infinite(vn_per_thread[kk]))	// if true, then there is a neighboring finite vertex
        {
  #ifdef myDebugComments
          printf(" kk: %i, neigh[k]: %i\n", kk, neigh[k_mt]);
  #endif
          auto neighbor_id = vn_per_thread[kk]->info(); 	// id of neighbor
#ifndef NDEBUG
    assert(neighbor_id < ncells);
#endif
  #ifdef myDebugComments
    printf("k: %i, neighbor: %i, nFV: %i, vneigh.size(): %i \n", k, neighbor, nFV, vneigh.size());
          std::cout << "neighbor: " << neighbor << std::endl;
  #endif
          vneigh[k_mt][nFV] = neighbor_id; 				// store name of neighbor
          //***
          // weighted radius of neighboring cell
          double rkk = (type[neighbor_id]->Get_extension_coeff())*r[neighbor_id];	
          double dd = Distance( k_mt,neighbor_id ); 			// distance between current cell and neigbor
          // *********** check for debugging
        #ifndef NDEBUG // this check is only performed in the Debug build, other wise it is not present which saves time
          if(dd != dd || dd <= 0)
          {
            std::cout << k_mt << "-th cell, neighbor " << neighbor_id << ", undefined distance " << dd << std::endl;
          }
        #endif
          // *********** end check for debugging
          vdist[k_mt][nFV] = dd;
          double sqr_dd = SQR(dd);

          if( dd < (rk+rkk) ) // here we compute the contact area
          {
            
            vcsurf[k_mt][nFV] = -PI*(sqr_dd-SQR(rk-rkk))*(sqr_dd-SQR(rk+rkk))/(4*sqr_dd);
            if( vcsurf[k_mt][nFV] < 0 )
            {
              vcsurf[k_mt][nFV] = 0;
              gnk[k_mt][nFV] = 0;
            }
            else
            {
              gnk[k_mt][nFV] = vcsurf[k_mt][nFV]/dd;
            }
            contact_surf[k_mt] += vcsurf[k_mt][nFV];
          }
          else
          {
            vcsurf[k_mt][nFV] = 0.;
          }
          // NOTE is saves time to this in one if
          // 	if(vcsurf[k][nFV] > 0)						// geometric factor
          // 	{
          // 		gnk[k][nFV] = vcsurf[k][nFV]/dd;
          // 	}
          // 	else
          // 	{
          // 		gnk[k][nFV] = 0;
          // 	}
                //***
        nFV++;	// here we increase the counter of finite vertices
        }
        else
        {
    isonCH[k_mt]=true;	// if any one of the adjacent vertices is infinite 
                  // then cell is on convex hull
        }
      }
    }

#if 0 // not possible in paralle section
    if( nFV != neigh[k] )	// here we check the number of vertices finished and if this is different from the total number of vertices
						// a resize of the vectors is made
    {
      neigh[k] = nFV;
      vneigh[k].resize(neigh[k]);						// resizing del vettore dei nomi dei vicini
      vcsurf[k].resize(neigh[k]);						// resizing del vettore delle aree delle superfici di contatto
      vdist[k].resize(neigh[k]);						// resizing del vettore delle distanze
      gnk[k].resize(neigh[k]);						// resizing del vettore dei fattori geometrici
    }
#endif 
    // *** inizio calcolo fattore geometrico ambientale
    // calcolo della superficie di contatto con l'ambiente: si assume che tutta la superficie della cellula che non e' a
    // contatto con le cellule adiacenti sia in contatto con l'ambiente: questo calcolo viene fatto comunque, ma il 
    // contatto con l'ambiente c'e' in realta' solo le la cellula sta sull'alpha shape

    env_surf[k_mt] = 0.;								// normalmente la sup. di contatto con l'ambiente e' nulla
    g_env[k_mt] = 0.;									// e naturalmente anche il fattore geometrico associato e' nullo


    //			env_surf[k] = surface[k] - contact_surf[k];		// qui si calcola la superficie esposta all'ambiente
    //			if( env_surf[k] < 0 )
    //				env_surf[k] = 0.;							// in linea di principio (causa algoritmo approssimato) questa superficie esposta 
												// puo' essere negativa, e in questo caso la si annulla


    // nuovo calcolo della superficie esposta all'ambiente: si calcola il numero di vicini + l'ambiente e si assume che la superficie surface[k]
    // sia equamente suddivisa tra vicini e ambiente
    env_surf[k_mt] = surface[k_mt]/(neigh[k_mt]+1);			// qui si calcola la superficie esposta all'ambiente

    g_env[k_mt] = env_surf[k_mt]/r[k_mt];					// fattore geometrico verso l'ambiente

    // *** fine del calcolo del fattore geometrico con l'ambiente

    vn_per_thread.clear();											// si ripulisce la lista dei vicini in preparazione del prossimo vertice
    //arr_of_vn_pointers[k]->clear();
  }//end for( k = 0;k<ncells; ++k) end PARALLEL
}//#pragma omp parallel


  myContainer->clear();
  // compute fixed alpha shape with ALPHAVALUE defined in sim.h
  if( ncells < 5 )	// if there are less than 5 cells they are all on alphashape
  {
    for(k=0; k<ncells; k++)
    {
      isonAS[k] = true;
    }
  }
  else 
  {
    Fixed_alpha_shape_3 as(v.begin(),v.end(),ALPHAVALUE);
    std::vector<Vertex_handle> vertices_on_alpha_shape;
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::REGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::SINGULAR);
    // cout << "there are " << vertices_on_alpha_shape.size() << " vertices on alpha shape " << endl;
    // initialize boolean indicator variable
    for(k=0; k<ncells; k++)
    {
	  isonAS[k] = false;
    }
    // set boolean true for vertices on alpha shape
    for(k=0; k<vertices_on_alpha_shape.size(); k++)
    {
      isonAS[vertices_on_alpha_shape[k]->info()] = true;
    }
  }
  // questo loop azzera i g_env di tutte le cellule che non stanno sull'alpha shape
  for(k=0; k<ncells; k++)
  {
    if( !isonAS[k] ) g_env[k]=0.;
  }
    
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
