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
#include <CGAL/natural_neighbor_coordinates_3.h>
#include <chrono>

Triangulation_3 *DelTri;

#if VBL_USE_TUMORCODE


// void CellsSystem::interpolate_O2_uptake_to_tumorcode_2(CellBasedO2Uptake &o2_uptake_model, std::vector<double> &O2Rates)
// {
//   auto t1 = std::chrono::high_resolution_clock::now();
//   float n_cells = (float) Get_ncells();
//   std::vector<double> x = Get_x();
//   std::vector<double> y = Get_y();
//   std::vector<double> z = Get_z();
//   #pragma omp parallel
//   {
//     BOOST_FOREACH(const BBox3 &bbox, o2_uptake_model.mtboxes->getCurrentThreadRange())
//     {
//       for(int i=0; i<x.size();++i)
//       {
//         float offset = 1020.0;
//         Float3 pos(x[i]+offset,y[i]+offset,z[i]+offset);
//         //AddSmoothDelta(o2_uptake_model.o2_consumption_field, bbox, o2_uptake_model.grid->ld, o2_uptake_model.grid->Dim(), pos, (float)O2Rates[i]);
//       }
//     }
//   }
//   auto t2 = std::chrono::high_resolution_clock::now();
//   auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
//   std::cout << "interpolation took: " << duration << " ms." << std::endl;
// }
// 
// void interpolate_O2_uptake_to_tumorcode(CellBasedO2Uptake &o2_uptake_model, std::vector<double> &O2Rates)
// {
//   typedef std::vector< std::pair< Triangulation_3::Vertex_handle, K::FT  > >                                                      Vertex_handle_vector;
//   
//   
// #define inter_parallel
//   
// #if 0
//   Vertex_handle_vector coords;
//   double norm_coeff;
//   Int3 p(3,3,3);
//   Point current_lattice_point_in_world_coordinates = Point(p[0]*o2_uptake_model.grid->Spacing(), p[1]*o2_uptake_model.grid->Spacing(), p[2]*o2_uptake_model.grid->Spacing());
//   CGAL::Triple<std::back_insert_iterator<Vertex_handle_vector>,K::FT, bool> result = CGAL::sibson_natural_neighbor_coordinates_3(*DelTri, current_lattice_point_in_world_coordinates,std::back_inserter(coords),norm_coeff);
// #endif
//   /** this actually worked!!!!
//    * but this way I only find the 4 cells next to single grid point,
//    * which neglects all the other cells
//    * I think this approximation is too rude
//    */
//   
//   //Float3 offset = o2_uptake_model.grid->ld.GetOriginPosition();
//   Float3 offset(-980.,-1020.0, -880.0);
//   float scale = 40.0;
//   //printf("offset 1: %f, offset 2: %f, offset 3: %f\n", offset[0], offset[1], offset[2]);
// #ifdef inter_parallel
//   #pragma omp parallel
//   {
// #endif
//     BOOST_FOREACH(const BBox3 &bbox, o2_uptake_model.mtboxes->getCurrentThreadRange())
//     {
//       FOR_BBOX3(p, bbox)
//       {
//         Vertex_handle_vector coords;
//         double norm_coeff;
//         //Point current_lattice_point_in_world_coordinates = Point(p[0]*o2_uptake_model.grid->Spacing()-offset[0], p[1]*o2_uptake_model.grid->Spacing()-offset[1], p[2]*o2_uptake_model.grid->Spacing()-offset[2]);
//         Point current_lattice_point_in_cell_coordinates = Point(p[0]*scale+offset[0], p[1]*scale+offset[1], p[2]*scale+offset[2]);
// 
//         Float3 pospos(current_lattice_point_in_cell_coordinates[0],current_lattice_point_in_cell_coordinates[1],current_lattice_point_in_cell_coordinates[2]);
//         //AddSmoothDelta(o2_uptake_model.o2_consumption_field, bbox, o2_uptake_model.grid->ld, o2_uptake_model.grid->Dim(), pospos, (float) 42.);
//         
// #if 0 // proper interpolation takes for ages, I try the average nearest neighbor
//         /**
//         * execute sibson_natural neighbors coordinates
//         */
//         CGAL::Triple<std::back_insert_iterator<Vertex_handle_vector>,K::FT, bool> result = CGAL::sibson_natural_neighbor_coordinates_3(*DelTri, current_lattice_point_in_world_coordinates,std::back_inserter(coords),norm_coeff);
//         if(!result.third)
//         {
//           //std::cout << "The coordinate computation was not successful." << std::endl;
//           //std::cout << "The point (" <<current_lattice_point_in_world_coordinates << ") lies outside the convex hull."<< std::endl;
//           o2_uptake_model.fieldLastSources(p) = 0;
//         }
//         else
//         {
//           K::FT interpolated_value(0);
//           for(auto vertexAndFTPair : coords )
//           {
//             auto cellId = vertexAndFTPair.first->info();
//             float o2_consumption_at_cell = O2Rates[cellId];
//             interpolated_value+= o2_consumption_at_cell * vertexAndFTPair.second;
//           }
//           o2_uptake_model.fieldLastSources(p) = interpolated_value / norm_coeff;
//         }
// #endif
//         //locate cell with grid point in it
//         Triangulation_3::Cell_handle c = DelTri->locate(current_lattice_point_in_cell_coordinates);
//         float localOxygenConsumption= 0.0;
//         if(! DelTri->is_infinite(c) )
//         {
//           for(int i=0;i<4;++i)
//           {
//             localOxygenConsumption += O2Rates[c->vertex(i)->info()];
//           }
//           //localOxygenConsumption = localOxygenConsumption/4.;
//           localOxygenConsumption = 4000.;
//           //std::cout << "oxygen: " << localOxygenConsumption << std::endl;
//         }
//         //o2_uptake_model.o2_consumption_field(p) = localOxygenConsumption;
//         //o2_uptake_model.o2_consumption_field(p) = 42.0;
//         
//         //Int3 ip; Float3 q;
//         //boost::tie(ip, q) = o2_uptake_model.grid->ld.WorldToFractionalCoordinate(pos);
// //           float last = fieldLastSources(p);
// //           float src = sources(p);
// //           float diff = src-last;
// //           if(std::abs(diff)>1.0e-2)
// //           {
// //             {
// //               _conv_addfunc f; f.s = diff;
// //               array_apply_offset3d<float, float,_conv_addfunc>(
// //                 gffield,
// //                 gf_lut,
// //                 f,
// //                 p-(gf_lut.size()/2),
// //                 0
// //               );
// //             }
// //             fieldLastSources(p) = src;
// //           }
//         //std::cout << current_lattice_point_in_world_coordinates << std::endl;
//       }
//     }
// #ifdef inter_parallel
//   }
// #endif
//   
//   //I want to check with the center
// //   std::cout << "interpolate here" << std::endl;
// //   std::cout << o2_uptake_model.fieldLastSources(0,0,0) << std::endl;
// //   std::cout << "change" << std::endl;
// //   o2_uptake_model.fieldLastSources(0,0,0) = 42;
// //   std::cout << o2_uptake_model.fieldLastSources(0,0,0) << std::endl;
// 
// }


#endif


void CellsSystem::Geometry()
{

	unsigned long k;
	
	// cout << "Triangulation OK" << endl;
	
// insertion of cells' centres into the vector of Point structures
// the vector of Points is defined in CellsSystem.h, and the initial reserve is set in 
// CellsSystem.cpp

	v.clear();	

// #pragma omp parallel for
	for(k=0; k<params.ncells; k++)
		v.push_back( std::make_pair(Point(x[k],y[k],z[k]),k) ); 
	
	// cout << "v size: " << v.size() << endl;

	// cout << "Setup punti OK" << endl;
	// cout << "Primo punto: {" << v[0] << "} = " << x[0] << ", " << y[0] << ", " << z[0] << endl;

//Delaunay triangulation
	Triangulation_3 DelTri( v.begin(),v.end() );

// next we find the list of connections and we estimate contact areas

	std::vector<Vertex_handle> vn; // vector of Vertex_handle's of neighbors

// loop over finite vertices (k is the name of current cell)
	for (Finite_vertices_iterator vit = DelTri.finite_vertices_begin(); vit != DelTri.finite_vertices_end(); ++vit)
			{
			
			k = vit->info();
			
			DelTri.incident_vertices(vit, back_inserter(vn));	// list of neighbors
			neigh[k] = vn.size();								// number of all neighbors
			
			//vneigh[k].resize(neigh[k]);							// init names of neighbors
			//vcsurf[k].resize(neigh[k]);							// init contact areas
			
			//vdist[k].resize(neigh[k]);							// allocate distance vect
			//gnk[k].resize(neigh[k]);							// allocate geom factors
			isonCH[k]=false;									// bool true if ON convex 
																// hull
			env_surf[k] = 0.;									// contact area with 
																// environment
			g_env[k] = 0.;										// geom factor with 
																// environment
			
			// weighted radius of current cell
			double rk = (type[k]->Get_extension_coeff())*r[k];
			
			contact_surf[k] = 0.;								// init total contact area
			
			int nFV = 0;										// init number of finite 
																// vertices

// in this loop we prepare the list of neighbors and we compute contact areas 
			if(neigh[k] > 0)
			{
			for(int kk=0; kk < neigh[k] ; kk++)					
				{
				
				if(!DelTri.is_infinite(vn[kk]))	// if true, then there is a neighboring 
												// finite vertex
					{
					
					int neighbor = vn[kk]->info(); // name of neighbor
					vneigh[k][nFV] = neighbor; // store name of neighbor
					
					//***

					// weighted radius of neighboring cell
					double rkk = (type[neighbor]->Get_extension_coeff())*r[neighbor];	
					double dd = Distance( k,neighbor ); // distance between current cell 
														// and neighbor
					
					// *********** check for debugging
					if(dd != dd || dd <= 0)
						{
						std::cout << k << "-th cell, neighbor " << neighbor << ", undefined distance " << dd << std::endl;
						}
					// *********** end check for debugging
					
					vdist[k][nFV] = dd;
					
					if( dd < (rk+rkk) ) // here we compute the contact area
						{
						vcsurf[k][nFV] = -PI*(SQR(dd)-SQR(rk-rkk))*(SQR(dd)-SQR(rk+rkk))/(4*SQR(dd));
						if( vcsurf[k][nFV] < 0 ) vcsurf[k][nFV] = 0;
						contact_surf[k] += vcsurf[k][nFV];
						}
					else
						vcsurf[k][nFV] = 0.;

					if(vcsurf[k][nFV] > 0)						// geometric factor
						gnk[k][nFV] = vcsurf[k][nFV]/dd;
					else
						gnk[k][nFV] = 0;

					//***

					nFV++;	// here we increase the counter of finite vertices
					}
				else
					isonCH[k]=true;	// if any one of the adjacent vertices is infinite 
									// then cell is on convex hull
				}
			}


			if( nFV != neigh[k] )	// qui si controlla il numero di vertici finiti e se questo e' diverso dal numero totale di vertici
									// si fa un resize dei vettori
				{
				neigh[k] = nFV;
				//vneigh[k].resize(neigh[k]);						// resizing del vettore dei nomi dei vicini
				//vcsurf[k].resize(neigh[k]);						// resizing del vettore delle aree delle superfici di contatto
				//vdist[k].resize(neigh[k]);						// resizing del vettore delle distanze
				//gnk[k].resize(neigh[k]);						// resizing del vettore dei fattori geometrici
				}
			
		// *** inizio calcolo fattore geometrico ambientale
			// calcolo della superficie di contatto con l'ambiente: si assume che tutta la superficie della cellula che non e' a
			// contatto con le cellule adiacenti sia in contatto con l'ambiente: questo calcolo viene fatto comunque, ma il 
			// contatto con l'ambiente c'e' in realta' solo le la cellula sta sull'alpha shape
			
			env_surf[k] = 0.;								// normalmente la sup. di contatto con l'ambiente e' nulla
			g_env[k] = 0.;									// e naturalmente anche il fattore geometrico associato e' nullo
			
			
//			env_surf[k] = surface[k] - contact_surf[k];		// qui si calcola la superficie esposta all'ambiente
//			if( env_surf[k] < 0 )
//				env_surf[k] = 0.;							// in linea di principio (causa algoritmo approssimato) questa superficie esposta 
															// puo' essere negativa, e in questo caso la si annulla


// nuovo calcolo della superficie esposta all'ambiente: si calcola il numero di vicini + l'ambiente e si assume che la superficie surface[k]
// sia equamente suddivisa tra vicini e ambiente
			env_surf[k] = surface[k]/(neigh[k]+1);			// qui si calcola la superficie esposta all'ambiente

			g_env[k] = env_surf[k]/r[k];					// fattore geometrico verso l'ambiente
			
		// *** fine del calcolo del fattore geometrico con l'ambiente
			
			
			vn.clear();											// si ripulisce la lista dei vicini in preparazione del prossimo vertice

			}
	
	// compute fixed alpha shape with ALPHAVALUE defined in sim.h

	if( params.ncells < 5 )	// if there are less than 5 cells they are all on alphashape
  {
    std::cout << "less than 5 cells" << std::endl;
		for(k=0; k<params.ncells; k++)
			isonAS[k] = true;
  }
	else 
  {
    //auto o2Rates = Get_O2Rate();
		Fixed_alpha_shape_3 as(v.begin(),v.end(),ALPHAVALUE);
		std::vector<Vertex_handle> vertices_on_alpha_shape;
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::REGULAR);
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Fixed_alpha_shape_3::SINGULAR);
	
		// cout << "there are " << vertices_on_alpha_shape.size() << " vertices on alpha shape " << endl;

		// initialize boolean indicator variable
		for(k=0; k<params.ncells; k++)
			isonAS[k] = false;
		
		// set boolean true for vertices on alpha shape
		for(k=0; k<vertices_on_alpha_shape.size(); k++)
    {
			isonAS[vertices_on_alpha_shape[k]->info()] = true;
      //std::cout << vertices_on_alpha_shape[k]->info() << " -> true " << std::endl;
    }
		
  }
	
	// this loop resets the g_env of all cells that are not on the alpha shape
	for(k=0; k<params.ncells; k++)
		if( !isonAS[k] ) 
      g_env[k]=0.;
		
	// the following loop identifies those cells that are in contact with blood vessels
#pragma omp parallel
{
  std::array<double,3> cellpos; // store the cell coordinates in a 3-vector
  unsigned long k_mt;
#pragma omp for
  for(k_mt=0; k_mt<params.ncells; k_mt++)
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
        isonBV[k_mt] = nvessel+1;		// note the shift, which is done to use the value 0 as false and otherwise use the true value to store the blood vessel position in the blood vessel vector
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
  for(unsigned long k=0; k<params.ncells; k++)
  {
    env_surf[k] = surface[k];				// qui si calcola la superficie esposta all'ambiente
    g_env[k] = env_surf[k]/r[k];			// fattore geometrico verso l'ambiente
    contact_surf[k] = 0;					// la sup. di contatto con le altre cellule e' nulla nel caso di cellule disperse
    g_bv[k] = 0;							// no contact with blood vessels
  }
}

int CellsSystem::checkNeighbourhood_consistency(std::string atPlace)
{
  for(int i = 0; i<params.ncells; ++i)
  {
    for(int k =0; k<neigh[i]; ++k)
    {
      if( neigh[i] > params.ncells )
      {
	std::cout << "checking consistency at: " << atPlace << std::endl;
	std::cout << "non plausible neigh: " << neigh[i] << " ncells: " << params.ncells << std::endl;
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
