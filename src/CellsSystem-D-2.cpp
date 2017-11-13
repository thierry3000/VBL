/*
 *  CellsSystem-D-2.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  This file takes care of the geometry via CGAL
 *
 */
#include "CellsSystem.h"
using namespace vbl;

// ************ Part of the CGAL interface ************ //
// 
// WARNING: CGAL methods are used only in this section of the program, which has been "isolated" from the rest to avoid
// interferences with other parts of the code. 


// The main method for calculating the structure of Delaunay
// This method isolates all calls to CGAL
//
void CellsSystem::Geometry()
{

  unsigned long k;

  // Creation of Delaunay's triangulation structure
  /** segmentation fault, when running in multiple threads... what does that mean?
   */
  Dt DelTri;

  // cout << "Triangolazione OK" << endl;

  // vettore dei punti
  // vector<Point> v( ncells );	

  // Inserting cell cells into the carrier
//#pragma omp parallel for
  for(k=0; k<ncells; k++)
  {
    v[k] = Point( x[k], y[k], z[k] );
  }
  // cout << "Setup punti OK" << endl;
  // cout << "Primo punto: {" << v[0] << "} = " << x[0] << ", " << y[0] << ", " << z[0] << endl;


  // modified point insertion in triangulation (equivalent to sequential insertion)
  build_triangulation_with_indices(v.begin(),v.end(),DelTri);

  // Below is the list of links and the areas of contact surfaces are estimated

  std::vector<Vertex_handle> vn;// vector of neighbors vertex_handle

  // the loop on the finished vertices (k is the name of the current cell)
  for (Finite_vertices_iterator vit = DelTri.finite_vertices_begin(); vit != DelTri.finite_vertices_end(); ++vit)
  {
    k = vit->info();
    DelTri.incident_vertices(vit, back_inserter(vn));		// here is the list of neighbors
    neigh[k] = vn.size();								// number of neighbors (finite and infinite)
			
    vneigh[k].resize(neigh[k]);							// Initializing the carrier of neighbor names
    vcsurf[k].resize(neigh[k]);							// Initialize the vector of contact surface areas
			
    vdist[k].resize(neigh[k]);							// allocation of distances carrier
    gnk[k].resize(neigh[k]);							// allocation of the vector of the geometric factors
    isonCH[k]=false;									// variable that says if the cell is on the convex hull
    env_surf[k] = 0.;									// surface of contact with the environment
    g_env[k] = 0.;										// geometric factor for contact with the environment
			
    // weighed beam of the considered cell
    double rk = (type[k]->Get_extension_coeff())*r[k];
			
    contact_surf[k] = 0.;								// Initialization total contact surface calculation
			
    int nFV = 0;										// initialization of the number of adjacent finite vertices
    for(int kk=0; kk < neigh[k] ; kk++)					// in this loop you prepare the list of neighbors and calculate the shoulders. of contact
    {
      if(!DelTri.is_infinite(vn[kk]))	// if there is an adjacent end vertex
      {
        int neighbor = vn[kk]->info();				// neighbor's name
        vneigh[k][nFV] = neighbor;					// memorizing the neighbor's name

        // weighed radius of the adjacent cell
        double rkk = (type[neighbor]->Get_extension_coeff())*r[neighbor];	
        double dd = Distance( k,neighbor );			// distance between the two cells
					
        // *********** controllo per debugging
        if(dd != dd || dd <= 0)
        {
          std::cout << "cellula " << k << "-esima, vicina " << neighbor << "-esima, distanza indefinita" << std::endl;
        }
        // *********** fine controllo per debugging

        vdist[k][nFV] = dd;

        if( dd < (rk+rkk) ) 						// calculation of contact surface
        {
          vcsurf[k][nFV] = -PI*(SQR(dd)-SQR(rk-rkk))*(SQR(dd)-SQR(rk+rkk))/(4*SQR(dd));
          if( vcsurf[k][nFV] < 0 )
          {
            vcsurf[k][nFV] = 0;
          }
          contact_surf[k] += vcsurf[k][nFV];
        }
        else
        {
          vcsurf[k][nFV] = 0.;
        }

        if(vcsurf[k][nFV] > 0)						// geometric factor
          gnk[k][nFV] = vcsurf[k][nFV]/dd;
        else
          gnk[k][nFV] = 0;

					//***
        // here the counter of the finite vertices is increased
        nFV++;
        }
      else
      {
        isonCH[k]=true;// if one of the adjacent vertices is infinite then the vertex is on the convex hull
      }
    }

/* here you control the number of finite vertices and if this is different from the total number of vertices
 * you make a resize of the carriers
 */
    if( nFV != neigh[k] )
    {
      neigh[k] = nFV;
      vneigh[k].resize(neigh[k]);						// resizing del vettore dei nomi dei vicini
      vcsurf[k].resize(neigh[k]);						// resizing del vettore delle aree delle superfici di contatto
      vdist[k].resize(neigh[k]);						// resizing del vettore delle distanze
      gnk[k].resize(neigh[k]);						// resizing del vettore dei fattori geometrici
    }
/** start calculating geometric environmental factor
 * calculating the surface of contact with the environment: 
 * it is assumed that the entire surface of the cell that is not in contact 
 * with adjacent cells with the environment: this calculation is done however, but the
 * contact with the environment is in fact only the cell is on the alpha shape
 */
    env_surf[k] = 0.;								// normally the shoulder. contacting the environment is nothing
    g_env[k] = 0.;									// and of course the associated geometric factor is also null
			
			
//			env_surf[k] = surface[k] - contact_surf[k];		// qui si calcola la superficie esposta all'ambiente
//			if( env_surf[k] < 0 )
//				env_surf[k] = 0.;							// in linea di principio (causa algoritmo approssimato) questa superficie esposta 
															// puo' essere negativa, e in questo caso la si annulla

/** new calculation of the surface exposed to the environment: 
 * calculate the number of neighbors + the environment and assume that the surface surface [k]
 * be equally divided between neighbors and the environment
 */
    env_surf[k] = surface[k]/(neigh[k]+1);			// here we calculate the surface exposed to the environment
    g_env[k] = env_surf[k]/r[k];					// geometric factor towards the environment
			
// *** the end of the geometric factor calculation with the environment
			
			vn.clear();											// the neighbors list is cleared in preparation for the next summit

  }
	
	// Calculation of alpha shape as, with ALPHA defined in sim.h

	if( ncells < 5 )	// if there are less than 5 cells, these are certainly on the alpha shape
  {
    for(k=0; k<ncells; k++)
    {
      isonAS[k] = true;
    }
  }
	else 
  {
    Alpha_shape_3 as(DelTri);
    as.set_mode(Alpha_shape_3::GENERAL);

    // Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
    // std::cout << "Smallest alpha value to get a solid through data points is " << scientific << alpha_solid << std::endl;

    std::vector<Vertex_handle> vertices_on_alpha_shape;
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Alpha_shape_3::REGULAR,ALPHAVALUE);
    as.get_alpha_shape_vertices(std::back_inserter(vertices_on_alpha_shape),Alpha_shape_3::SINGULAR,ALPHAVALUE);
    
    
    // cout << "there are " << vertices_on_alpha_shape.size() << " vertices on alpha shape " << endl;

    // by default, the vertices are internal
    for(k=0; k<ncells; k++)
      isonAS[k] = false;
    
    for(k=0; k<vertices_on_alpha_shape.size(); k++)
      isonAS[vertices_on_alpha_shape[k]->info()] = true;
  }
	
	// this loop resets the g_env of all cells that are not on the alpha shape
	for(k=0; k<ncells; k++)
		if( !isonAS[k] ) 
      g_env[k]=0.;
		
	// the following loop identifies those cells that are in contact with blood vessels 
	for(k=0; k<ncells; k++)
	{
	  //std::vector<double> cellpos(3); // store the cell coordinates in a 3-vector
    std::array<double,3> cellpos;
	  cellpos[0]=x[k];
	  cellpos[1]=y[k];
	  cellpos[2]=z[k];
        
	  isonBV[k] = 0; // by default, cells are not close to blood vessels 
	  g_bv[k] = 0; // by default, there is no contact term with blood vessels
		
	  for( int nvessel=0; nvessel<nbv; nvessel++ ) // loop over all blood vessels
	  {
	    double x0[3]; // position on blood vessel axis closest to cell
	    double dbv = BloodVesselVector[nvessel].DistanceFromVessel( cellpos, x0 ); // distance between cell and blood vessel
	    //double dbv = bloodVesselMap[nvessel].DistanceFromVessel( cellpos, x0 );
	    // if the cell's center and the blood vessel axis are closer than the sum of the radii, then there is contact
	    if( dbv < BloodVesselVector[nvessel].GetBloodVesselR() + r[k] )
	    //if( dbv < bloodVesselMap[nvessel].GetBloodVesselR() + r[k] )
	    {
	      /** note the shift, which is done to use the value 0
	       * as false and otherwise use the true value 
	       * to store the blood vessel position in the blood vessel vector
	       */
	      isonBV[k] = nvessel+1;
	      bv_surf[k] = surface[k];
	      /** this are bold statements; 
	       * here we assume that the cell behaves like a disk (not a sphere), 
	       * and that half of the surface area faces the blood vessel
	       */
	      g_bv[k] = bv_surf[k]/r[k];
	      /** the last approximation should probably be mitigated with the inclusion 
	       * of a surface-modulating parameter 
	       * or even better by an improved geometrical modelling
	       */
	      break;	// blood vessel found, we jump out of the blood vessel loop	
	    }
	  }
	  // g_bv[k] = 0; // **** TEMPORARY, used only to eliminate blood vessels from calculations !!! 
	}
}

// Minimum calculations in case of dispersed cells
void CellsSystem::NoGeometry()
{
  for(unsigned long k=0; k<ncells; k++)
  {
    env_surf[k] = surface[k];				// Here we calculate the surface exposed to the environment
    g_env[k] = env_surf[k]/r[k];			// Geometric factor towards the environment
    contact_surf[k] = 0;					// Sup. Contacting the other cells is nothing in the case of dispersed cells
  }


}
