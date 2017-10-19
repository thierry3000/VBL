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
  Dt DelTri;

  // cout << "Triangolazione OK" << endl;

  // vettore dei punti
  // vector<Point> v( ncells );	

  // inserimento del centri delle cellule nel vettore
#pragma omp parallel for
  for(k=0; k<ncells; k++)
  {
    v[k] = Point( x[k], y[k], z[k] );
  }
  // cout << "Setup punti OK" << endl;
  // cout << "Primo punto: {" << v[0] << "} = " << x[0] << ", " << y[0] << ", " << z[0] << endl;


  // inserimento modificato dei punti nella triangolazione (equivalente all'inserimento sequenziale)
  build_triangulation_with_indices(v.begin(),v.end(),DelTri);


  // nel seguito si trova la lista dei collegamenti e si stimano le aree delle superfici di contatto

  std::vector<Vertex_handle> vn;// vettore delle vertex_handle dei vicini

  // si rifa' il loop sui vertici finiti (k e' il nome della cellula attuale)
  for (Finite_vertices_iterator vit = DelTri.finite_vertices_begin(); vit != DelTri.finite_vertices_end(); ++vit)
  {
    k = vit->info();
    DelTri.incident_vertices(vit, back_inserter(vn));		// qui si ottiene la lista dei vicini
    neigh[k] = vn.size();								// numero dei vicini (finiti E infiniti)
			
    vneigh[k].resize(neigh[k]);							// inizializzazione del vettore dei nomi dei vicini
    vcsurf[k].resize(neigh[k]);							// inizializzazione del vettore delle aree delle superfici di contatto
			
    vdist[k].resize(neigh[k]);							// allocazione del vettore delle distanze
    gnk[k].resize(neigh[k]);							// allocazione del vettore dei fattori geometrici
    isonCH[k]=false;									// variabile che dice se la cellula e' sul convex hull
    env_surf[k] = 0.;									// superficie di contatto con l'ambiente
    g_env[k] = 0.;										// fattore geometrico relativo al contatto con l'ambiente
			
    // raggio pesato della cellula considerata
    double rk = (type[k]->Get_extension_coeff())*r[k];
			
    contact_surf[k] = 0.;								// inizializzazione calcolo superficie di contatto totale
			
    int nFV = 0;										// inizializzazione del numero dei vertici finiti adiacenti	
    for(int kk=0; kk < neigh[k] ; kk++)					// in questo loop si prepara la lista dei vicini e si calcolano le sup. di contatto
    {
      if(!DelTri.is_infinite(vn[kk]))	// se si c'e' un vertice finito adiacente
      {
	int neighbor = vn[kk]->info();				// nome del vicino
	vneigh[k][nFV] = neighbor;					// memorizzazione del nome del vicino

	// raggio pesato della cellula adiacente
	double rkk = (type[neighbor]->Get_extension_coeff())*r[neighbor];	
	double dd = Distance( k,neighbor );			// distanza tra le due cellule
					
	// *********** controllo per debugging
	if(dd != dd || dd <= 0)
	{
	  std::cout << "cellula " << k << "-esima, vicina " << neighbor << "-esima, distanza indefinita" << std::endl;
	}
	// *********** fine controllo per debugging

	vdist[k][nFV] = dd;

	if( dd < (rk+rkk) ) 						// calcolo della superficie di contatto
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

					if(vcsurf[k][nFV] > 0)						// fattore geometrico
						gnk[k][nFV] = vcsurf[k][nFV]/dd;
					else
						gnk[k][nFV] = 0;

					//***

					nFV++;										// qui si incrementa il contatore dei vertici finiti
					}
				else
					isonCH[k]=true;								// se uno dei vertici adiacenti e' infinito allora il vertice si trova sul convex hull
				}


			if( nFV != neigh[k] )	// qui si controlla il numero di vertici finiti e se questo e' diverso dal numero totale di vertici
									// si fa un resize dei vettori
				{
				neigh[k] = nFV;
				vneigh[k].resize(neigh[k]);						// resizing del vettore dei nomi dei vicini
				vcsurf[k].resize(neigh[k]);						// resizing del vettore delle aree delle superfici di contatto
				vdist[k].resize(neigh[k]);						// resizing del vettore delle distanze
				gnk[k].resize(neigh[k]);						// resizing del vettore dei fattori geometrici
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
	
	// Calculation of alpha shape as, with ALPHA defined in sim.h

	if( ncells < 5 )	// se ci sono meno di 5 cellule, queste stanno certamente sull'alpha shape
		{
		for(k=0; k<ncells; k++)
			isonAS[k] = true;
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

		// per default i vertici sono interni
		for(k=0; k<ncells; k++)
			isonAS[k] = false;
		
		for(k=0; k<vertices_on_alpha_shape.size(); k++)
			isonAS[vertices_on_alpha_shape[k]->info()] = true;
		
		}
	
	// questo loop azzera i g_env di tutte le cellule che non stanno sull'alpha shape
	for(k=0; k<ncells; k++)
		if( !isonAS[k] ) g_env[k]=0.;
		
	// the following loop identifies those cells that are in contact with blood vessels 
	for(k=0; k<ncells; k++)
	{
	  std::vector<double> cellpos(3); // store the cell coordinates in a 3-vector
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
