/*
 *  CellsSystem-E.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  This file takes care of mechanical forces
 *
 */


// ***************************************************************
// Cluster mechanics
// ***************************************************************
//
// Cell dynamics (this method calculates the shift in a timestep)
//
// The algorithm used by this method corresponds to a semi-implicit Euler integration (explicit on forces, implied
// on the dissipative part; Should approximate approximately to an implicit verlet method)
// 
// Note that in the current version, Dynamics is not an interface for the GetForces and NewPositionsAndVelocities methods
//
void CellsSystem::Dynamics( )
{
// first part, initialization and calculation of the force
// the calculation of force is done once and for all, corresponds to the explicit part of the integration

	GetForces( );
	
	// calcolo di nuove velocità e posizioni
	NewPositionsAndVelocities( );
}

// Calculating force between cells in the current configuration
void CellsSystem::GetForces()
{
	// resizing vectors
	fx.resize(ncells);
	fy.resize(ncells);
	fz.resize(ncells);
	
	// initialize vectors with zeros
	fx.assign(ncells,0);
	fy.assign(ncells,0);
	fz.assign(ncells,0);
		
	// the cell-cell force is calculated only if there are at least two cells
	if(ncells>1)
  {
#pragma omp parallel for
		for(unsigned long n=0; n<ncells; n++)
    {
			
			// parameters related to the n-th cell
			double rn = (type[n]->Get_packing_factor())*r[n];			// beam of the weighed cell with the packing factor
			double yn = (type[n]->Get_YoungMod());								// modulo di Young della cellula, NOTE fixed value -> no divide by zero here
			double pn = (type[n]->Get_PoissonRatio());						// PoissonRatio della cellula
			
			double arn = (type[n]->Get_adhesion_range());	        // adhesion range of the n-th cell type
			double adn = (type[n]->Get_adhesion_decay());	        // decay rate of k-th cell adhesion

			int nneigh = neigh[n];				                        // number of cells adjacent to the n-th cell
			
			// loop on cells adjacent to the n-th cell (it is calculated only if there is at least one cell nearby)
			if(nneigh > 0)
      {
				for( int k=0; k< nneigh; k++)					
        {
          /* NOTE T.F.
           * I had some strange SIGFPE on GetForces, so I try to check the requriremts here
           * 
           * 1) ((1-pn*pn)/yn+(1-pk*pk)/yk)) == 0 ?
           * if I calculted correctly with p_n ==p_k and y_k==y_n this could only happen when p_n*p_n == 1 which never happens
           * 
           * 2) (rn+rk) == 0?
           * 
           * 3) dd == 0?
           */
					int neighbor = vneigh[n][k];	                                  // name of the k-th adjacent cell
					double dd = Distance(n,neighbor);	                              // distance between the n-th cell and the k-th nearest cell

					// parameters related to the near k-esima
					double rk = (type[neighbor]->Get_packing_factor())*r[neighbor];	// radius of the nearest k-th weighed
					double yk = (type[neighbor]->Get_YoungMod());			              // modulo di Young della cellula NOTE fixed value -> no divide by zero here
					double pk = (type[neighbor]->Get_PoissonRatio());		            // PoissonRatio della cellula

					double ark = (type[neighbor]->Get_adhesion_range());	          // range di adesione del tipo della vicina k-esima
					double adk = (type[neighbor]->Get_adhesion_decay());	          // decay rate dell'adesione della vicina k-esima

					double kC = sqrt(rn*rk)*(rn+rk)/(0.75*((1-pn*pn)/yn+(1-pk*pk)/yk));	// calculation of the "elastic constant" derived from the Hertz model
					
          /*
           * Here the force is calculated, but when it is repulsive saturated to a max.
           * Here you choose to define the max value on the base of cellular characteristics.
           * An alternative choice could be to define a priori the max value for all cells.
           * What is the right choice? Here I take the first, arbitrarily.
           */
					const double saturation_distance = 0.25;                       // saturation distance of the repulsive force (micron)
					// const double saturation_distance = 0.03;							       // saturation distance of the repulsive force
					
					//checks
					if(rn+rk == 0.)
          {
            std::cout << "very problematic: rn+rk == 0" << std::endl;
          }
          if( dd == 0.)
          {
            std::cout << "very problematic: dd == 0" << std::endl;
          }
					
					double fm = kC * pow(fabs((rn+rk-dd)/(rn+rk)),(double)1.5);		 // force form
					if(rn+rk-dd > saturation_distance) 
						fm = kC * pow( saturation_distance/(rn+rk), (double)1.5);		 // the repulsive force saturated at a rather small value ...
					
          // in the case where the distance is greater than the sum of the rays, the force is attractive
					if( dd > (rn+rk) )
						fm *= -0.5*( 1. - tanh( 0.5*(adn+adk)*( dd/(rn+rk)-(1.+0.5*(arn+ark)) ) ) );	// the force changes sign and is modulated
          
          /* NOTE T.F.
           * this looks rather thread safe, but what if dd is zero here?
           */
					fx[n] += fm*(x[n]-x[neighbor])/dd;								// the projection of fm is added to each of the components of the total force
					fy[n] += fm*(y[n]-y[neighbor])/dd;
					fz[n] += fm*(z[n]-z[neighbor])/dd;

          
          
		// sezione per il dump diagnostico, scommentare se serve
		/*
					cout << "\n*************************\n" << endl;
					cout << "forza di interazione tra le cellule: " << n << " e " << neighbor << endl;
					cout << "in posizione (" << scientific << x[n] << ", " << y[n] << ", " << z[n];
					cout << ");    (" << x[neighbor] << ", " << y[neighbor] << ", " << z[neighbor] << ");    " << endl;
					cout << "raggi cellulari: " << scientific << rn << "; " << rk << endl;
					cout << "range di adesione: " << scientific << arn << "; " << ark << endl;
					cout << "decay rate adesione: " << scientific << arn << "; " << ark << endl;
					cout << "distanza tra cellule: " << scientific << dd << endl;
					cout << "kC: " << scientific << kC << endl;
					cout << "modulo della forza: " << scientific << fm << " (" << sqrt( SQR(fx[n]) + SQR(fy[n]) + SQR(fz[n]) ) << ") " << "\n" << endl;
		*/
        }
      }
    }
  }
}


// Calculating the position and velocity of each cell (it is calculated only if there are at least two cells)
void CellsSystem::NewPositionsAndVelocities( )
{
	
	// inizializzazione dei vettori
	vxnew = vx;
	vynew = vy;
	vznew = vz;
	
	
	loop_count = 0;	// conteggio dei passaggi nel loop
	
	if( ncells > 1 )
  {
		while( true )	// The infinite loop is interrupted by a break when the stop condition is reached
    {
			double dvmax = 0;	// max modulo della diff di velocita' tra due iterazioni

//#pragma omp parallel for ordered schedule(dynamic)
#pragma omp parallel for
			for(unsigned long n=0; n<ncells; n++)	// loop on the cells
      {
				int nneigh = neigh[n];		// number of adjacent cells
				double rn = r[n];					// cell radius (now used to calculate the friction coefficient)
				
				// inizializzazioni
				
				// coefficients of friction
				// WARNING ... this statement should enter the phenotype
				double gamma = 6.*PI*VISCOSITY_ENV*rn;
				double gamma_int = 0.1e15;	// viscosita' interna in pg/(micron s)
				
				// elementi di matrice: questi elementi vengono calcolati al volo nel loop, ma potrebbero venire valutati prima del loop una volta per tutte
				// in questo modo si utilizzerebbe memoria, ma il programma diventerebbe piu' veloce ... potrebbe essere un'ottimizzazione da fare in futuro
				double axx = 0;			
				double axy = 0.;
				double axz = 0.;
				double ayy = 0;
				double ayz = 0.;
				double azz = 0;
				
				// vettore "costante"
				double bx = 0.;
				double by = 0.;
				double bz = 0.;
				
				double mn = mass[n];				// massa della cellula
				
				for( int k=0; k< nneigh; k++)							// loop sulle cellule adiacenti
        {
					int neighbor = vneigh[n][k];		// nome della k-esima cellula adiacente
					
					double drx = (x[n]-x[neighbor]);				// componenti del raggio vettore tra le cellule n e neighbor
					double dry = (y[n]-y[neighbor]);
					double drz = (z[n]-z[neighbor]);
					double dd2 = drx*drx+dry*dry+drz*drz;			// distanza al quadrato tra cellula n-esima e la k-esima cellula vicina
										
					double sp = vx[neighbor]*drx + vy[neighbor]*dry + vz[neighbor]*drz;
					
					bx += sp*drx/dd2;
					by += sp*dry/dd2;
					bz += sp*drz/dd2;
					
					axx += drx*drx/dd2;
					axy += drx*dry/dd2;
					axz += drx*drz/dd2;
					ayy += dry*dry/dd2;
					ayz += dry*drz/dd2;
					azz += drz*drz/dd2;
        }//check, all used variables are local --> openMP compatible!
				
				bx = vx[n] + (gamma_int*bx + fx[n])*dt/mn;
				by = vy[n] + (gamma_int*by + fy[n])*dt/mn;
				bz = vz[n] + (gamma_int*bz + fz[n])*dt/mn;
				
				axx = 1. + (gamma + gamma_int*axx)*dt/mn;
				axy = gamma_int*axy*dt/mn;
				axz = gamma_int*axz*dt/mn;
				ayy = 1. + (gamma + gamma_int*ayy)*dt/mn;
				ayz = gamma_int*ayz*dt/mn;
				azz = 1. + (gamma + gamma_int*azz)*dt/mn;
				
				double det = axx*ayy*azz + 2.*axy*axz*ayz - axx*ayz*ayz - ayy*axz*axz - azz*axy*axy;
				
				vxnew[n] = ( (ayy*azz-ayz*ayz)*bx + (axz*ayz-axy*azz)*by + (axy*ayz-axz*ayy)*bz )/det;
				vynew[n] = ( (axz*ayz-axy*azz)*bx + (axx*azz-axz*axz)*by + (axy*axz-ayz*axx)*bz )/det;
				vznew[n] = ( (axy*ayz-axz*ayy)*bx + (axy*axz-ayz*axx)*by + (axx*ayy-axy*axy)*bz )/det;

        /* T.F. 
         * I am not sure whether this is threadsafe, what happens if
         * one thread reads vxnew[n] while an other thread writes vxnew[m] ?
         * I try and see what happens.
         */
				double dv = sqrt((vx[n]-vxnew[n])*(vx[n]-vxnew[n])+(vy[n]-vynew[n])*(vy[n]-vynew[n])+(vz[n]-vznew[n])*(vz[n]-vznew[n]));
				if(dv > dvmax) dvmax = dv;

// statement per il debugging, scommentare se necessario
/*				
				cout << "\n------------------------" << endl;
				cout << "cellula: " << n << endl;
				cout << "massa: " << scientific << mn << endl;
				cout << "coefficiente di attrito con il mezzo: " << scientific << 6.*PI*VISCOSITY_ENV*rn << endl;
				cout << "matrice M:" << endl;
				cout << axx << "\t " << axy << "\t " << axz << endl;
				cout << axy << "\t " << ayy << "\t " << ayz << endl;
				cout << axz << "\t " << ayz << "\t " << azz << endl;
				cout << "\nmatrice inversa: " << endl;
				cout << (ayy*azz-ayz*ayz)/det << "\t " << (axz*ayz-axy*azz)/det << "\t " << (axy*ayz-axz*ayy)/det << endl;
				cout << (axz*ayz-axy*azz)/det << "\t " << (axx*azz-axz*axz)/det << "\t " << (axy*axz-ayz*axx)/det << endl;
				cout << (axy*ayz-axz*ayy)/det << "\t " << (axy*axz-ayz*axx)/det << "\t " << (axx*ayy-axy*axy)/det << endl;
				cout << "\nvettore b" << endl;
				cout << bx << "\t " << by << "\t " << bz << "\n " << endl;
				cout << "v:    (" << scientific << vx[n] << ", " << vy[n] << ", " << vz[n] << ") " << endl;
				cout << "vnew: (" << scientific << vxnew[n] << ", " << vynew[n] << ", " << vznew[n] << ") " << endl;
*/				
      }// end // loop on the cells
			
			loop_count++; 
			
			//if(loop_count > 100) exit(0);
			// if( dvmax < (1.e-11) ) break;	// la condizione di interruzione del loop e' che l'imprecisione nella determinazione della velocità sia minore di 0.01 nm/s
			if( dvmax < delta_vmax ) break;	// la condizione di interruzione del loop e' che l'imprecisione nella determinazione della velocità sia minore di 0.1 nm/s
#pragma omp parallel for
			for(unsigned long n=0; n<ncells; n++)	// qui si copiano in v i nuovi valori vnew
      {
				vx[n] = vxnew[n];
				vy[n] = vynew[n];
				vz[n] = vznew[n];
      }
    }// end infinite loop
  }//if( ncells > 1 )
	else 
	// solution if there is only one cell: the only coefficient of friction is that with the environment and the only index is' n = 0
  {
		double rn = r[0];					// raggio dell'unica cellula (ora serve al calcolo del coeff di attrito)
    double mn = mass[0];				// massa della cellula
    double gamma = 6.*PI*VISCOSITY_ENV*rn;
		
		vx[0] = vx[0]/(1.+gamma*dt/mn);
		vy[0] = vy[0]/(1.+gamma*dt/mn);
		vz[0] = vz[0]/(1.+gamma*dt/mn);
  }


	// calcolo delle nuove posizioni
	
	maxdr = 0;	// inizializzazione dello spostamento massimo nel sistema di cellule
	
#pragma omp parallel for
	for(unsigned long n=0; n<ncells; n++)				// loop sulle cellule
  {
		double rn = r[n];			// raggio della cellula

		double dx = vx[n]*dt;						// variazioni delle coordinate
		double dy = vy[n]*dt;
		double dz = vz[n]*dt;
			
		double dr = sqrt( dx*dx + dy*dy + dz*dz);	// spostamento totale
		if(dr>maxdr) maxdr = dr;						// update dello spostamento massimo
			
		if( dr > rn )		// management of possible errors in the calculation of the position
    {
			std::cout << "Errore in Cells::NewPositionsAndVelocities: " << std::endl;
			std::cout << "Al passo " << nstep;
			std::cout << " si e' verificato un problema nel calcolo della posizione della cellula " << n << "-esima" << std::endl;
			std::cout << "dx = " << std::scientific << dx << std::endl;
			std::cout << "dy = " << std::scientific << dy << std::endl;
			std::cout << "dz = " << std::scientific << dz << std::endl;
			std::cout << "stepsize: " << std::scientific << dr << std::endl;
			std::cout << "cell radius: " << std::scientific << r[n] << std::endl;
			std::cout << "packing factor: " << (type[n]->Get_packing_factor()) << std::endl;
			std::cout << "cellForcex = " << std::scientific << fx[n] << std::endl;
			std::cout << "cellForcey = " << std::scientific << fy[n] << std::endl;
			std::cout << "cellForcez = " << std::scientific << fz[n] << std::endl;
			exit(0);
    }
    // calcolo delle nuove coordinate
		x[n] += dx;
		y[n] += dy;
		z[n] += dz;
  }
}


// dinamica dummy
void CellsSystem::DummyDynamics( )
{
	maxdr = 0;	// spostamento massimo
}

