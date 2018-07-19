/*
 *  CellsSystem-C.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  2 methods that take care of cellular events
 *  CellEvents and CleanCellsSystem
 *
 */


//
// A method that handles cellular events
//
// The method returns the number of mitoses that occurred during the step
// 
//bool CellsSystem::CellEvents( )
int CellsSystem::CellEvents( )
{

	bool mitosis_flag = false;				// logical flag that says if at least one mitosis has occurred
	n_mitoses = 0;							// number of mitosis in this passage

	unsigned long ncells_now = params.ncells;		// here the initial number of cells is stored

	params.alive = 0;								// the live cell counter is zeroed
	
	
	for(unsigned long n=0; n<ncells_now; n++)	// the loop goes up to the initial number of cells
  {
		// increase in cellular age
		age[n] += (float)params.dt;
		phase_age[n] += (float)params.dt;

		/*
     * **** PHASE CHANGES *****
     */
				
		double ConcAcL = AcL[n]/volume[n];
		
		// instantaneous conditions of cell death (irreversible transition to the dead phase)
		if( (ready2start) && (phase[n] != dead) )
    {
			death_condition[n] = 0;											// indicator of the state of death
			
			double dose = dose_rateSignal.SignalIntegral(params.treal-params.dt,params.treal);	// radiation dose integrated during the step
			
			//if( ATPp[n] < ATPmin[n] )
			//	death_condition[n] += 1;
			if( pHi[n] < type[n]->Get_pHimin() )
				death_condition[n] += 2;
			if( ran2() < 1.-exp(-(type[n]->Get_a_R()) * ConcAcL*params.dt \
				  -(type[n]->Get_alpha_R( phase[n] ))*dose \
				  -(type[n]->Get_beta_R( phase[n] ))*pow(dose,2)) )
				death_condition[n] += 4;

			if( death_condition[n] != 0 )
      {
				phase[n] = dead;
				phase_age[n] = 0.;
			}
				
    }


	// here we see if there are phase changes for living cells
	
	// transizione G1m-G1p
		if(phase[n] == G1m_phase && ConcS[n] < (type[n]->Get_Thresh_S_start())*(type[n]->Get_ConcS_0()) ) 		
    {
			// test of possible cell death at the beginning of the G1p phase
			if( ATPp[n] < ATPmin[n] )
      {
				phase[n] = dead;
				death_condition[n] += 8;
				phase_age[n] = 0.;
			}
			// if the cell is vital then it changes phase
			else
      {
				phase[n] = G1p_phase;
				phase_age[n] = 0.;
      }
    }
	// transizione G1p-S
		else if(phase[n] == G1p_phase && ConcS[n] < (type[n]->Get_Thresh_S_stop())*(type[n]->Get_ConcS_0()) )	
    {
			// test of possible cell death at the beginning of the S phase
			if( ATPp[n] < ATPmin[n] )
      {
				phase[n] = dead;
				death_condition[n] += 16;
				phase_age[n] = 0.;
			}
			// if the cell is vital then it changes phase
			else
      {
				phase[n] = S_phase;
				phase_age[n] = 0.;
				cyclinD[n] = 0.;
				cyclinE[n] = 0.;
      }
    }
	// transizione S-G2
		else if(phase[n] == S_phase && DNA[n] >= 1 )											
    {
			// test di possibile morte cellulare all'inizio della fase G2
			if( ATPp[n] < ATPmin[n] )
      {
				phase[n] = dead;
				death_condition[n] += 32;
				phase_age[n] = 0.;
      }
			// se la cellula e' vitale allora cambia fase
			else
      {
				phase[n] = G2_phase;
				phase_age[n] = 0.;
      }
    }
	// transizione G2-M
		else if(phase[n] == G2_phase && cyclinX[n] > type[n]->Get_CycXThr() )						
    {
			// test di possibile morte cellulare all'inizio della fase M
			if( ATPp[n] < ATPmin[n] )
      {
				phase[n] = dead;
				death_condition[n] += 64;
				phase_age[n] = 0.;
      }
			// se la cellula e' vitale allora cambia fase
			else
      {
				phase[n] = M_phase;
				phase_age[n] = 0.;
      }
    }

	// ***** MITOSIS start! *****
		else if( phase[n] == M_phase && phase_age[n] > M_T[n] )								
    {
			// it is pointed out to the outside that there has been a mitosis
			mitosis_flag = true;
			n_mitoses++;
			
		// *** memorization of some important values
			
			// age of the mother
			age_mother[n] = age[n];
			
			// here the parameters of the mother cell that are used for the calculation of the subdivision are memorized
			double old_volume = volume[n];
			double old_M = M[n];
			double old_r = r[n];											
			
			double ConcG_extra = G_extra[n]/volume_extra[n];		// concentrations in the extracellular space
			double ConcA_extra = A_extra[n]/volume_extra[n];
			double ConcAcL_extra = AcL_extra[n]/volume_extra[n];

			double old_x = x[n];		// old center coordinates
			double old_y = y[n];
			double old_z = z[n];
				
			
		// *** first part of the cell's update (variables common to the two daughters)
			
			// phase and age reset
			phase[n] = G1m_phase;
			age[n] = 0.;
			phase_age[n] = 0.;
			
			// update of the number of mitosis
			n_mitosis[n]++;
			
			// reset delle cicline
			cyclinD[n] = 0.;
			cyclinE[n] = 0.;
			cyclinX[n] = 0.;
			
			// reset del DNA sintetizzato
			DNA[n] = 0.;			

			// reset di ConcS			
			ConcS[n] = type[n]->Get_ConcS_0();			


		// *** copy the structure of the mother cell into the new cell (which goes to the bottom of the list)
		// P.S. this instruction automatically increases the number of type instances
			if(params.treal < EventTime || (ran2() > pAlt || pAlt < 2 ))
      {
        ReplicateCell( n );
      }
      else
      {
        ReplicateCell( n, &CellTypeVector[1] );
        if(pAlt == 2) 
          pAlt = 0;
      }
		

      // at this point ncells has been increased by 1, and the newly inserted cell is in the ncells-1 position
		
		
      // *** second part of the cell update (different variables for the two daughters)
			
			// this part of the program serves to model the stochasticity produced by the unequal distribution of mitochondria
			// clusters_M is the number of clusters in which mitochondria gather: for example, if there are 170 mitochondria, with a ClusteringFactor = 15 
			// (ClusteringFactor is the average number of mitochondria per cluster) we obtain clusters_M = 11, rest_M = 5; 

			int ClusteringFactor = type[n]->Get_ClusteringFactor();
			
			int clusters_M = (int)floor(old_M/ClusteringFactor);
			int rest_M = (int)(old_M - clusters_M * ClusteringFactor);
			
			/**
       * after having defined the clusters, the following one divides the clusters in a binomial way,
       * then find the number of mitochondria going into the n-th cell and those going into newcell.
       * The rest is assigned randomly. For example, if the if is true, the remains are assigned to the n-th cell,
       * otherwise no; in this second case the assignment to newcell is automatic,
       * because newcell contains everything that the n-th cell does not contain
       */ 
			
			double Mit = bnldev(0.5, clusters_M, params.idum) * ClusteringFactor;	// distribuzione binomiale
			if(ran2() > 0.5) Mit += rest_M;
			
			M[n] = Mit;										// questi statements definiscono M e ATPmin
			M[params.ncells-1] = old_M - Mit;
			
			
			// calculations related to the volume
			double Vmin = type[n]->Get_Vmin();
			double C2 = type[n]->Get_C2();
			// here the cytoplasmic volume is defined on the basis of the subdivision of the number of mitochondria
      // note that at the time of the subdivision there are 2 nuclei and therefore the total cytoplasmic volume is old_volume-2 * Vmin
			double volume_C = (old_volume-2.*Vmin-C2*old_M) * Mit/old_M;
			double volume_C2 = (old_volume-2.*Vmin-C2*old_M) * (old_M-Mit)/old_M;
						
			// total volumes and associated quantities (in this calculation we assume DNA = 0 since the cells are in phase G1m)
			volume[n] =  Vmin + C2*M[n] + volume_C;
			r[n] = pow(3.*volume[n]/(4.*PI), (double)1./3.); 
			surface[n] = 4.*PI*r[n]*r[n]; 
			mass[n] = type[n]->density * volume[n];
			volume_extra[n] = surface[n]*(type[n]->extvolume_thickness)*(type[n]->extvolume_fraction);
			ATPp[n] = (volume[n] - type[n]->C2 * M[n] - type[n]->Vmin)/type[n]->C1;
			ATPmin[n] = (type[n]->fATPmin)*(type[n]->C2 * M[n])/type[n]->C1;
			
			volume[params.ncells-1] = Vmin + C2*M[params.ncells-1] + volume_C2;
			r[params.ncells-1] = pow(3.*volume[params.ncells-1]/(4.*PI), (double)1./3.); 
			surface[params.ncells-1] = 4.*PI*r[params.ncells-1]*r[params.ncells-1]; 
			mass[params.ncells-1] = type[params.ncells-1]->density * volume[params.ncells-1];
			volume_extra[params.ncells-1] = surface[params.ncells-1]*(type[params.ncells-1]->extvolume_thickness)*(type[params.ncells-1]->extvolume_fraction);
			ATPp[params.ncells-1] = (volume[params.ncells-1] - type[params.ncells-1]->C2 * M[params.ncells-1] - type[params.ncells-1]->Vmin)/type[params.ncells-1]->C1;
			ATPmin[params.ncells-1] = (type[params.ncells-1]->fATPmin)*(type[params.ncells-1]->C2 * M[params.ncells-1])/type[params.ncells-1]->C1;
			
			// ratio of volumes
			double volume_ratio = volume[n]/old_volume;
			double volume_ratio2 = volume[params.ncells-1]/old_volume;
			
			// the concentration of the substances that spread in the extracellular spaces is the same as that of the mother cell
			// (note that this does not preserve the amount of substance, in other words it is assumed that the diffusion is sufficient
      // fast to make this assumption)
			G_extra[n] = ConcG_extra*volume_extra[n];
			A_extra[n] = ConcA_extra*volume_extra[n];
			AcL_extra[n] = ConcAcL_extra*volume_extra[n];
			
			G_extra[params.ncells-1] = ConcG_extra*volume_extra[params.ncells-1];
			A_extra[params.ncells-1] = ConcA_extra*volume_extra[params.ncells-1];
			AcL_extra[params.ncells-1] = ConcAcL_extra*volume_extra[params.ncells-1];
			
			// setup of the metabolic variables inherited from the mother cell
      // (initially the values ​​are equal in the two cells and therefore are taken only by CellVector [n])
			
			double Gnow = G[n];								// glucosio
			G[n] = Gnow*volume_ratio;
			G[params.ncells-1] = Gnow*volume_ratio2;
			
			double G6Pnow = G6P[n];							// G6P
			G6P[n] = G6Pnow*volume_ratio;
			G6P[params.ncells-1] = G6Pnow*volume_ratio2;
			
			double Anow = A[n];								// glutammina
			A[n] = Anow*volume_ratio;
			A[params.ncells-1] = Anow*volume_ratio2;
			
			double storenow = store[n];						// store
			store[n] = storenow*volume_ratio;
			store[params.ncells-1] = storenow*volume_ratio2;
			
															
			ATPstart[n] = (double)ATPp[n];									// inizializzazione dei contatori dell'ATP
			ATPstart[params.ncells-1] = (double)ATPp[params.ncells-1];
			ATPprod[n] = 0.;
			ATPprod[params.ncells-1] = 0.;
			ATPcons[n] = 0.;
			ATPcons[params.ncells-1] = 0.;
			
			double AcLnow = AcL[n];							// AcL
			AcL[n] = AcLnow*volume_ratio;
			AcL[params.ncells-1] = AcLnow*volume_ratio2;
			
			// double H = CellVector[n].Get_H();								// H
			// CellVector[n].Set_H( H*volume_ratio );
			// newcell.Set_H( H*volume_ratio2 );
			
			double O2now = O2[n];
			O2[n] = O2now*volume_ratio;
			O2[params.ncells-1] = O2now*volume_ratio2;
			
			// double CO2 = CellVector[n].Get_CO2();
			// CellVector[n].Set_CO2( CO2*volume_ratio );
			// newcell.Set_CO2( CO2*volume_ratio2 );

			// assume proteins are widespread throughout the cell, including the nucleus, and therefore splitting as for other substances
			double proteinnow = protein[n];
			protein[n] = proteinnow*volume_ratio;
			protein[params.ncells-1] = proteinnow*volume_ratio2;

			// the splitting algorithm of the pRb is based on the idea that it sticks to the nuclear matrix
			double splitting_fraction = (double) bnldev(0.5, type[n]->Get_NUCLEAR_OBJ(), params.idum); // distribuzione binomiale
			double splitting_ratio = splitting_fraction / type[n]->Get_NUCLEAR_OBJ();
			double splitting_ratio2 = 1.-splitting_ratio;
			
			double pRbnow = pRb[n];
			pRb[n] = pRbnow*splitting_ratio;
			pRb[params.ncells-1] = pRbnow*splitting_ratio2;


			// definizione di DNA_spread e M_T
			DNA_spread[n] = type[n]->Get_DNA_MAX_SPREAD() * (2.*ran2()-1.);
			M_T[n] = type[n]->Get_M_T_MEAN() * (1.+ type[n]->Get_PHASE_SPREAD() * (2.*ran2()-1.));

			DNA_spread[params.ncells-1] = type[params.ncells-1]->Get_DNA_MAX_SPREAD() * (2.*ran2()-1.);
			M_T[params.ncells-1] = type[params.ncells-1]->Get_M_T_MEAN() * (1.+ type[params.ncells-1]->Get_PHASE_SPREAD() * (2.*ran2()-1.));
			

			// cells are temporally created isolated
			isonCH[n] = true;
			isonAS[n] = true;
			neigh[n] = 0 ;
			contact_surf[n] = 0.;
			env_surf[n] = surface[n];
			g_env[n] = env_surf[n]/r[n];
			
			isonCH[params.ncells-1] = true;
			isonAS[params.ncells-1] = true;
			neigh[params.ncells-1] = 0 ;
			contact_surf[params.ncells-1] = 0.;
			env_surf[params.ncells-1] = surface[params.ncells-1];
			g_env[params.ncells-1] = env_surf[params.ncells-1]/r[params.ncells-1];
						

		// *** if the simulation is ready then the new coordinates are calculated
			if( ready2start )
      {
				// geometria della mitosi (solo nel caso di simulazione Full3D)
				if( params.sim_type == Full3D ) 
        {
          //number in [-1 1]
					double xr = 1.-2.*ran2(); 
					double yr = 1.-2.*ran2(); 
					double zr = 1.-2.*ran2(); 
					double len = sqrt(xr*xr + yr*yr + zr*zr);
					xr /= len;
					yr /= len;
					zr /= len;	// at this point {xr, yr, zr} is a random vector of the sphere

					double r_1 = r[n];							// the new values ​​of the ray coming from the metabolism routine
					double r_2 = r[params.ncells-1];
					
          double moving_factor_1 = old_r -r_1;
					double x_1 = old_x + xr*(moving_factor_1);		// calculation of the new coordinates of the cell centers
					double y_1 = old_y + yr*(moving_factor_1);
					double z_1 = old_z + zr*(moving_factor_1);
					
          double moving_factor_2 = old_r -r_2;
					double x_2 = old_x - xr*(moving_factor_2);
					double y_2 = old_y - yr*(moving_factor_2);
					double z_2 = old_z - zr*(moving_factor_2);
					
					x[n] = x_1;										//Here the positions of the centers are stored
					y[n] = y_1;
					z[n] = z_1;

					x[params.ncells-1] = x_2;
					y[params.ncells-1] = y_2;
					z[params.ncells-1] = z_2;
					
					vx[params.ncells-1] = vx[n];	// the velocity of the center of the new cell is initially the same as that of the mother cell
					vy[params.ncells-1] = vy[n];
					vz[params.ncells-1] = vz[n];
					}
				}


			// in the case of startup or conf. fixed, we forget the last cell ... otherwise the live cell counter is updated
			if( !ready2start || params.sim_type == FixedConfig )
      {
				params.ncells--;						// decremento del numero di cellule, si torna al valore precedente, l'ultima cellula viene dimenticata
				type[n]->Delete_instance();		// cancellazione dell'instance di type che era stato aggiunto da ReplicateCell
			}
			else 
				params.alive++;
    } // fine del trattamento della mitosi

/* T.F.
 * comment out because of multithread errros on snowden
 */
// 		// controlli di consistenza di alcune variabili
		if(phase[n] != dead) 
    {
			int code = CheckMVA( n );
			if(code < 0) 
        errorlog_file << "Errore " << code << " alla fine di CellsSystem::CellEvents nel controllo di consistenza per la cellula " << n << "\n" << std::endl;
		}

	
		if(phase[n] != dead) 
      params.alive++;								// if the cell is not dead, then it's alive ...
  } // fine del loop sulle cellule
	//return mitosis_flag;
  return n_mitoses;
}



//
// metodo che elimina le cellule morte troppo piccole
// 
void CellsSystem::CleanCellsSystem( )
{
	//bool cleaned = false;
	
	for(unsigned long n=params.ncells; n>1; n--)				
		{
		if( (phase[n-1] == dead) && (volume[n-1] < volume_extra[n-1]) )
				{
				// if(!cleaned) cout << "CleanCellsSystem: " << ncells << " cells at start of call" << endl;
				RemoveCell(n-1);
				//cleaned = true;
				}
		}
		
	//if(cleaned) cout << "CleanCellsSystem: " << ncells << " cells at end of call" << endl;
}
