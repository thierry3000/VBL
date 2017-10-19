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

#include "CellsSystem.h"
using namespace vbl;

extern double EventTime;
extern bool eventON;
extern double pAlt;

//
// A method that handles cellular events
//
// The method returns the number of mitoses that occurred during the step
// 
bool CellsSystem::CellEvents( )
{

	bool mitosis_flag = false;				// flag logica che dice se e' avvenuta almeno una mitosi
	n_mitoses = 0;							// numero di mitosi in questo passo

	unsigned long ncells_now = ncells;		// qui si immagazzina il numero iniziale di cellule

	alive = 0;								// si azzera il contatore delle cellule vive
	
	
	for(unsigned long n=0; n<ncells_now; n++)	// il loop va fino al numero iniziale di cellule				
		{
		
		// incremento dell'eta' cellulare
		age[n] += (float)dt;
		phase_age[n] += (float)dt;

		
	// ***** cambiamenti di fase *****
				
		double ConcAcL = AcL[n]/volume[n];
		
		// condizioni istantanee di morte cellulare (transizione irreversibile alla fase dead)
		if( (ready2start) && (phase[n] != dead) )
			{
			
			death_condition[n] = 0;											// indicatore dello stato di morte
			
			double dose = dose_rateSignal.SignalIntegral(treal-dt,treal);	// dose di radiazione integrata durante il passo
			
			//if( ATPp[n] < ATPmin[n] )
			//	death_condition[n] += 1;
			if( pHi[n] < type[n]->Get_pHimin() )
				death_condition[n] += 2;
			if( ran2(idum) < 1.-exp(-(type[n]->Get_a_R()) * ConcAcL*dt \
				  -(type[n]->Get_alpha_R( phase[n] ))*dose \
				  -(type[n]->Get_beta_R( phase[n] ))*pow(dose,2)) )
				death_condition[n] += 4;

			if( death_condition[n] != 0 )
				{
				phase[n] = dead;
				phase_age[n] = 0.;
				}
				
			}


	// qui si vede se ci sono cambiamenti di fase per cellule vive
	
	// transizione G1m-G1p
		if(phase[n] == G1m_phase && ConcS[n] < (type[n]->Get_Thresh_S_start())*(type[n]->Get_ConcS_0()) ) 		
			{
			
			// test di possibile morte cellulare all'inizio della fase G1p
			if( ATPp[n] < ATPmin[n] )
				{
				phase[n] = dead;
				death_condition[n] += 8;
				phase_age[n] = 0.;
				}
			// se la cellula e' vitale allora cambia fase
			else
				{
				phase[n] = G1p_phase;
				phase_age[n] = 0.;
				}
			}
	// transizione G1p-S
		else if(phase[n] == G1p_phase && ConcS[n] < (type[n]->Get_Thresh_S_stop())*(type[n]->Get_ConcS_0()) )	
			{
			// test di possibile morte cellulare all'inizio della fase S
			if( ATPp[n] < ATPmin[n] )
				{
				phase[n] = dead;
				death_condition[n] += 16;
				phase_age[n] = 0.;
				}
			// se la cellula e' vitale allora cambia fase
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

	// ***** inizio MITOSI! *****
		else if( phase[n] == M_phase && phase_age[n] > M_T[n] )								
			{
                
			// si segnala all'esterno che e' avvenuta una mitosi ... 
			mitosis_flag = true;
			n_mitoses++;
			
		// *** memorizzazione di alcuni valori importanti
			
			// eta' della madre
			age_mother[n] = age[n];
			
			// qui si memorizzano i parametri della cellula madre che servono al calcolo della suddivisione
			double old_volume = volume[n];
			double old_M = M[n];
			double old_r = r[n];											
			
			double ConcG_extra = G_extra[n]/volume_extra[n];		// concentrazioni nello spazio extracellulare
			double ConcA_extra = A_extra[n]/volume_extra[n];
			double ConcAcL_extra = AcL_extra[n]/volume_extra[n];

			double old_x = x[n];		// vecchie coordinate del centro
			double old_y = y[n];
			double old_z = z[n];
				
			
		// *** prima parte dell'update della cellula (variabili comuni alle due figlie)
			
			// reset della fase e dell'eta'
			phase[n] = G1m_phase;
			age[n] = 0.;
			phase_age[n] = 0.;
			
			// update del numero di mitosi
			n_mitosis[n]++;
			
			// reset delle cicline
			cyclinD[n] = 0.;
			cyclinE[n] = 0.;
			cyclinX[n] = 0.;
			
			// reset del DNA sintetizzato
			DNA[n] = 0.;			

			// reset di ConcS			
			ConcS[n] = type[n]->Get_ConcS_0();			


		// *** si copia la struttura della cellula madre nella nuova cellula (che va in fondo alla lista)
		// P.S. questa istruzione aumenta automaticamente il numero di instances di type
			if(treal < EventTime || (ran2(idum) > pAlt || pAlt < 2 ))
                ReplicateCell( n );
            else
                {
                ReplicateCell( n, &CellTypeVector[1] );
                if(pAlt == 2) pAlt = 0;
		}
		

		// a questo punto ncells e' stato incrementato di 1, e la cellula appena inserita si trova in posizione ncells-1
		
		
		// *** seconda parte dell'update della cellula (variabili diverse per le due figlie)
			
			// questa parte del programma serve ad modellizzare la stocasticita' prodotta dalla ripartizione ineguale dei mitocondri
			// clusters_M e' il numero di clusters in cui si riuniscono i mitocondri: ad esempio se ci sono 170 mitocondri, con un ClusteringFactor = 15 
			// (ClusteringFactor e' il numero medio di mitocondri per cluster) si ottiene clusters_M = 11, rest_M = 5; 

			int ClusteringFactor = type[n]->Get_ClusteringFactor();
			
			int clusters_M = (int)floor(old_M/ClusteringFactor);
			int rest_M = (int)(old_M - clusters_M * ClusteringFactor);
			
			// dopo aver definito i clusters, quello che segue suddivide i cluster in modo binomiale, quindi trova il numero di mitocondri che vanno nella cellula n-esima e 
			// quelli che vanno in newcell. Il resto e' assegnato casualmente. Ad esempio se l'if e' vero si assegnano i resti alla cellula n-esima, altrimenti no; in questo secondo 
			// caso l'assegnazione a newcell e' automatica, perche' newcell contiene tutto quello che la cellula n-esima non contiene
			
			double Mit = bnldev(0.5, clusters_M, idum) * ClusteringFactor;	// distribuzione binomiale
			if(ran2(idum) > 0.5) Mit += rest_M;
			
			M[n] = Mit;										// questi statements definiscono M e ATPmin
			M[ncells-1] = old_M - Mit;
			
			
			// calcoli relativi al volume
			double Vmin = type[n]->Get_Vmin();
			double C2 = type[n]->Get_C2();
			// qui si definisce il volume citoplasmatico sulla base della suddivisione del numero di mitocondri
			// si noti che al momento della suddivisione ci sono 2 nuclei e quindi il volume citoplasmatico totale è old_volume-2*Vmin
			double volume_C = (old_volume-2.*Vmin-C2*old_M) * Mit/old_M;
			double volume_C2 = (old_volume-2.*Vmin-C2*old_M) * (old_M-Mit)/old_M;
						
			// volumi totali e quantità associate (in questo calcolo si assume DNA = 0 visto che le cellule sono in fase G1m)
			volume[n] =  Vmin + C2*M[n] + volume_C;
			r[n] = pow(3.*volume[n]/(4.*PI), (double)1./3.); 
			surface[n] = 4.*PI*r[n]*r[n]; 
			mass[n] = type[n]->density * volume[n];
			volume_extra[n] = surface[n]*(type[n]->extvolume_thickness)*(type[n]->extvolume_fraction);
			ATPp[n] = (volume[n] - type[n]->C2 * M[n] - type[n]->Vmin)/type[n]->C1;
			ATPmin[n] = (type[n]->fATPmin)*(type[n]->C2 * M[n])/type[n]->C1;
			
			volume[ncells-1] = Vmin + C2*M[ncells-1] + volume_C2;
			r[ncells-1] = pow(3.*volume[ncells-1]/(4.*PI), (double)1./3.); 
			surface[ncells-1] = 4.*PI*r[ncells-1]*r[ncells-1]; 
			mass[ncells-1] = type[ncells-1]->density * volume[ncells-1];
			volume_extra[ncells-1] = surface[ncells-1]*(type[ncells-1]->extvolume_thickness)*(type[ncells-1]->extvolume_fraction);
			ATPp[ncells-1] = (volume[ncells-1] - type[ncells-1]->C2 * M[ncells-1] - type[ncells-1]->Vmin)/type[ncells-1]->C1;
			ATPmin[ncells-1] = (type[ncells-1]->fATPmin)*(type[ncells-1]->C2 * M[ncells-1])/type[ncells-1]->C1;
			
			// rapporto dei volumi
			double volume_ratio = volume[n]/old_volume;
			double volume_ratio2 = volume[ncells-1]/old_volume;
			
			// la concentrazione delle sostanze che diffondono negli spazi extracellulari e' la stessa di quella della cellula madre
			// (si noti che questo non conserva la quantita' di sostanza, in altre parole si assume che la diffusione sia sufficientemente
			// veloce da poter fare questa assunzione)
			G_extra[n] = ConcG_extra*volume_extra[n];
			A_extra[n] = ConcA_extra*volume_extra[n];
			AcL_extra[n] = ConcAcL_extra*volume_extra[n];
			
			G_extra[ncells-1] = ConcG_extra*volume_extra[ncells-1];
			A_extra[ncells-1] = ConcA_extra*volume_extra[ncells-1];
			AcL_extra[ncells-1] = ConcAcL_extra*volume_extra[ncells-1];
			
			// setup delle variabili metaboliche ereditate dalla cellula madre
			// (inizialmente i valori sono uguali nelle due cellule e quindi vengono presi solo da CellVector[n] )
			
			double Gnow = G[n];								// glucosio
			G[n] = Gnow*volume_ratio;
			G[ncells-1] = Gnow*volume_ratio2;
			
			double G6Pnow = G6P[n];							// G6P
			G6P[n] = G6Pnow*volume_ratio;
			G6P[ncells-1] = G6Pnow*volume_ratio2;
			
			double Anow = A[n];								// glutammina
			A[n] = Anow*volume_ratio;
			A[ncells-1] = Anow*volume_ratio2;
			
			double storenow = store[n];						// store
			store[n] = storenow*volume_ratio;
			store[ncells-1] = storenow*volume_ratio2;
			
															
			ATPstart[n] = (double)ATPp[n];									// inizializzazione dei contatori dell'ATP
			ATPstart[ncells-1] = (double)ATPp[ncells-1];
			ATPprod[n] = 0.;
			ATPprod[ncells-1] = 0.;
			ATPcons[n] = 0.;
			ATPcons[ncells-1] = 0.;
			
			double AcLnow = AcL[n];							// AcL
			AcL[n] = AcLnow*volume_ratio;
			AcL[ncells-1] = AcLnow*volume_ratio2;
			
			// double H = CellVector[n].Get_H();								// H
			// CellVector[n].Set_H( H*volume_ratio );
			// newcell.Set_H( H*volume_ratio2 );
			
			double O2now = O2[n];
			O2[n] = O2now*volume_ratio;
			O2[ncells-1] = O2now*volume_ratio2;
			
			// double CO2 = CellVector[n].Get_CO2();
			// CellVector[n].Set_CO2( CO2*volume_ratio );
			// newcell.Set_CO2( CO2*volume_ratio2 );

			// si assume le proteine siano diffuse in tutta la cellula, compreso il nucleo, e quindi lo splitting come per le altre sostanze
			double proteinnow = protein[n];
			protein[n] = proteinnow*volume_ratio;
			protein[ncells-1] = proteinnow*volume_ratio2;

			// l'algoritmo di splitting della pRb e' basato sull'idea che si appiccichi alla matrice nucleare
			double splitting_fraction = (double) bnldev(0.5, type[n]->Get_NUCLEAR_OBJ(), idum); // distribuzione binomiale
			double splitting_ratio = splitting_fraction / type[n]->Get_NUCLEAR_OBJ();
			double splitting_ratio2 = 1.-splitting_ratio;
			
			double pRbnow = pRb[n];
			pRb[n] = pRbnow*splitting_ratio;
			pRb[ncells-1] = pRbnow*splitting_ratio2;


			// definizione di DNA_spread e M_T
			DNA_spread[n] = type[n]->Get_DNA_MAX_SPREAD() * (2.*ran2(idum)-1.);
			M_T[n] = type[n]->Get_M_T_MEAN() * (1.+ type[n]->Get_PHASE_SPREAD() * (2.*ran2(idum)-1.));

			DNA_spread[ncells-1] = type[ncells-1]->Get_DNA_MAX_SPREAD() * (2.*ran2(idum)-1.);
			M_T[ncells-1] = type[ncells-1]->Get_M_T_MEAN() * (1.+ type[ncells-1]->Get_PHASE_SPREAD() * (2.*ran2(idum)-1.));
			

			// le cellule sono temporanemente create isolate
			isonCH[n] = true;
			isonAS[n] = true;
			neigh[n] = 0 ;
			contact_surf[n] = 0.;
			env_surf[n] = surface[n];
			g_env[n] = env_surf[n]/r[n];
			
			isonCH[ncells-1] = true;
			isonAS[ncells-1] = true;
			neigh[ncells-1] = 0 ;
			contact_surf[ncells-1] = 0.;
			env_surf[ncells-1] = surface[ncells-1];
			g_env[ncells-1] = env_surf[ncells-1]/r[ncells-1];
						

		// *** se la simulazione e' pronta allora si calcolano le nuove coordinate
			if( ready2start )
				{
				
				// geometria della mitosi (solo nel caso di simulazione Full3D)
				if( sim_type == Full3D ) 
					{
					double xr = 1.-2.*ran2(idum); 
					double yr = 1.-2.*ran2(idum); 
					double zr = 1.-2.*ran2(idum); 
					double len = sqrt(xr*xr + yr*yr + zr*zr);
					xr /= len;
					yr /= len;
					zr /= len;	// a questo punto {xr, yr, zr} e' un versore casuale

					double r_1 = r[n];							// i nuovi valori del raggio che vengono dalla routine del metabolismo
					double r_2 = r[ncells-1];
					
					double x_1 = old_x + xr*(old_r - r_1);		// calcolo delle nuove coordinate dei centri delle cellule
					double y_1 = old_y + yr*(old_r - r_1);
					double z_1 = old_z + zr*(old_r - r_1);
					
					double x_2 = old_x - xr*(old_r - r_2);
					double y_2 = old_y - yr*(old_r - r_2);
					double z_2 = old_z - zr*(old_r - r_2);
					
					x[n] = x_1;										// qui si immagazzinano le posizioni dei centri
					y[n] = y_1;
					z[n] = z_1;

					x[ncells-1] = x_2;
					y[ncells-1] = y_2;
					z[ncells-1] = z_2;
					
					vx[ncells-1] = vx[n];	// la velocita' del centro della nuova cellula è inizialmente uguale a quella della cellula madre
					vy[ncells-1] = vy[n];
					vz[ncells-1] = vz[n];
					}
				}


			// nel caso di startup o conf. fissa, si dimentica l'ultima cellula ... altrimenti si aggiorna il contatore delle cellule vive			
			if( !ready2start || sim_type == FixedConfig )
				{
				ncells--;						// decremento del numero di cellule, si torna al valore precedente, l'ultima cellula viene dimenticata
				type[n]->Delete_instance();		// cancellazione dell'instance di type che era stato aggiunto da ReplicateCell
				}
			else 
				alive++;

			
			} // fine del trattamento della mitosi

		
		// controlli di consistenza di alcune variabili
		if(phase[n] != dead) 
			{
			int code = CheckMVA( n );
			if(code < 0) errorlog_file << "Errore " << code << " alla fine di CellsSystem::CellEvents nel controllo di consistenza per la cellula " << n << "\n" << std::endl;
			}

	
		if(phase[n] != dead) alive++;								// se la cellula non e' morta, allora e' viva ...	

		
		} // fine del loop sulle cellule
	

	return mitosis_flag;

}



//
// metodo che elimina le cellule morte troppo piccole
// 
void CellsSystem::CleanCellsSystem( )
{
	//bool cleaned = false;
	
	for(unsigned long n=ncells; n>1; n--)				
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
