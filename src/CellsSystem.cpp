/*
 *  Cells.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 18/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */


#include "../include/vbl.h"


using namespace vbl;

#include <chrono>
#include <boost/format.hpp>


/** 
 * STRUCTURE of main loop:
 * 
 * A) Dynamics
 * B) Diff
 * C) CellEvents
 * D) Geometry
 * 
 * W_timing:	if defined, we store the runtimes of the functions
 * 		NOTE that this is machine dependent
 * 
 * while:	the loop runs as long as there are cells (Get_alive()>0)
 * 		and other ... criterias are met
 */
unsigned int CellsSystem::runMainLoop(boost::optional<double> endtime)
{
  bool active_run = true;	// Boolean variable that becomes false when a condition stops running
  int nthreads, tid;
#ifdef W_timing
  myTiming.reset();
#endif
  
  while(active_run)
  {
#ifdef W_timing
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#endif

    // The calculation of geometry and dynamic is only done if 3D calculation is selected and the system is ready to start
    if(Get_ready2start() && Get_sim_type() == Full3D )	
    {
      /** Calculation of mechanical interactions
       * calls GetForces();
       */
#ifdef ENABLE_RELEASE_DEBUG_OUTPUT
      checkNeighbourhood_consistency(std::string("before dynamics"));
#endif
      Dynamics( );
    // CellsSystem.Print2logfile("Cellule dopo una chiamata a Dynamics");				
    }
#ifdef W_timing
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    myTiming.time_diff = end-begin;
    myTiming.dynamics = myTiming.dynamics + myTiming.time_diff.count();
#endif
      

    CPU_timer(Start_intertime);			// start intertempo del CPU timer
#ifdef W_timing
    begin = std::chrono::steady_clock::now();
#endif
#ifdef ENABLE_RELEASE_DEBUG_OUTPUT
    checkNeighbourhood_consistency(std::string("before Diff"));
#endif
    
    Diff();								// metabolismo e diffusione
    
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff = end-begin;
    myTiming.diff = myTiming.diff + myTiming.time_diff.count();
#endif
    CPU_timer(Stop_intertime);			// stop intertempo del CPU timer


#ifdef W_timing
    begin = std::chrono::steady_clock::now();
#endif
    //bool mitosis = CellEvents( );		// Cellular events
    checkNeighbourhood_consistency(std::string("before CellEvents"));
    int no_mitosis = CellEvents( );		// Cellular events
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff = end-begin;
    myTiming.cellEvents = myTiming.cellEvents + myTiming.time_diff.count();
#endif
    
    // The calculation of geometry and dynamic is only done if 3D calculation is selected and if the system is ready to go
#ifdef W_timing
    begin = std::chrono::steady_clock::now();
#endif
    if(Get_ready2start() && Get_sim_type() == Full3D )	
    {
      /**
       * The calculation of the triangulation is done if there has been at least one mitosis or if the system has moved significantly
       * (maximum cell displacement is 0.5 micron) or, however, if at least 250 seconds have elapsed
       * since the last call to CGAL
       * The triangulation calculation is done every time if the timestep is bigger than 10 s
       */
      if( no_mitosis>0 || (Get_maxdr() > 0.5e-6) || Get_time_from_CGAL() > 250. || Get_dt() > 10.) 
      {
        //the geometry part assumes that there are at least 5 cells!
        if(Get_ncells() > 10)
        {
//           if(Get_ncells() > 5)
//           {
//             CleanCellsSystem( );		// First you do the cleaning of the memory (elimination of dead cells now too small)
//             Geometry( );				// Calculation of cluster geometry
//           }
//           else
//           {
//             CleanCellsSystem( );		// First you do the cleaning of the memory (elimination of dead cells now too small)
//             Geometry_serial( );				// Calculation of cluster geometry
//           }
#ifdef ENABLE_RELEASE_DEBUG_OUTPUT
	  checkNeighbourhood_consistency(std::string("in Geometry"));
#endif
          CleanCellsSystem();
	  try
	  {
	    Geometry();
	  }
	  catch(int e)
	  {
	    std::cout << " error: " << e << std::endl;
	    std::cout << " in geometry " << std::endl;
	  }
	  catch (const std::exception& e) 
	  { // caught by reference to base
	    std::cout << " a standard exception was caught, with message: \n" << e.what() << "\n";
	    exit(EXIT_FAILURE);
	  }
          
          Set_time_from_CGAL(0.);		// Timer reset from last call to CGAL
        }
        else 
        {
          NoGeometry( );
        }
      // CellsSystem.Print2logfile("Cellule dopo una chiamata a Geometry");
      }// end mitosis

      //CellsSystem.Dynamics( );						// calcolo delle interazioni meccaniche
      // CellsSystem.Print2logfile("Cellule dopo una chiamata a Dynamics");				
    }
    // if the cells are not yet ready or the case of disrupted cells is being considered
    // only a limited number of geometric properties are calculated
    else 
    {
      NoGeometry( );
    }
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff = end-begin;
    myTiming.geometry = myTiming.geometry + myTiming.time_diff.count();
#endif

    StepStat( false );// calcolo delle statistiche passo-passo
	      
    // if( CellsSystem.Get_ready2start() ) CellsSystem.PrintLog(0);


    //
    // *** Printing summary statistics every nscreen step, or at each step in the case of slowmotion ***
    //
#ifdef W_timing
    begin = std::chrono::steady_clock::now();
#endif
    
    if(Get_nstep()%Get_nscreen() == 0 || Get_slow_motion() )		
    // if(CellsSystem.Get_nstep()%CellsSystem.Get_nscreen() == 0 || mitosis)		
    {
      CPU_timer(Restart_timer);		// Updating the CPU timer
      Timing(false);					// Updating the real time timer
      
      Printout( );
      Print2logfile("Cell status during the loop");		// Comment on normal operation
      //if( CellsSystem.Get_ready2start() )
      //	CellsSystem.PrintAll2logfile("Stato delle cellule durante il loop");	// COMMENTARE PER IL FUNZIONAMENTO NORMALE
																			      // PRODUCE UN OUTPUT ENORME ... !!!
      CPU_timer(Clear_intertime);		// Interrupt timer reset

      if( Get_ready2start() ) Print2file();
														      
      StepStat( true );				// Reinitializing statistics (also resets the convergence_fail vector)
    }
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff = end-begin;
    myTiming.writeToFile = myTiming.writeToFile + myTiming.time_diff.count();
#endif


//     //
//     // *** Printing the configurations every apoint steps, or at each step in the case of slowmotion ***
//     //
//     if(Get_nstep()%Get_nprint() == 0 || Get_slow_motion() )		
//     {
//       if(Get_ready2start() && Get_sim_type() == Full3D )
//       {
//         #pragma omp parallel sections
//         {
//           #pragma omp section
//           {
//             PrintConfiguration(true);		// Printed on a binary file
//             //CellsSystem.PrintConfiguration(false);	// printout su file ascii
//           }// end #pragma omp section
// 
//           #pragma omp section
//           {
//             if(Get_ncells()>2)
//               /* calculation of flow only makes sense for more than one cell!!!*/
//               PrintFlows();					// printout dei flussi extracellulari su file binario
//           }//end #pragma omp section
//         }// end #pragma omp parallel sections
//         
//         Step_nconfiguration();	// print extracellular streams on a binary file
//       }// end CellsSystem.Get_ready2start()
//     }// end CellsSystem.Get_nstep()
		      
	      
    active_run = TimersAdvanceUntil( endtime );		// advancing timers (including the last call timer to CGAL)

    if(Get_alive() == 0) 
    {
      std::cout << "\nRun break, there are no more live cells" << std::endl;
      return 1;
    }
#ifdef ENABLE_RELEASE_DEBUG_OUTPUT
    checkNeighbourhood_consistency(std::string("before loop end"));
#endif

  }// end while
// #ifdef W_timing
//   timing_file.close();
// #endif
  return 0;
}


// Allocation / redeployment of the dynamic reserve
void CellsSystem::Set_reserve(const int reserve)
{

	name.reserve(reserve);
	mark.reserve(reserve);
	type.reserve(reserve);
	
	Temperature.reserve(reserve);
	
	phase.reserve(reserve);
	
	death_condition.reserve(reserve);
	age.reserve(reserve);
	phase_age.reserve(reserve);
	age_mother.reserve(reserve);
	n_mitosis.reserve(reserve);
	
	x.reserve(reserve);
	y.reserve(reserve);
	z.reserve(reserve);
	
	vx.reserve(reserve);
	vy.reserve(reserve);
	vz.reserve(reserve);
	
	vxnew.reserve(reserve);
	vynew.reserve(reserve);
	vznew.reserve(reserve);
	
	fx.reserve(reserve);
	fy.reserve(reserve);
	fz.reserve(reserve);
	
	v.reserve(reserve);
	
	r.reserve(reserve);
	surface.reserve(reserve);
	volume.reserve(reserve);
	mass.reserve(reserve);
	
	volume_extra.reserve(reserve);
	
	neigh.reserve(reserve);
	vneigh.reserve(reserve);
	vdist.reserve(reserve);
	vcsurf.reserve(reserve);
	gnk.reserve(reserve);
	contact_surf.reserve(reserve);
	
	isonCH.reserve(reserve);
	isonAS.reserve(reserve);
	isonBV.reserve(reserve);
	env_surf.reserve(reserve);
	bv_surf.reserve(reserve);
	g_env.reserve(reserve);
	g_bv.reserve(reserve);
	
	M.reserve(reserve);
	
	G.reserve(reserve);
	G6P.reserve(reserve);
	O2.reserve(reserve);
	store.reserve(reserve);
	A.reserve(reserve);
	AcL.reserve(reserve);
	
	h.reserve(reserve);
	pHi.reserve(reserve);
	
	protein.reserve(reserve);
	prot_rate.reserve(reserve);
	DNA.reserve(reserve);
	DNA_rate.reserve(reserve);
	
	GAbsRate.reserve(reserve);
	GConsRate.reserve(reserve);
	AAbsRate.reserve(reserve);
	AConsRate.reserve(reserve);
	StoreFillRate.reserve(reserve);
	StoreConsRate.reserve(reserve);
	AcLRate.reserve(reserve);
	AcLOutRate.reserve(reserve);
	O2Rate.reserve(reserve); //*************************** new for O2 rate ******* april 2018
	
	ATP_St.reserve(reserve);
	ATP_Ox.reserve(reserve);
	ATP_NOx.reserve(reserve);
	ATP2.reserve(reserve);
	ATP3.reserve(reserve);
	ConsATP.reserve(reserve);
	ConsATP_1.reserve(reserve);
	ConsATP_2.reserve(reserve);
	ConsATP_3.reserve(reserve);
	ConsATP_4.reserve(reserve);
	ConsATP_5.reserve(reserve);
	ATPtot.reserve(reserve);
	ATPp.reserve(reserve);
	ATPmin.reserve(reserve);
	
	ATPstart.reserve(reserve);
	ATPprod.reserve(reserve);
	ATPcons.reserve(reserve);
	
	G_extra.reserve(reserve);
	A_extra.reserve(reserve);
	AcL_extra.reserve(reserve);
	
	pH.reserve(reserve);
	SensO2.reserve(reserve);
	ConsO.reserve(reserve);
	
	DNA_spread.reserve(reserve);
	
	M_T.reserve(reserve);
	pRb.reserve(reserve);
	
	ConcS.reserve(reserve);
	
	cyclinD.reserve(reserve);
	cyclinE.reserve(reserve);
	cyclinX.reserve(reserve);
	
	NpRbk.reserve(reserve);
	
	volumeOld.reserve(reserve);
	volumeNew.reserve(reserve);
	volume_extraOld.reserve(reserve);
	volume_extraNew.reserve(reserve);
	
	MitOld.reserve(reserve);
	MitNew.reserve(reserve);
	
	pHiOld.reserve(reserve);
	pHiNew.reserve(reserve);
	pHOld.reserve(reserve);
	pHNew.reserve(reserve);
	
	mGinOld.reserve(reserve);
	mGinNew.reserve(reserve);
	mGextOld.reserve(reserve);
	mGextNew.reserve(reserve);
	
	mG6POld.reserve(reserve);
	mG6PNew.reserve(reserve);
	
	mO2Old.reserve(reserve);
	mO2New.reserve(reserve);
	
	StoreOld.reserve(reserve);
	StoreNew.reserve(reserve);
	
	mAinOld.reserve(reserve);
	mAinNew.reserve(reserve);
	mAextOld.reserve(reserve);
	mAextNew.reserve(reserve);
	
	mAcLinOld.reserve(reserve);
	mAcLinNew.reserve(reserve);
	mAcLextOld.reserve(reserve);
	mAcLextNew.reserve(reserve);
	
	ATPpOld.reserve(reserve);
	ATPpNew.reserve(reserve);
	
	proteinNew.reserve(reserve);
	pRbNew.reserve(reserve);
	delta_protein.reserve(reserve);
	ConcSNew.reserve(reserve);
	DNANew.reserve(reserve);
	
}


// assignment/reassignment of dynamic reserve to blood vessel vector
void CellsSystem::Set_BV_reserve(const int reserve_bv)
{
	BloodVesselVector.reserve(reserve_bv);
//   BloodVesselVector = new BloodVessel[reserve_bv];
}


// aggiunta di nuove cellule non inizializzate
//  1. incrementa ncells 
//  2. ridimensiona i vettori
void CellsSystem::AddCells( const int newcells )
{

	params.ncells += newcells;				// aggiornamento del numero di cellule
	
	name.resize(params.ncells);
	mark.resize(params.ncells);
	type.resize(params.ncells);
	
	Temperature.resize(params.ncells);
	
	phase.resize(params.ncells);
	
	death_condition.resize(params.ncells);
	age.resize(params.ncells);
	phase_age.resize(params.ncells);
	age_mother.resize(params.ncells);
	n_mitosis.resize(params.ncells);
	
	x.resize(params.ncells);
	y.resize(params.ncells);
	z.resize(params.ncells);
	
	vx.resize(params.ncells);
	vy.resize(params.ncells);
	vz.resize(params.ncells);
	
	vxnew.resize(params.ncells);
	vynew.resize(params.ncells);
	vznew.resize(params.ncells);
	
	fx.resize(params.ncells);
	fy.resize(params.ncells);
	fz.resize(params.ncells);
	
	v.resize(params.ncells);
	
	r.resize(params.ncells);
	surface.resize(params.ncells);
	volume.resize(params.ncells);
	mass.resize(params.ncells);
	
	volume_extra.resize(params.ncells);
	
	neigh.resize(params.ncells);
	vneigh.resize(params.ncells);
	vdist.resize(params.ncells);
	vcsurf.resize(params.ncells);
	gnk.resize(params.ncells);
	contact_surf.resize(params.ncells);
	
	isonCH.resize(params.ncells);
	isonAS.resize(params.ncells);
	isonBV.resize(params.ncells);
	env_surf.resize(params.ncells);
	bv_surf.resize(params.ncells);
	g_env.resize(params.ncells);
	g_bv.resize(params.ncells);
	
	M.resize(params.ncells);
	
	G.resize(params.ncells);
	G6P.resize(params.ncells);
	O2.resize(params.ncells);
	store.resize(params.ncells);
	A.resize(params.ncells);
	AcL.resize(params.ncells);
	
	h.resize(params.ncells);
	pHi.resize(params.ncells);
	
	protein.resize(params.ncells);
	prot_rate.resize(params.ncells);
	DNA.resize(params.ncells);
	DNA_rate.resize(params.ncells);
	
	GAbsRate.resize(params.ncells);
	GConsRate.resize(params.ncells);
	AAbsRate.resize(params.ncells);
	AConsRate.resize(params.ncells);
	StoreFillRate.resize(params.ncells);
	StoreConsRate.resize(params.ncells);
	AcLRate.resize(params.ncells);
	AcLOutRate.resize(params.ncells);
	O2Rate.resize(params.ncells); //*************************** new for O2 rate ******* april 2018
	
	ATP_St.resize(params.ncells);
	ATP_Ox.resize(params.ncells);
	ATP_NOx.resize(params.ncells);
	ATP2.resize(params.ncells);
	ATP3.resize(params.ncells);
	ConsATP.resize(params.ncells);
	ConsATP_1.resize(params.ncells);
	ConsATP_2.resize(params.ncells);
	ConsATP_3.resize(params.ncells);
	ConsATP_4.resize(params.ncells);
	ConsATP_5.resize(params.ncells);
	ATPtot.resize(params.ncells);
	ATPp.resize(params.ncells);
	ATPmin.resize(params.ncells);
	
	ATPstart.resize(params.ncells);
	ATPprod.resize(params.ncells);
	ATPcons.resize(params.ncells);
	
	G_extra.resize(params.ncells);
	A_extra.resize(params.ncells);
	AcL_extra.resize(params.ncells);
	
	pH.resize(params.ncells);
	SensO2.resize(params.ncells);
	ConsO.resize(params.ncells);
	
	DNA_spread.resize(params.ncells);
	
	M_T.resize(params.ncells);
	pRb.resize(params.ncells);
	
	ConcS.resize(params.ncells);
	
	cyclinD.resize(params.ncells);
	cyclinE.resize(params.ncells);
	cyclinX.resize(params.ncells);
	
	NpRbk.resize(params.ncells);


	volumeOld.resize(params.ncells);
	volumeNew.resize(params.ncells);
	volume_extraOld.resize(params.ncells);
	volume_extraNew.resize(params.ncells);
	
	MitOld.resize(params.ncells);
	MitNew.resize(params.ncells);
	
	pHiOld.resize(params.ncells);
	pHiNew.resize(params.ncells);
	pHOld.resize(params.ncells);
	pHNew.resize(params.ncells);
	
	mGinOld.resize(params.ncells);
	mGinNew.resize(params.ncells);
	mGextOld.resize(params.ncells);
	mGextNew.resize(params.ncells);
	
	mG6POld.resize(params.ncells);
	mG6PNew.resize(params.ncells);
	
	mO2Old.resize(params.ncells);
	mO2New.resize(params.ncells);
	
	StoreOld.resize(params.ncells);
	StoreNew.resize(params.ncells);
	
	mAinOld.resize(params.ncells);
	mAinNew.resize(params.ncells);
	mAextOld.resize(params.ncells);
	mAextNew.resize(params.ncells);
	
	mAcLinOld.resize(params.ncells);
	mAcLinNew.resize(params.ncells);
	mAcLextOld.resize(params.ncells);
	mAcLextNew.resize(params.ncells);
	
	ATPpOld.resize(params.ncells);
	ATPpNew.resize(params.ncells);
	
	proteinNew.resize(params.ncells);
	pRbNew.resize(params.ncells);
	delta_protein.resize(params.ncells);
	ConcSNew.resize(params.ncells);
	DNANew.resize(params.ncells);
	
}



// aggiunta di una nuova cellula inizializzata
void CellsSystem::AddInitializedCell(int& idum, CellType* cType, Environment* cEnv)
{

	// qui si incrementa il numero di cellule nel sistema (incrementa anche il contatore ncells)
	AddCells( 1 );
	
	// posizione della nuova cellula 
	unsigned long k = params.ncells - 1;
	
	// nome della nuova cellula
	name[k] = lastname;
	lastname++;
	
	// assegnazioni definite dall'utente
	type[k] = cType;
	type[k]->New_instance();										// incrementa il numero di instances del tipo cellulare dato
	
	
	// inizializzazione standard
	
	mark[k] = 0;													// cellule non marcate
	Temperature[k] = 37.;											// temperatura della cellula in gradi C
	
	x[k] = y[k] = z[k] = 0.;										// la cellula e' posizionata per default nell'origine 
	vx[k] = vy[k] = vz[k] = 0.;										// ... e ha velocita' nulla ... 

	M[k] = 100;														// numero iniziale di mitocondri
	DNA[k] = 0.;													// quantita' iniziale di DNA prodotto per la replicazione
	ATPmin[k] = (type[k]->fATPmin)*(type[k]->C2 * M[k])/type[k]->C1;	// calcolo di ATPmin
	ATPp[k] = 5; // inizializzazione dell'ATPp ad un valore sufficientemente alto (pg)

	volume[k] = type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1;	// volume iniziale
	mass[k] = type[k]->density * volume[k];							// massa iniziale
	r[k] = pow(3.*volume[k]/(4.*PI),(double)1./3.);			// calcolo del raggio corrispondente
	surface[k] = 4.*PI*r[k]*r[k];									// calcolo della superficie

	neigh[k] = 0;													// per default si crea una cellula isolata
	contact_surf[k] = 0.;
	isonCH[k] = true;
	isonAS[k] = true;
	isonBV[k] = false;
	env_surf[k] = surface[k];										// nel caso di una cellula isolata tutta la superficie e' a contatto con l'ambiente
	g_env[k] = env_surf[k]/r[k];									// fattore geometrico nel caso della cellula isolata
	g_bv[k] = 0.;
	
	
	G[k] = 0.;														// glucosio
	G6P[k] = 0.4 * G[k] * volume[k];								// G6P
	O2[k] = (cEnv->GetEnvironmentO2()/cEnv->GetEnvironmentvolume0()) * volume[k];					// ossigeno
	store[k] = STOCK_MAX * (M[k]*0.01);								// contenuto della riserva proporzionale al numero di mitocondri
	A[k] = A_CELL;													// quantità iniziale degli altri nutrienti
	AcL[k] = 0.;													// acido lattico

    h[k] = 0;                                                       // 02-dependent glucose-transport efficiency
	pHi[k] = pHi_STANDARD;
	// H = pow((double)10.,-pHi)*(volume);
	// CO2 = (cEnv->CO2/cEnv->volume0) * volume;					// anidride carbonica

	GAbsRate[k] = 0.;
	GConsRate[k] = 0.;
	AAbsRate[k] = 0.;
	AConsRate[k] = 0.;
	StoreFillRate[k] = 0.;
	StoreConsRate[k] = 0.;
	AcLRate[k] = 0.;
	AcLOutRate[k] = 0.;
	O2Rate[k] = 0.;
	

	ATP_St[k] = 0.;
	ATP_Ox[k] = 0.;
	ATP_NOx[k] = 0.;
	ATP2[k] = 0.;
	ATP3[k] = 0.;
	ConsATP[k] = 0.;
	ConsATP_1[k] = 0.;
	ConsATP_2[k] = 0.;
	ConsATP_3[k] = 0.;
	ConsATP_4[k] = 0.;
	ConsATP_5[k] = 0.;
	ATPtot[k] = 0.;

	ATPstart[k] = 0.;
	ATPprod[k] = 0.;
	ATPcons[k] = 0.;
		
	// spazio extracellulare
	volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);		// si inizializza il vol. extracell.
	
	G_extra[k] = (cEnv->GetEnvironmentG()/cEnv->GetEnvironmentvolume0()) * volume_extra[k];			// si inizializzano le quantita' a un valore corrispondente a quello ambientale
	A_extra[k] = (cEnv->GetEnvironmentA()/cEnv->GetEnvironmentvolume0()) * volume_extra[k];
	AcL_extra[k] = (cEnv->GetEnvironmentAcL()/cEnv->GetEnvironmentvolume0()) * volume_extra[k];
	// CO2_extra = (cEnv->CO2/cEnv->volume0) * volume_extra;
	
	pH[k] =  7.5443-(AcL_extra[k]/volume_extra[k])/BufCapEnv;
	// H_extra = pow((double)10.,-pH)*(volume_extra);

	SensO2[k] = 1.;     // efficienza di consumo dell'ossigeno
	ConsO[k] = 0.;      // consumo ossigeno inizializzato a zero
	// ProdCO2 = 0.;	// produzione anidride carbonica inizializzata a zero
	

	phase[k] = G1m_phase;											// fase cellulare iniziale
	//phase[k] = G0_phase;											// fase cellulare iniziale
	
	death_condition[k] = 0.;										// inizialmente la cellula e' viva e quindi non si registra nessuna condizione di morte
	age[k] = 0;														// eta' della cellula (appena nata, ueeeeeee ....)
	phase_age[k] = 0;												// eta' dello stato cellulare attuale
	age_mother[k] = 0.;												// eta' della madre al momento della mitosi (inizializzata a 0, che vuol dire che non c'e' madre per questa cellula)
	n_mitosis[k] = 0;												// numero di mitosi dall'inizio della simulazione

	DNA_spread[k] = type[k]->DNA_MAX_SPREAD*(2*ran2()-1.);		// fluttuazione della velocita' di sintesi del DNA
	M_T[k] = type[k]->M_T_MEAN;										// durata fase M
		
	ConcS[k] = type[k]->ConcS_0;									// inizializzazione della concentrazione di S
	
	pRb[k] = pRb_STANDARD;											// massa iniziale della pRb in kg
	
	protein[k] = 0.;												// inizialmente non ci sono le altre proteine
	prot_rate[k] = 0.;
	DNA_rate[k] = 0.;
	cyclinD[k] = 0.;
	cyclinE[k] = 0.;
	cyclinX[k] = 0.;
	

	NpRbk[k] = 0.;

}

// metodo di copia: copia tutto meno le caratteristiche geometriche-topologiche e le variabili temporanee (che restano non-inizializzate)
int CellsSystem::CopyCell( const unsigned long int k, const unsigned long int kstart, const unsigned long int kstop)
{

	// controllo di consistenza dell'indice k
	if( kstart > kstop || kstop > params.ncells-1)
	 return -1;
	
	// copia del contenuto della cellula k-esima in tutte le celle da kstart a kstop (estremi inclusi)
	for(unsigned long int kk=kstart; kk <= kstop; kk++)
		{
		// il nome della cellula viene definito automaticamente
		name[kk] = lastname;
		lastname++;
		
		mark[kk] = mark[k];
		type[kk] = type[k];
		type[kk]->New_instance();		// incrementa il numero di instances del tipo cellulare dato
		
		Temperature[kk] = Temperature[k];
		
		phase[kk] = phase[k];
		
		death_condition[kk] = death_condition[k];
		age[kk] = age[k];
		phase_age[kk] = phase_age[k];
		age_mother[kk] = age_mother[k];
		n_mitosis[kk] = n_mitosis[k];
		
		x[kk] = x[k];
		y[kk] = y[k];
		z[kk] = z[k];
		
		vx[kk] = vx[k];
		vy[kk] = vy[k];
		vz[kk] = vz[k];
		
		r[kk] = r[k];
		surface[kk] = surface[k];
		volume[kk] = volume[k];
		mass[kk] = mass[k];
		
		volume_extra[kk] = volume_extra[k];
		
		neigh[kk] = neigh[k];
		vneigh[kk] = vneigh[k];
		vdist[kk] = vdist[k];
		vcsurf[kk] = vcsurf[k];
		gnk[kk] = gnk[k];
		contact_surf[kk] = contact_surf[k];
		
		isonCH[kk] = isonCH[k];
		isonAS[kk] = isonAS[k];
		isonBV[kk] = isonBV[k];
		env_surf[kk] = env_surf[k];
		g_env[kk] = g_env[k];
		bv_surf[kk] = bv_surf[k];
		g_bv[kk] = g_bv[k];
		
		M[kk] = M[k];
		
		G[kk] = G[k];
		G6P[kk] = G6P[k];
		O2[kk] = O2[k];
		store[kk] = store[k];
		A[kk] = A[k];
		AcL[kk] = AcL[k];
		
        h[kk] = h[k];
		pHi[kk] = pHi[k];
		
		protein[kk] = protein[k];
		prot_rate[kk] = prot_rate[k];
		DNA[kk] = DNA[k];
		DNA_rate[kk] = DNA_rate[k];
		
		GAbsRate[kk] = GAbsRate[k];
		GConsRate[kk] = GConsRate[k];
		AAbsRate[kk] = AAbsRate[k];
		AConsRate[kk] = AConsRate[k];
		StoreFillRate[kk] = StoreFillRate[k];
		StoreConsRate[kk] = StoreConsRate[k];
		AcLRate[kk] = AcLRate[k];
		AcLOutRate[kk] = AcLOutRate[k];
		O2Rate[kk] = O2Rate[k];
		
		ATP_St[kk] = ATP_St[k];
		ATP_Ox[kk] = ATP_Ox[k];
		ATP_NOx[kk] = ATP_NOx[k];
		ATP2[kk] = ATP2[k];
		ATP3[kk] = ATP3[k];
		ConsATP[kk] = ConsATP[k];
		ConsATP_1[kk] = ConsATP_1[k];
		ConsATP_2[kk] = ConsATP_2[k];
		ConsATP_3[kk] = ConsATP_3[k];
		ConsATP_4[kk] = ConsATP_4[k];
		ConsATP_5[kk] = ConsATP_5[k];
		ATPtot[kk] = ATPtot[k];
		ATPp[kk] = ATPp[k];
		ATPmin[kk] = ATPmin[k];
		
		ATPstart[kk] = ATPstart[k];
		ATPprod[kk] = ATPprod[k];
		ATPcons[kk] = ATPcons[k];
		
		G_extra[kk] = G_extra[k];
		A_extra[kk] = A_extra[k];
		AcL_extra[kk] = AcL_extra[k];
		
		pH[kk] = pH[k];
		SensO2[kk] = SensO2[k];
		ConsO[kk] = ConsO[k];
		
		DNA_spread[kk] = DNA_spread[k];
		
		M_T[kk] = M_T[k];
		pRb[kk] = pRb[k];
		
		ConcS[kk] = ConcS[k];
		
		cyclinD[kk] = cyclinD[k];
		cyclinE[kk] = cyclinE[k];
		cyclinX[kk] = cyclinX[k];
		
		NpRbk[kk] = NpRbk[k];
		}
	
	return kstop-kstart+1;

}

// Copy method: copy all the geometric-topological features and temporary variables (which remain non-initialized), also modifies the cell type
int CellsSystem::CopyCell( const unsigned long int k, const unsigned long int kstart, const unsigned long int kstop, CellType* newtype)
{
    
    // controllo di consistenza dell'indice k
    if( kstart > kstop || kstop > params.ncells-1)
        return -1;
    
    // copia del contenuto della cellula k-esima in tutte le celle da kstart a kstop (estremi inclusi)
    for(unsigned long int kk=kstart; kk <= kstop; kk++)
    {
        // il nome della cellula viene definito automaticamente
        name[kk] = lastname;
        lastname++;
        
        mark[kk] = mark[k];
        type[kk] = newtype;
        type[kk]->New_instance();		// incrementa il numero di instances del tipo cellulare dato
        
        Temperature[kk] = Temperature[k];
        
        phase[kk] = phase[k];
        
        death_condition[kk] = death_condition[k];
        age[kk] = age[k];
        phase_age[kk] = phase_age[k];
        age_mother[kk] = age_mother[k];
        n_mitosis[kk] = n_mitosis[k];
        
        x[kk] = x[k];
        y[kk] = y[k];
        z[kk] = z[k];
        
        vx[kk] = vx[k];
        vy[kk] = vy[k];
        vz[kk] = vz[k];
        
        r[kk] = r[k];
        surface[kk] = surface[k];
        volume[kk] = volume[k];
        mass[kk] = mass[k];
        
        volume_extra[kk] = volume_extra[k];
        
        neigh[kk] = neigh[k];
        vneigh[kk] = vneigh[k];
        vdist[kk] = vdist[k];
        vcsurf[kk] = vcsurf[k];
        gnk[kk] = gnk[k];
        contact_surf[kk] = contact_surf[k];
        
        isonCH[kk] = isonCH[k];
        isonAS[kk] = isonAS[k];
        isonBV[kk] = isonBV[k];
        env_surf[kk] = env_surf[k];
        g_env[kk] = g_env[k];
        bv_surf[kk] = bv_surf[k];
        g_bv[kk] = g_bv[k];
        
        M[kk] = M[k];
        
        G[kk] = G[k];
        G6P[kk] = G6P[k];
        O2[kk] = O2[k];
        store[kk] = store[k];
        A[kk] = A[k];
        AcL[kk] = AcL[k];
        
        h[kk] = h[k];
        pHi[kk] = pHi[k];
        
        protein[kk] = protein[k];
        prot_rate[kk] = prot_rate[k];
        DNA[kk] = DNA[k];
        DNA_rate[kk] = DNA_rate[k];
        
        GAbsRate[kk] = GAbsRate[k];
        GConsRate[kk] = GConsRate[k];
        AAbsRate[kk] = AAbsRate[k];
        AConsRate[kk] = AConsRate[k];
        StoreFillRate[kk] = StoreFillRate[k];
        StoreConsRate[kk] = StoreConsRate[k];
        AcLRate[kk] = AcLRate[k];
        AcLOutRate[kk] = AcLOutRate[k];
        O2Rate[kk] = O2Rate[k];
        
        ATP_St[kk] = ATP_St[k];
        ATP_Ox[kk] = ATP_Ox[k];
        ATP_NOx[kk] = ATP_NOx[k];
        ATP2[kk] = ATP2[k];
        ATP3[kk] = ATP3[k];
        ConsATP[kk] = ConsATP[k];
        ConsATP_1[kk] = ConsATP_1[k];
        ConsATP_2[kk] = ConsATP_2[k];
        ConsATP_3[kk] = ConsATP_3[k];
        ConsATP_4[kk] = ConsATP_4[k];
        ConsATP_5[kk] = ConsATP_5[k];
        ATPtot[kk] = ATPtot[k];
        ATPp[kk] = ATPp[k];
        ATPmin[kk] = ATPmin[k];
        
        ATPstart[kk] = ATPstart[k];
        ATPprod[kk] = ATPprod[k];
        ATPcons[kk] = ATPcons[k];
        
        G_extra[kk] = G_extra[k];
        A_extra[kk] = A_extra[k];
        AcL_extra[kk] = AcL_extra[k];
        
        pH[kk] = pH[k];
        SensO2[kk] = SensO2[k];
        ConsO[kk] = ConsO[k];
        
        DNA_spread[kk] = DNA_spread[k];
        
        M_T[kk] = M_T[k];
        pRb[kk] = pRb[k];
        
        ConcS[kk] = ConcS[k];
        
        cyclinD[kk] = cyclinD[k];
        cyclinE[kk] = cyclinE[k];
        cyclinX[kk] = cyclinX[k];
        
        NpRbk[kk] = NpRbk[k];
    }
    
    return kstop-kstart+1;
    
}


// metodo di copia: replica la cellula k-esima inserendone una copia alla fine (tutto tranne le caratteristiche geometriche-topologiche)
int CellsSystem::ReplicateCell( const unsigned long int k )
{
	
	// controllo di consistenza dell'indice k
	//if( k < 0 || k > ncells-1 )
	if( k > params.ncells-1 )
		return -1;
	
	// AddCells inserisce una cellula non inizializzata ed incrementa di 1 il contatore del numero di cellule ncells
	AddCells( 1 );								// in questa istruzione ncells -> ncells+1
	int ic = CopyCell( k, params.ncells-1, params.ncells-1);	// si copia la cellula k-esima in posizione ncells-1
	
	if(ic < 1) 
		return ic;
	else
		return 1;

}

// metodo di copia: replica la cellula k-esima inserendone una copia alla fine (tutto tranne le caratteristiche geometriche-topologiche), inoltre questa versione cambia il tipo cellulare
int CellsSystem::ReplicateCell( const unsigned long int k, CellType* newtype )
{
    
    // controllo di consistenza dell'indice k
    //if( k < 0 || k > ncells-1 )
    if( k > params.ncells-1 )
        return -1;
    
    // AddCells inserisce una cellula non inizializzata ed incrementa di 1 il contatore del numero di cellule ncells
    AddCells( 1 );								// in questa istruzione ncells -> ncells+1
    int ic = CopyCell( k, params.ncells-1, params.ncells-1, newtype);	// si copia la cellula k-esima in posizione ncells-1, con un diverso tipo cellulare
    
    if(ic < 1)
        return ic;
    else
        return 1;
    
}


// metodo di copia: replica la cellula k-esima inserendo n copie alla fine (tutto tranne le caratteristiche geometriche-topologiche)
int CellsSystem::ReplicateCell( const unsigned long int k, const unsigned long int n )
{

	//if( k < 0 || k > ncells-1 )
	if( k > params.ncells-1 )
		return -1;
		
	if( n < 1 )
		return 0;
	 
	int kstart = params.ncells;
	int kstop = params.ncells + n - 1;
	
	AddCells( n );
	CopyCell( k, kstart, kstop);
	
	return n;

}



// stampa su file di una singola cellula

void CellsSystem::PrintCell( std::ostream& stream, const unsigned long int k )
{
	stream << std::setprecision(10);
	stream << "name: " << name[k] << std::endl;
	stream << "mark: " << mark[k] << std::endl;
	stream << "type: " << type[k] << std::endl; // stampa del tipo cellulare
	stream << std::endl;
	
	stream << "temperature: " << Temperature[k] << std::endl;
	stream << std::endl;
	
	stream << "phase: " <<  phase[k] << std::endl;
	stream << "death_condition: " <<  death_condition[k] << std::endl;
	stream << "age: " <<  age[k] << std::endl;
	stream << "phase_age: " <<  phase_age[k] << std::endl;
	stream << "age_mother: " <<  age_mother[k] << std::endl;
	stream << "n_mitosis: " <<  n_mitosis[k] << std::endl;
	stream << std::endl;
	
	stream << "x: " <<  x[k] << std::endl;
	stream << "y: " <<  y[k] << std::endl;
	stream << "z: " <<  z[k] << std::endl;
	stream << "vx: " <<  vx[k] << std::endl;
	stream << "vy: " <<  vy[k] << std::endl;
	stream << "vz: " <<  vz[k] << std::endl;
	stream << std::endl;
	stream << "r: " <<  r[k] << std::endl;	
	stream << "surface: " <<  surface[k] << std::endl;
	stream << "volume: " <<  volume[k] << std::endl;
	stream << "mass: " <<  mass[k] << std::endl;
	stream << std::endl;
	stream << "volume_extra: " << volume_extra[k] << std::endl;
	stream << std::endl;
	stream << "neigh: " <<  neigh[k] << std::endl;
	stream << "- lista di vicini e relative superfici di contatto, distanze, fattori geometrici:" << std::endl;
	for(int n=0; n<neigh[k]; n++)
		stream << "- \t"<< vneigh[k][n] << "\t" << vcsurf[k][n] << "\t" << vdist[k][n] << "\t" << gnk[k][n] << std::endl;
	stream << "contact_surf: " <<  contact_surf[k] << std::endl;	
	stream << std::endl;
	stream << "isonCH: " <<  isonCH[k] << std::endl;
	stream << "isonAS: " <<  isonAS[k] << std::endl;
	stream << "isonBV: " <<  isonBV[k] << std::endl;
	stream << "env_surf: " <<  env_surf[k] << std::endl;
	stream << "g_env: " <<  g_env[k] << std::endl;
	stream << "bv_surf: " <<  bv_surf[k] << std::endl;
	stream << "g_bv: " <<  g_bv[k] << std::endl;
	stream << std::endl;
	stream << "M: " <<  M[k] << std::endl;
	stream << std::endl;
	stream << "G: " <<  G[k] << std::endl;
	stream << "G6P: " <<  G6P[k] << std::endl;
	stream << "O2: " <<  O2[k] << std::endl;
	stream << "store: " <<  store[k] << std::endl;
	stream << "A: " <<  A[k] << std::endl;
	stream << "AcL: " <<  AcL[k] << std::endl;
	stream << std::endl;
	stream << "h: " <<  h[k] << std::endl;
	stream << "pHi: " <<  pHi[k] << std::endl;
	// stream << "H: " <<  H[k] << endl;
	// stream << "CO2: " <<  CO2[k] << endl;
	stream << std::endl;
	stream << "protein: " <<  protein[k] << std::endl;
	stream << "prot_rate: " << prot_rate[k] << std::endl;
	stream << "DNA: " <<  DNA[k] << std::endl;
	stream << "DNA_rate: " <<  DNA_rate[k] << std::endl;
	stream << std::endl;
	stream << "GAbsRate: " <<  GAbsRate[k] << std::endl;
	stream << "GConsRate: " <<  GConsRate[k] << std::endl;
	stream << "AAbsRate: " <<  AAbsRate[k] << std::endl;
	stream << "AConsRate: " <<  AConsRate[k] << std::endl;
	stream << "StoreFillRate: " <<  StoreFillRate[k] << std::endl;
	stream << "StoreConsRate: " <<  StoreConsRate[k] << std::endl;
	stream << "AcLRate: " << AcLRate[k] << std::endl;
	stream << "AcLOutRate: " <<  AcLOutRate[k] << std::endl;
	stream << std::endl;
	stream << "ATP_St: " <<  ATP_St[k] << std::endl;
	stream << "ATP_Ox: " <<  ATP_Ox[k] << std::endl;
	stream << "ATP_NOx: " <<  ATP_NOx[k] << std::endl;
	stream << "ATP2: " <<  ATP2[k] << std::endl;
	stream << "ATP3: " <<  ATP3[k] << std::endl;
	stream << "ConsATP: " <<  ConsATP[k] << std::endl;
	stream << "ConsATP_1: " <<  ConsATP_1[k] << std::endl;
	stream << "ConsATP_2: " <<  ConsATP_2[k] << std::endl;
	stream << "ConsATP_3: " <<  ConsATP_3[k] << std::endl;
	stream << "ConsATP_4: " <<  ConsATP_4[k] << std::endl;
	stream << "ConsATP_5: " <<  ConsATP_5[k] << std::endl;
	stream << "ATPtot: " <<  ATPtot[k] << std::endl;
	stream << "ATPp: " <<  ATPp[k] << std::endl;
	stream << "ATPmin: " <<  ATPmin[k] << std::endl;
	stream << std::endl;
	stream << "ATPstart: " <<  ATPstart[k] << std::endl;
	stream << "ATPprod: " <<  ATPprod[k] << std::endl;
	stream << "ATPcons: " <<  ATPcons[k] << std::endl;
	stream << std::endl;
	stream << "G_extra: " <<  G_extra[k] << std::endl;
	stream << "A_extra: " <<  A_extra[k] << std::endl;
	stream << "AcL_extra: " <<  AcL_extra[k] << std::endl;
	stream << std::endl;
	stream << "pH: " <<  pH[k] << std::endl;
	// stream << "H_extra: " <<  H_extra[k] << endl;
	// stream << "CO2_extra: " <<  CO2_extra[k] << endl;
	stream << std::endl;
	stream << "SensO2: " <<  SensO2[k] << std::endl;
	stream << "ConsO: " <<  ConsO[k] << std::endl;
	// stream << "ProdCO2: " <<  ProdCO2[k] << endl;
	stream << std::endl;
	stream << "DNA_spread: " <<  DNA_spread[k] << std::endl;
	stream << std::endl;
	stream << "M_T: " <<  M_T[k] << std::endl;
	stream << "pRb: " <<  pRb[k] << std::endl;
	stream << std::endl;
	stream << "ConcS: " <<  ConcS[k] << std::endl;
	stream << std::endl;
	stream << "cyclinD: " <<  cyclinD[k] << std::endl;
	stream << "cyclinE: " <<  cyclinE[k] << std::endl;
	stream << "cyclinX: " <<  cyclinX[k] << std::endl;
	stream << std::endl;
	stream << "NpRbk: " <<  NpRbk[k] << std::endl;
	stream << std::endl;
	
}


// stampa dei dati sotto forma di un'unica stringa leggibile da un programma spreadsheet
//
void CellsSystem::PrintCellData( const unsigned long int k, std::ofstream& stream, long int nrec )
{

	static bool first_print=true;

// stampa degli header nel caso che questa sia la prima volta che si stampa
	if(first_print)
		{
		first_print = false;

	stream << "n \t" \
			<< "name \t" \
			<< "mark \t" \
			<< "T \t" \
			<< "phase \t" \
			<< "death_condition \t" \
			<< "age \t" \
			<< "phase_age \t" \
			<< "age_mother \t" \
			<< "n_mitosis \t" \
			<< "x \t" \
			<< "y \t" \
			<< "z \t" \
			<< "vx \t" \
			<< "vy \t" \
			<< "vz \t" \
			<< "r \t" \
			<< "surface \t" \
			<< "volume \t" \
			<< "mass \t" \
			<< "volume_extra \t" \
			<< "neigh \t" \
			<< "contact_surf \t" \
			<< "isonCH \t" \
			<< "isonAS \t" \
			<< "isonBV \t" \
			<< "env_surf \t" \
			<< "g_env \t" \
			<< "bv_surf \t" \
			<< "g_bv \t" \
			<< "M \t" \
			<< "G \t" \
			<< "G6P \t" \
			<< "O2 \t" \
			<< "store \t" \
			<< "A \t" \
			<< "AcL \t" \
			<< "h \t" \
			<< "pHi \t" \
			<< "protein \t" \
			<< "prot_rate \t" \
			<< "DNA \t" \
			<< "DNA_rate \t" \
			<< "GAbsRate \t" \
			<< "GConsRate \t" \
			<< "AAbsRate \t" \
			<< "AConsRate \t" \
			<< "StoreFillRate \t" \
			<< "StoreConsRate \t" \
			<< "AcLRate \t" \
			<< "AcLOutRate \t" \
			<< "O2Rate \t" \
			<< "ATP_St \t" \
			<< "ATP_Ox \t" \
			<< "ATP_NOx \t" \
			<< "ATP2 \t" \
			<< "ATP3 \t" \
			<< "ConsATP \t" \
			<< "ConsATP_1 \t" \
			<< "ConsATP_2 \t" \
			<< "ConsATP_3 \t" \
			<< "ConsATP_4 \t" \
			<< "ConsATP_5 \t" \
			<< "ATPtot \t" \
			<< "ATPp \t" \
			<< "ATPmin \t" \
			<< "ATPstart \t" \
			<< "ATPprod \t" \
			<< "ATPcons \t" \
			<< "G_extra \t" \
			<< "A_extra \t" \
			<< "AcL_extra \t" \
			<< "pH \t" \
			<< "SensO2 \t" \
			<< "ConsO \t" \
			<< "DNA_spread \t" \
			<< "M_T \t" \
			<< "pRb \t" \
			<< "ConcS \t" \
			<< "cyclinD \t" \
			<< "cyclinE \t" \
			<< "cyclinX \t" \
			<< "NpRbk" << std::endl;		

/*
righe temporanemente eliminate
			<< "H \t" \
			<< "CO2 \t" \
			<< "ProdCO2 \t" \
			<< "H_extra \t" \
			<< "CO2_extra \t" \


*/
		
		}


	stream << std::setprecision(8) << std::scientific \
			<< nrec << "\t" \
			<< name[k] << "\t" \
			<< mark[k] << "\t" \
			<< Temperature[k] << "\t" \
			<< phase[k] << "\t" \
			<< death_condition[k] << "\t" \
			<< age[k] << "\t" \
			<< phase_age[k] << "\t" \
			<< age_mother[k] << "\t" \
			<< n_mitosis[k] << "\t" \
			<< x[k] << "\t" \
			<< y[k] << "\t" \
			<< z[k] << "\t" \
			<< vx[k] << "\t" \
			<< vy[k] << "\t" \
			<< vz[k] << "\t" \
			<< r[k] << "\t" \
			<< surface[k] << "\t" \
			<< volume[k] << "\t" \
			<< mass[k] << "\t" \
			<< volume_extra[k] << "\t" \
			<< neigh[k] << "\t" \
			<< contact_surf[k] << "\t" \
			<< isonCH[k] << "\t" \
			<< isonAS[k] << "\t" \
			<< isonBV[k] << "\t" \
			<< env_surf[k] << "\t" \
			<< g_env[k] << "\t" \
			<< bv_surf[k] << "\t" \
			<< g_bv[k] << "\t" \
			<< M[k] << "\t" \
			<< G[k] << "\t" \
			<< G6P[k] << "\t" \
			<< O2[k] << "\t" \
			<< store[k] << "\t" \
			<< A[k] << "\t" \
			<< AcL[k] << "\t" \
			<< h[k] << "\t" \
			<< pHi[k] << "\t" \
			<< protein[k] << "\t" \
			<< prot_rate[k] << "\t" \
			<< DNA[k] << "\t" \
			<< DNA_rate[k] << "\t" \
			<< GAbsRate[k] << "\t" \
			<< GConsRate[k] << "\t" \
			<< AAbsRate[k] << "\t" \
			<< AConsRate[k] << "\t" \
			<< StoreFillRate[k] << "\t" \
			<< StoreConsRate[k] << "\t" \
			<< AcLRate[k] << "\t" \
			<< AcLOutRate[k] << "\t" \
			<< O2Rate[k] << "\t" \
			<< ATP_St[k] << "\t" \
			<< ATP_Ox[k] << "\t" \
			<< ATP_NOx[k] << "\t" \
			<< ATP2[k] << "\t" \
			<< ATP3[k] << "\t" \
			<< ConsATP[k] << "\t" \
			<< ConsATP_1[k] << "\t" \
			<< ConsATP_2[k] << "\t" \
			<< ConsATP_3[k] << "\t" \
			<< ConsATP_4[k] << "\t" \
			<< ConsATP_5[k] << "\t" \
			<< ATPtot[k] << "\t" \
			<< ATPp[k] << "\t" \
			<< ATPmin[k] << "\t" \
			<< ATPstart[k] << "\t" \
			<< ATPprod[k] << "\t" \
			<< ATPcons[k] << "\t" \
			<< G_extra[k] << "\t" \
			<< A_extra[k] << "\t" \
			<< AcL_extra[k] << "\t" \
			<< pH[k] << "\t" \
			<< SensO2[k] << "\t" \
			<< ConsO[k] << "\t" \
			<< DNA_spread[k] << "\t" \
			<< M_T[k] << "\t" \
			<< pRb[k] << "\t" \
			<< ConcS[k] << "\t" \
			<< cyclinD[k] << "\t" \
			<< cyclinE[k] << "\t" \
			<< cyclinX[k] << "\t" \
			<< NpRbk[k] << std::endl;			

/*
righe temporanemente eliminate
			<< H << "\t" \
			<< CO2 << "\t" \
			<< ProdCO2 << "\t" \
			<< "H_extra \t" \
			<< "CO2_extra \t" \


*/

}

//Removing a cell corresponding to the nth position
// 1. Scales the vectors using the erase method (this is not the most efficient method, but its functionality is clear,
// and for the moment I use this)
// 2. decreases ncells
void CellsSystem::RemoveCell( const unsigned long int n )
{

	name.erase(name.begin()+n);
	mark.erase(mark.begin()+n);
	type.erase(type.begin()+n);
	
	Temperature.erase(Temperature.begin()+n);
	
	phase.erase(phase.begin()+n);
	
	death_condition.erase(death_condition.begin()+n);
	age.erase(age.begin()+n);
	phase_age.erase(phase_age.begin()+n);
	age_mother.erase(age_mother.begin()+n);
	n_mitosis.erase(n_mitosis.begin()+n);
	
	x.erase(x.begin()+n);
	y.erase(y.begin()+n);
	z.erase(z.begin()+n);
	
	vx.erase(vx.begin()+n);
	vy.erase(vy.begin()+n);
	vz.erase(vz.begin()+n);
	
	vxnew.erase(vxnew.begin()+n);
	vynew.erase(vynew.begin()+n);
	vznew.erase(vznew.begin()+n);
	
	fx.erase(fx.begin()+n);
	fy.erase(fy.begin()+n);
	fz.erase(fz.begin()+n);
	
	v.erase(v.begin()+n);
	
	r.erase(r.begin()+n);
	surface.erase(surface.begin()+n);
	volume.erase(volume.begin()+n);
	mass.erase(mass.begin()+n);
	
	volume_extra.erase(volume_extra.begin()+n);
	
	neigh.erase(neigh.begin()+n);
	vneigh.erase(vneigh.begin()+n);
	vdist.erase(vdist.begin()+n);
	vcsurf.erase(vcsurf.begin()+n);
	gnk.erase(gnk.begin()+n);
	contact_surf.erase(contact_surf.begin()+n);
	
	isonCH.erase(isonCH.begin()+n);
	isonAS.erase(isonAS.begin()+n);
	isonBV.erase(isonBV.begin()+n);
	env_surf.erase(env_surf.begin()+n);
	g_env.erase(g_env.begin()+n);
	bv_surf.erase(bv_surf.begin()+n);
	g_bv.erase(g_bv.begin()+n);
	
	M.erase(M.begin()+n);
	
	G.erase(G.begin()+n);
	G6P.erase(G6P.begin()+n);
	O2.erase(O2.begin()+n);
	store.erase(store.begin()+n);
	A.erase(A.begin()+n);
	AcL.erase(AcL.begin()+n);
	
	h.erase(h.begin()+n);
	pHi.erase(pHi.begin()+n);
	
	protein.erase(protein.begin()+n);
	prot_rate.erase(prot_rate.begin()+n);
	DNA.erase(DNA.begin()+n);
	DNA_rate.erase(DNA_rate.begin()+n);
	
	GAbsRate.erase(GAbsRate.begin()+n);
	GConsRate.erase(GConsRate.begin()+n);
	AAbsRate.erase(AAbsRate.begin()+n);
	AConsRate.erase(AConsRate.begin()+n);
	StoreFillRate.erase(StoreFillRate.begin()+n);
	StoreConsRate.erase(StoreConsRate.begin()+n);
	AcLRate.erase(AcLRate.begin()+n);
	AcLOutRate.erase(AcLOutRate.begin()+n);
	O2Rate.erase(O2Rate.begin()+n); //*************************** new for O2 rate ******* april 2018
	
	ATP_St.erase(ATP_St.begin()+n);
	ATP_Ox.erase(ATP_Ox.begin()+n);
	ATP_NOx.erase(ATP_NOx.begin()+n);
	ATP2.erase(ATP2.begin()+n);
	ATP3.erase(ATP3.begin()+n);
	ConsATP.erase(ConsATP.begin()+n);
	ConsATP_1.erase(ConsATP_1.begin()+n);
	ConsATP_2.erase(ConsATP_2.begin()+n);
	ConsATP_3.erase(ConsATP_3.begin()+n);
	ConsATP_4.erase(ConsATP_4.begin()+n);
	ConsATP_5.erase(ConsATP_5.begin()+n);
	ATPtot.erase(ATPtot.begin()+n);
	ATPp.erase(ATPp.begin()+n);
	ATPmin.erase(ATPmin.begin()+n);
	
	ATPstart.erase(ATPstart.begin()+n);
	ATPprod.erase(ATPprod.begin()+n);
	ATPcons.erase(ATPcons.begin()+n);
	
	G_extra.erase(G_extra.begin()+n);
	A_extra.erase(A_extra.begin()+n);
	AcL_extra.erase(AcL_extra.begin()+n);
	
	pH.erase(pH.begin()+n);
	SensO2.erase(SensO2.begin()+n);
	ConsO.erase(ConsO.begin()+n);
	
	DNA_spread.erase(DNA_spread.begin()+n);
	
	M_T.erase(M_T.begin()+n);
	pRb.erase(pRb.begin()+n);
	
	ConcS.erase(ConcS.begin()+n);
	
	cyclinD.erase(cyclinD.begin()+n);
	cyclinE.erase(cyclinE.begin()+n);
	cyclinX.erase(cyclinX.begin()+n);
	
	NpRbk.erase(NpRbk.begin()+n);


	volumeOld.erase(volumeOld.begin()+n);
	volumeNew.erase(volumeNew.begin()+n);
	volume_extraOld.erase(volume_extraOld.begin()+n);
	volume_extraNew.erase(volume_extraNew.begin()+n);
	
	MitOld.erase(MitOld.begin()+n);
	MitNew.erase(MitNew.begin()+n);
	
	pHiOld.erase(pHiOld.begin()+n);
	pHiNew.erase(pHiNew.begin()+n);
	pHOld.erase(pHOld.begin()+n);
	pHNew.erase(pHNew.begin()+n);
	
	mGinOld.erase(mGinOld.begin()+n);
	mGinNew.erase(mGinNew.begin()+n);
	mGextOld.erase(mGextOld.begin()+n);
	mGextNew.erase(mGextNew.begin()+n);
	
	mG6POld.erase(mG6POld.begin()+n);
	mG6PNew.erase(mG6PNew.begin()+n);
	
	mO2Old.erase(mO2Old.begin()+n);
	mO2New.erase(mO2New.begin()+n);
	
	StoreOld.erase(StoreOld.begin()+n);
	StoreNew.erase(StoreNew.begin()+n);
	
	mAinOld.erase(mAinOld.begin()+n);
	mAinNew.erase(mAinNew.begin()+n);
	mAextOld.erase(mAextOld.begin()+n);
	mAextNew.erase(mAextNew.begin()+n);
	
	mAcLinOld.erase(mAcLinOld.begin()+n);
	mAcLinNew.erase(mAcLinNew.begin()+n);
	mAcLextOld.erase(mAcLextOld.begin()+n);
	mAcLextNew.erase(mAcLextNew.begin()+n);
	
	ATPpOld.erase(ATPpOld.begin()+n);
	ATPpNew.erase(ATPpNew.begin()+n);
	
	proteinNew.erase(proteinNew.begin()+n);
	pRbNew.erase(pRbNew.begin()+n);
	delta_protein.erase(delta_protein.begin()+n);
	ConcSNew.erase(ConcSNew.begin()+n);
	DNANew.erase(DNANew.begin()+n);
	
	params.ncells--;							// aggiornamento del numero di cellule
	
}

void CellsSystem::Add_BloodVesselVector(BloodVessel NewBV)
{
  
  //Vector()[index] = NewBV;
  BloodVesselVector[nbv] = NewBV;
  nbv++;
  //BloodVesselVector[index] = NewBV; /*cout << "New blood vessel in CellsSystem" << endl;*/
  //bloodVesselMap[index] = NewBV;
}
void CellsSystem::clean_BloodVesselVector()
{
  this->nbv=0;
  //std::fill(std::begin(BloodVesselVector),std::end(BloodVesselVector), BloodVessel());
  //BloodVesselVector.resize(Get_ncells());
}
//  ******************** Timing ********************
//
// funzione che restituisce il tempo trascorso in secondi
//
// 
double CellsSystem::Timing( bool reset )
{

	if(reset) 
		{
		time(&time_old);
		timing = 0.;
		}
	else
		{
		time(&time_now);
		timing = difftime(time_now, time_old);
		time_old = time_now;
		}
		
	return timing;

}
// double Get_dt() { return dt; };
// double Get_dt_sm() { return dt_sm; };
// bool Get_slow_motion() { return slow_motion; };
// double Get_t() { return t; };
// double Get_treal() { return treal; };
// double Get_tmax() { return tmax; };
bool CellsSystem::Get_slow_motion() 
{ 
  return params.slow_motion; 
};
double CellsSystem::Get_dt()
{
  return params.dt;
}
void CellsSystem::Set_dt( double newdt ) 
{ 
  params.dt = newdt;
};
double CellsSystem::Get_time_from_CGAL()
{ 
  return time_from_CGAL; 
};
unsigned long CellsSystem::Get_nstep() 
{ 
  return params.nstep;
};
unsigned long CellsSystem::Get_nstep_start() 
{ 
  return params.nstep_start; 
};
unsigned long CellsSystem::Get_nmax() 
{ 
  return params.nmax; 
};
void CellsSystem::Set_idum ( int newidum ) 
{ 
  params.idum = newidum; 
};

int CellsSystem::Get_idum() 
{ 
  return params.idum; 
};
unsigned long CellsSystem::Get_nprint() 
{ 
  return params.nprint; 
};
unsigned long CellsSystem::Get_nscreen() 
{ 
  return params.nscreen; 
};
unsigned long CellsSystem::Get_nconfiguration() 
{ 
  return params.nconfiguration; 
};
void CellsSystem::Set_nconfiguration( unsigned long newnconfiguration ) 
{ 
  params.nconfiguration = newnconfiguration; 
};
void CellsSystem::Step_nconfiguration() 
{ 
  params.nconfiguration++; 
};
double CellsSystem::Get_eps() 
{ 
  return params.eps;
};
double CellsSystem::Get_delta_vmax()
{ 
  return params.delta_vmax; 
};
void CellsSystem::Set_eps( double neweps ) 
{ 
  params.eps = neweps; 
};
void CellsSystem::Set_delta_vmax( double newdelta_vmax ) 
{ 
  params.delta_vmax = newdelta_vmax; 
};
unsigned long CellsSystem::Get_ncells() 
{ 
  return params.ncells; 
};
unsigned long CellsSystem::Get_alive() 
{ 
  return params.alive; 
};
unsigned long CellsSystem::Get_ntypes() 
{ 
  return params.ntypes; 
};
bool CellsSystem::Get_flowON() 
{ 
  return params.flowON; 
};
bool CellsSystem::Get_doseON() 
{ 
  return params.doseON;
};

boost::property_tree::ptree vbl::ReadInParameters::as_ptree() const
{
  boost::property_tree::ptree pt;
#define DOPT(name_buffer) pt.put(#name_buffer, name_buffer);
  DOPT(sim_type)
  DOPT(run)
  DOPT(dt)
  DOPT(dt_sm)
  DOPT(t)
  DOPT(t_ini)
  DOPT(treal)
  DOPT(tmax)
  DOPT(tsm_start)
  DOPT(tsm_stop)
  DOPT(slow_motion)
  DOPT(t_CPU_max)
  DOPT(nstep)
  DOPT(nstep_start)
  DOPT(nmax)
  DOPT(idum)
  DOPT(nprint)
  DOPT(nscreen)
  DOPT(nconfiguration)
  DOPT(eps)
  DOPT(delta_vmax)
  DOPT(ncells)
  DOPT(alive)
  DOPT(ntypes)
  DOPT(flowON)
  DOPT(doseON)
#undef DOPT
  return pt;
}
void vbl::ReadInParameters::assign(const boost::property_tree::ptree &pt)
{
  #define DOPT(name) boost::property_tree::get_from_ptree(name, #name, pt);
  DOPT(sim_type)
  DOPT(run)
  DOPT(dt)
  DOPT(dt_sm)
  DOPT(t)
  DOPT(t_ini)
  DOPT(treal)
  DOPT(tmax)
  DOPT(tsm_start)
  DOPT(tsm_stop)
  DOPT(slow_motion)
  DOPT(t_CPU_max)
  DOPT(nstep)
  DOPT(nstep_start)
  DOPT(nmax)
  DOPT(idum)
  DOPT(nprint)
  DOPT(nscreen)
  DOPT(nconfiguration)
  DOPT(eps)
  DOPT(delta_vmax)
  DOPT(ncells)
  DOPT(alive)
  DOPT(ntypes)
  DOPT(flowON)
  DOPT(doseON)
  #undef DOPT
}
boost::property_tree::ptree CellsSystem::as_ptree() const
{
  boost::property_tree::ptree big_pt;
  boost::format my_string_template("type_%i");
  for(int k = 0; k<params.ntypes; k++)
  {
    std::string current_type_name = boost::str(my_string_template %  k);
    big_pt.put_child(current_type_name, CellTypeVector[k].as_ptree());
  }
  big_pt.put_child("Environment", Env.as_ptree());
  big_pt.put_child("Environment_0", Env_0.as_ptree());
  
  big_pt.put_child("flowSignal" , flowSignal.as_ptree());
  big_pt.put_child("dose_rateSignal", dose_rateSignal.as_ptree());
  return big_pt;
}

void CellsSystem::assign(const boost::property_tree::ptree &big_pt)
{
  Env.assign(big_pt.get_child("Environment"));
  Env_0.assign(big_pt.get_child("Environment_0"));
  flowSignal.assign(big_pt.get_child("flowSignal"));
  dose_rateSignal.assign(big_pt.get_child("dose_rateSignal"));
  
  boost::format my_string_template("type_%i");
  for(int k = 0; k<big_pt.get<int>("ntypes"); k++)
  {
    std::string current_type_name = boost::str(my_string_template %  k);
    
    //big_pt.put_child(current_type_name, CellTypeVector[k].as_ptree());
  }
}
std::vector<unsigned long> CellsSystem::get_CellTypeIndexVector()
{
  std::cout<< " resized CellTypeIndexVector to: " << params.ncells << std::endl;
  CellTypeIndexVector.resize(params.ncells);
  for(unsigned long n = 0; n<params.ncells; n++)
  {
		unsigned long k = (unsigned long)(type[n]-&CellTypeVector[0]);
    CellTypeIndexVector[n] = k;
		//stream.write( (char*)(&k) , sizeof( unsigned long int ) );
		// cout << k << endl;
  }
  return CellTypeIndexVector;
}
void CellsSystem::set_CellTypeFromIndexVector(std::vector<unsigned long> &cellIndexVector)
{
  type.resize(cellIndexVector.size());
  for(int i=0;i<cellIndexVector.size();i++)
  {
    type[i]=&(CellTypeVector[i]);
  }
}
void CellsSystem::set_CellPhaseFromIntVector(std::vector<int> &int_buffer)
{
  for(int i=0;i<int_buffer.size(); i++)
  {
    Set_phase(i, CellPhase(i));
  }
}
std::vector<int> CellsSystem::Get_phase_int()
{
  std::vector<int> buffer;
  buffer.resize(phase.size());
  for(int i = 0; i< phase.size(); i++)
  {
    buffer[i] = (int) phase[i];
  }
  return buffer;
}

#include "CellsSystem-A.cpp"
#include "CellsSystem-B.cpp"
#include "CellsSystem-C.cpp"
#include "CellsSystem-D.cpp"
#include "CellsSystem-E.cpp"
#include "CellsSystem-F.cpp"

