/*
 *  CellsSystem-B.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */


// se questo parametro e' true allora vengono stampati anche messaggi d'errore aggiuntivi
#define	EXTENDED_ERRORLOG	false
#include "sim.h"
#include "InputFromFile.h"
#include "CellType.h"
#include "Environment.h"
#include "EnvironmentalSignals.h"
#include "BloodVessel.h"
#include "Utilities.h"
//#include "geometry.h"

#include "CellsSystem.h"

using namespace vbl;
// scommentare il seguente statement per includere anche i messaggi d'errore per il metodo della secante
// #define SECANT_IN_ERRORLOG


//
// metodo che gestisce la diffusione nel network di cellule
//
// i nuovi valori vengono calcolati per mezzo di equazioni nonlineari che provengono da una soluzione generale con il metodo BDF 
// (in pratica Eulero implicito); le equazioni non lineari vengono a loro volta risolte per mezzo di una versione del metodo di Newton
// 
void CellsSystem::Diff()
{

	
	double rATPprod=0; // rate di produzione e consumo totale 
	double rATPcons=0;

	// assegnazioni preliminari fondamentali per il loop del calcolo di diffusione e metabolismo (vengono fatte comunque, anche per le cellule morte)

/* TF
 * I hope we do not do error ;-)
 * anyway this broke often when running in large number of threads.
 * Therefore I comment it out.
 */
// #pragma omp parallel for
// 	for(unsigned long n=0; n<ncells; n++)						
//   {
// 		
// 		// controllo preliminare di consistenza
// 		if( phase[n] != dead ) 
// 		{
// 			int code = CheckMVA(n);
// 			if(code < 0) errorlog_file << "Errore " << code << " all'inizio di CellsSystem::Diff nel controllo di consistenza per la cellula " << n << "\n" << std::endl;
// 		}
// 	}

	// volume
	volumeOld = volumeNew = volume;
	volume_extraOld = volume_extraNew = volume_extra;
	
	// mitocondri
	MitOld = MitNew = M;
	
	// pH
	pHiOld = pHiNew = pHi;
	pHOld = pHNew = pH;
			
	// glucosio
	mGinOld = mGinNew = G;
	mGextOld = mGextNew = G_extra;
	
	// G6P
	mG6POld = mG6PNew = G6P;
	
	// O2
	mO2Old = mO2New = O2;
	
	// Store
	StoreOld = StoreNew = store;
	
	// glutammina
	mAinOld = mAinNew = A;
	mAextOld = mAextNew = A_extra;
	
	// AcL
	mAcLinOld = mAcLinNew = AcL;
	mAcLextOld = mAcLextNew = AcL_extra;
	
	// ATP
	ATPpOld = ATPpNew = ATPp;
	
	// protein
	proteinNew = protein;
	
	// delta_protein (inizialmente la variazione e' nulla)
	delta_protein.assign(ncells,0.);
	
	// pRb
	pRbNew = pRb;
	
	// ConcS
	ConcSNew = ConcS;
	
	// DNA 
	DNANew = DNA;
	
	


	
	// valori ambientali
	
	double volume_envOld = Env.GetEnvironmentvolume();		// volume libero dell'ambiente
	double volume_envNew = volume_envOld;
	
	double mG_envOld = Env.GetEnvironmentG();				// massa di glucosio nell'ambiente
	double mG_envNew = mG_envOld;
	
	double mO2_envOld = Env.GetEnvironmentO2();			// massa di O2 nell'ambiente
	double mO2_envNew = mO2_envOld;
	
	double mA_envOld = Env.GetEnvironmentA();				// massa di glutammina nell'ambiente
	double mA_envNew = mA_envOld;
	
	double mAcL_envOld = Env.GetEnvironmentAcL();			// massa di acido lattico nell'ambiente
	double mAcL_envNew = mAcL_envOld;
	

	// * loop per la soluzione delle equazioni nonlineari *
	
	bool isOK;													// variabile che ferma il loop
	ncalls++;													// update del numero di volte che si inizializza il loop
  min_nrepeats = 10.*pow(ncells,1/3.);                        // numero minimo di ripetizioni dell'algoritmo di diffusione
  if(min_nrepeats > MAXREPEATS/2 ) min_nrepeats = MAXREPEATS/2; 
	nrepeats=0;													// inizializzazione del numero di ripetizioni del loop
	
	double prec=eps;										// precisione ottenuta

	do// inizio del loop
  {
    isOK = true;
    
    // prima si calcolano i nuovi valori ambientali 
    
    double s0_env = 0.;									// somme che servono per il calcolo della diffusione nell'ambiente
    double sG_env = 0.;
    double sO2_env = 0.;
    double sA_env = 0.;
    double sAcL_env = 0.;
    double cell_volume = 0.;								// volume totale delle cellule 

#ifdef W_timing
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //begin loop 1
#endif
//#pragma omp parallel for ordered schedule(dynamic)
#pragma omp parallel for
    for(unsigned long n=0; n<ncells; n++)						// loop sulle cellule per il calcolo delle somme
    {
      // volume totale occupato dalle cellule
      cell_volume += volumeOld[n];
      
      // le somme per il calcolo della diffusione nell'ambiente sono non-zero solo per cellule che appartengono all'alpha shape
      if( isonAS[n] )
      {
        s0_env += g_env[n];
        sG_env += g_env[n]*(mGextOld[n]/volume_extraOld[n]);
        sO2_env += g_env[n]*(mO2Old[n]/volumeOld[n]);
        sA_env += g_env[n]*(mAextOld[n]/volume_extraOld[n]);
        sAcL_env += g_env[n]*(mAcLextOld[n]/volume_extraOld[n]);
      }
        
    }
#ifdef W_timing
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    myTiming.time_diff_loop = end-begin;
    myTiming.diff_loop_1 = myTiming.diff_loop_1 + myTiming.time_diff_loop.count();
#endif
	
    volume_envNew = Env_0.GetEnvironmentvolume() - cell_volume;

    if( flowON && ready2start )				// se si lavora con il bioreattore allora le equazioni per l'ambiente devono tenere conto di termini in piu' ... 
    {
      mG_envNew = ( Env.GetEnvironmentG() + dt*( Diff_W_G*sG_env + Env_0.GetEnvironmentG()/Env.GetEnvironmentvolume() * flowSignal.SignalValue(treal) ) )\
        /( 1. + dt*(Diff_W_G*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
      
      // la concentrazione di ossigeno nell'ambiente varia solo se oxygenflowON è acceso
          if( oxygenflowON )
              mO2_envNew = ( Env.GetEnvironmentO2() + dt*( Diff_W_O2*sO2_env + Env_0.GetEnvironmentO2()/Env.GetEnvironmentvolume() * flowSignal.SignalValue(treal) ) )\
              /( 1. + dt*(Diff_W_O2*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
          else
              mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;
      
      mA_envNew = ( Env.GetEnvironmentA() + dt*( Diff_W_A*sA_env + Env_0.GetEnvironmentA()/Env.GetEnvironmentvolume() * flowSignal.SignalValue(treal) ) )\
        /( 1. + dt*(Diff_W_A*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
      
      // l'acido lattico nella soluzione nutritiva pulita e' sempre zero
      mAcL_envNew = ( Env.GetEnvironmentAcL() + dt*( Diff_W_AcL*sAcL_env ) )\
        /( 1. + dt*(Diff_W_AcL*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
    }
    else if( !flowON && ready2start )		// ... altrimenti si utilizzano le equazioni senza i termini aggiuntivi (questo e' vero anche per l'inizializzazione)
    {
      mG_envNew = ( Env.GetEnvironmentG() + dt*( Diff_W_G*sG_env ) )\
        /( 1. + dt*(Diff_W_G*s0_env)/volume_envOld );
      
          if( oxygenflowON )
              mO2_envNew = ( Env.GetEnvironmentO2() + dt*( Diff_W_O2*sO2_env ) )\
              /( 1. + dt*(Diff_W_O2*s0_env)/volume_envOld );
          else
              mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;
      
      mA_envNew = ( Env.GetEnvironmentA() + dt*( Diff_W_A*sA_env ) )\
        /( 1. + dt*(Diff_W_A*s0_env)/volume_envOld );
      
      // l'acido lattico nella soluzione nutritiva pulita e' sempre zero
      mAcL_envNew = ( Env.GetEnvironmentAcL() + dt*( Diff_W_AcL*sAcL_env ) )\
        /( 1. + dt*(Diff_W_AcL*s0_env)/volume_envOld );
    }
    else if ( !ready2start )				// infine, durante la fase di inizializzazione si tengono fisse le concentrazioni ambientali
    {
      mG_envNew = (Env_0.GetEnvironmentG()/Env_0.GetEnvironmentvolume()) * volume_envOld;
      
      mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;
      
      mA_envNew = (Env_0.GetEnvironmentA()/Env_0.GetEnvironmentvolume()) * volume_envOld;
      
      mAcL_envNew = (Env_0.GetEnvironmentAcL()/Env_0.GetEnvironmentvolume()) * volume_envOld;
    }

// P.S. ... note that in the current implementation the bioreactor works only by introducing clean culture liquid, equal to the initial one
// if you want to introduce different doses (for example higher glutamine values than the initial value), then you must modify this part
// also note that the starting time (phase 0) of the pulses that vary over time refers to the initialization time


// Here we move to the calculation of the new values for cells and extracellular spaces

    O2Flow = 0.;	// inizializzazione della somma per il calcolo del rate di assorbimento dell'ossigeno

#ifdef W_timing
    begin = std::chrono::steady_clock::now();
    // begin loop 2
#endif
    double time_run_bico = 0.0;
    
#pragma omp parallel for
    for(unsigned long n=0; n<ncells; n++)// loop on the cells
    {
      // coefficienti (queste assegnazioni non sono strettamente necessarie, ma servono a rendere il codice piu' leggibile)
      // get local thread variables!
      double vmax2 = type[n]->Get_VMAX_2();
      double vmax22 = type[n]->Get_VMAX_22();
      double vmaxP = type[n]->Get_VMAX_P();
      double vmaxP_A = type[n]->Get_VMAX_P_A();
      double vmaxP_ATP = type[n]->Get_VMAX_P_ATP();
      double vmaxDNA = type[n]->Get_VMAX_DNA();
      double v_WORK = type[n]->Get_v_WORK();
      double C1 = type[n]->Get_C1();		
      double C2 = type[n]->Get_C2();		
      double Vmin = type[n]->Get_Vmin();		
      double vmaxDNA_A = type[n]->Get_VMAX_DNA_A();
      double vmaxDNA_ATP = type[n]->Get_VMAX_DNA_ATP();
      double vmaxM = type[n]->Get_VMAX_M();
      double vmaxM_A = type[n]->Get_VMAX_M_A();
      double vmaxM_ATP = type[n]->Get_VMAX_M_ATP();
      double Km1 = type[n]->Get_Km1();
      double Km2 = type[n]->Get_Km2();
      double Km22 = type[n]->Get_Km22();
      double Kmc = type[n]->Get_Kmc();
      double Kmd = type[n]->Get_Kmd();
      double Ka = type[n]->Get_Ka();
      double KmA = type[n]->Get_KmA();
      double KmAL = type[n]->Get_KmAL();
      double Kmp = type[n]->Get_Kmp();
      double KmDNA = type[n]->Get_KmDNA();
      double KmM = type[n]->Get_KmM();
      double ATPSt = type[n]->Get_ATPSt();

      double coeffg1 = type[n]->Get_coeffg1();
      double coeffg2 = type[n]->Get_coeffg2();
      double coeffg3 = type[n]->Get_coeffg3();
      double coeffr1 = type[n]->Get_coeffr1();
      
      
      // quantita' derivate (attenzione, l'ordine e' importante)
      double tpH = 0.5 * ( 1. + tanh( type[n]->Get_tph_slope()*(pHiOld[n] - type[n]->Get_tph_thr()) ) );
      double tp11 = 0.5 * ( 1. + tanh( type[n]->Get_tp11_slope()*(pHiOld[n] - type[n]->Get_tp11_thr()) ) );
      double a2c = 0.5 * ( 1. + tanh( type[n]->Get_a2c_slope()*(pHOld[n] - type[n]->Get_a2c_thr()) ) );
      double c2a = 0.5 * ( 1. + tanh( type[n]->Get_c2a_slope()*(pHiOld[n] - type[n]->Get_c2a_thr()) ) );
      double a2cA = 0.5 * ( 1. + tanh( type[n]->Get_a2cA_slope()*(pHOld[n] - type[n]->Get_a2cA_thr()) ) );
      double c2aA = 0.5 * ( 1. + tanh( type[n]->Get_c2aA_slope()*(pHiOld[n] - type[n]->Get_c2aA_thr()) ) );
      double a2cAcL = 2. - tanh( type[n]->Get_a2cAcL_slope()*( pHOld[n]-type[n]->Get_a2cAcL_thr() ) );
      double c2aAcL = 2. - tanh( type[n]->Get_c2aAcL_slope()*( pHiOld[n]-type[n]->Get_c2aAcL_thr() ) );
      
      // derivata di a2cAcL rispetto la massa di acido lattico nel volume extracellulare: e' l'unica derivata che si deve calcolare al momento
      // ma questo dovrebbe cambiare nel momento in cui si inserira' il nuovo modello di acidita' cellulare
      double da2cAcL = (1./(volume_extraOld[n]*BufCapEnv))/pow(cosh( type[n]->Get_a2cAcL_slope()*( pHOld[n]-type[n]->Get_a2cAcL_thr() ) ),2);
      
      
      if( phase[n] != dead )
        pHiNew[n] = pHi_STANDARD;	// *** al momento il pH interno e' fisso ***
                      // *** questo significa anche che le derivate rispetto pHi sono nulle !!! ***
      else
        pHiNew[n] = pHOld[n];		// se la cellula è morta ha un pH interno uguale a quello esterno ... 
      
      pHNew[n] = 7.5443-( mAcLextOld[n]/volume_extraOld[n])/BufCapEnv;
      
      
      // update di raggio e superficie
      r[n] = pow(3.*volumeOld[n]/(4.*PI), (double)1./3.);
      surface[n] = 4.*PI*r[n]*r[n];

      
      if( phase[n] != dead )
      {
        // if the cell is alive the volume changes
        volumeNew[n] = Vmin * (1+DNANew[n]) + C2*MitOld[n] + C1*ATPpOld[n];		
        volume_extraNew[n] = surface[n]*type[n]->Get_extvolume_thickness()*type[n]->Get_extvolume_fraction();
      }
      else
      {
        // if the cell is dead the volume remains fixed (it is reduced only at the end of this method)
        volumeNew[n] = volumeOld[n];
        volume_extraNew[n] = volume_extraOld[n];
      }
      
      ATP_St[n] = ATPSt;
      double dSensO2_O2 = 0.;	// *** REV *** does it make sense to include these two vars in CellsSystem ???
      double SensATP = 0.;
      
      if( phase[n] != dead )
      {
        double buffer_fast = volumeOld[n]*type[n]->Get_KmO2();
        double buffer_fast2 = StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n]);
        double buffer_fast3 = mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]);
        SensO2[n] = mO2Old[n]/(buffer_fast + mO2Old[n]);
        dSensO2_O2 = buffer_fast/pow(buffer_fast + mO2Old[n],2);
        ATP_Ox[n] = 30.*(PM_ATP/PM_G)*SensO2[n] * ( coeffg2 + coeffg3*buffer_fast2 ) * mG6POld[n];
        SensATP = 0.5*( 1 - tanh( 100.*(ATP_Ox[n]/ATP_St[n]-1.) ) );
        ATP_NOx[n] = 2.*(PM_ATP/PM_G)*tpH*( coeffg1*mG6POld[n] + coeffr1*buffer_fast2 );
        ATP2[n] = SensO2[n]*SensATP*(ATP_St[n] - ATP_Ox[n])*buffer_fast2 ;
        ATP3[n] = (14./5.)*SensO2[n]*SensATP*(ATP_St[n] - ATP_Ox[n]) * buffer_fast3;
        ConsATP[n] = 2.*(PM_ATP/PM_G) * tp11*tpH * coeffr1 * (1.-SensATP) * buffer_fast3;
      }
      else
        ATP_Ox[n] = ATP_NOx[n] = ATP2[n] = ConsATP[n] = 0.;
      
      double vp, dvp_A, dvp_ATPp;					// produzione di proteine senza il termine vmax
      double buffer_ATPpOld_n_times_mAinOld_n_ = ATPpOld[n]*mAinOld[n];
      // la produzione di proteine e' bloccata in fase G0
      if(phase[n] != G0_phase && phase[n] != dead)	
      {
        vp = buffer_ATPpOld_n_times_mAinOld_n_/(pow(volumeOld[n],2)*Kmp+buffer_ATPpOld_n_times_mAinOld_n_);
        dvp_A = ATPpOld[n]*pow(volumeOld[n],2)*Kmp/pow(pow(volumeOld[n],2)*Kmp+buffer_ATPpOld_n_times_mAinOld_n_,2);
        dvp_ATPp = pow(volumeOld[n],2)*Kmp*mAinOld[n]/pow(pow(volumeOld[n],2)*Kmp+buffer_ATPpOld_n_times_mAinOld_n_,2);
      }
      else
        vp= dvp_A = dvp_ATPp = 0.;
      
      double vDNA, dvDNA_A, dvDNA_ATPp;			// produzione di DNA senza il termine vmax
      if(phase[n] == S_phase)	// la produzione di DNA procede solo in fase S
      {
        vDNA = buffer_ATPpOld_n_times_mAinOld_n_/(pow(volumeOld[n],2)*KmDNA+buffer_ATPpOld_n_times_mAinOld_n_);
        dvDNA_A = ATPpOld[n]*pow(volumeOld[n],2)*KmDNA/pow(pow(volumeOld[n],2)*KmDNA+buffer_ATPpOld_n_times_mAinOld_n_,2);
        dvDNA_ATPp = pow(volumeOld[n],2)*KmDNA*mAinOld[n]/pow(pow(volumeOld[n],2)*KmDNA+buffer_ATPpOld_n_times_mAinOld_n_,2);
      }
      else
        vDNA = dvDNA_A = dvDNA_ATPp = 0.;
      DNA_rate[n] = vmaxDNA*vDNA;
      
      // produzione di mtDNA: termine comune a mitocondri, glutammina e ATPp senza le vmax, e sue derivate rispetto A e ATPp
      double vM, dvM_A, dvM_ATPp;
      if(phase[n] != dead)
      {
        vM = buffer_ATPpOld_n_times_mAinOld_n_/(pow(volumeOld[n],2)*KmM+buffer_ATPpOld_n_times_mAinOld_n_);
        dvM_A = ATPpOld[n]*pow(volumeOld[n],2)*KmM/pow(pow(volumeOld[n],2)*KmM+buffer_ATPpOld_n_times_mAinOld_n_,2);
        dvM_ATPp = pow(volumeOld[n],2)*KmM*mAinOld[n]/pow(pow(volumeOld[n],2)*KmM+buffer_ATPpOld_n_times_mAinOld_n_,2);
      }
      else
        vM = dvM_A = dvM_ATPp = 0.;


      // qui si definisce il valore dell' O2-dependent glucose-transport efficiency h
      // ma solo se la cellula e' viva ... 
      if( phase[n] != dead )
      {
        double buffer_quick = 1.-mO2Old[n]/(volumeOld[n]*O2st);
        double buffer_quick2 = tanh(100.*(buffer_quick));
        h[n] = 0.5*(1.3*(buffer_quick)+1.)*(1.+buffer_quick2) + 0.5*(1-buffer_quick2);
      }
      
      double v1max = type[n]->Get_VMAX_1()*surface[n]*h[n];
      double vmaxA = type[n]->Get_VMAX_A()*surface[n];
      double vmaxAcL = type[n]->Get_VmaxAL0() * surface[n];
      
  // ***** contribution of blood vessels, start init

      // blood vessel concentrations (cell is in contact with at most one blood vessel)
      // note that blood vessel concentrations are kept constant in this version, unlike
      // those in the environment

      double rhoG_bv = 0.;
      double rhoO2_bv = 0.;
      double rhoA_bv = 0.;
      double rhoAcL_bv = 0.;

      if(isonBV[n])//n for until numerOfCells
      {
        BloodVessel BV = BloodVesselVector[isonBV[n]-1];	// extract blood vessel
        //vbl::BloodVessel BV = bloodVesselMap[isonBV[n]-1];
          
        rhoG_bv = BV.GetBloodVesselG();			// glucose concentration in BV	
        rhoO2_bv = 0.5*(BV.GetBloodVesselO2start()+BV.GetBloodVesselO2end());			
                            // oxygen concentration in BV
        rhoA_bv = BV.GetBloodVesselA();			// other nutrients concentration in BV
        rhoAcL_bv = BV.GetBloodVesselAcL();		// lactate concentration in BV 
        
        // cout << rhoG_bv << "\t" << rhoO2_bv << "\t" << rhoA_bv << "\t" << rhoAcL_bv << endl;
      }

  // ***** contribution of blood vessels, end init
      
      
      // sums
      double s0 = 0.;
      double sG = 0.;
      double sO2 = 0.;
      double sA = 0.;
      double sAcL = 0.;
      
      // calculation of the sums for the diffusion in the cluster of cells
      // WARNING nested loop --> what values does neigh[n] have?
      for(int nv=0; nv<neigh[n]; nv++)		
      {
        int name = vneigh[n][nv];
        s0 += gnk[n][nv];
        sG += gnk[n][nv]*mGextOld[name]/volume_extraOld[name];			
        sO2 += gnk[n][nv]*mO2Old[name]/volumeOld[name];			
        sA += gnk[n][nv]*mAextOld[name]/volume_extraOld[name];		
        sAcL += gnk[n][nv]*mAcLextOld[name]/volume_extraOld[name];			
      }

      // *** glucose
      
      double M_Gnow, DM_Gnow, T_Gnow, DT_Gnow_in, DT_Gnow_ext;	// valori che servono al calcolo con il metodo di Newton
      double M_Gnow_0, T_Gnow_0in, T_Gnow_0ext;					// valori aggiuntivi per il calcolo con il metodo della secante
    
      //if( true ) // environment for local speed up variables!
      {
        double buffer_1 = vmax2*mGinOld[n]*mGinOld[n];
        double buffer_2 = vmax22*mGinOld[n]*mGinOld[n];
        double buffer_3 = volumeOld[n]*Km2+mGinOld[n];
        double buffer_4 = volumeOld[n]*Ka+mGinOld[n];
        double buffer_5 = volumeOld[n]*Km22+mGinOld[n];
        // funzione metabolica	
        M_Gnow = ((int)(phase[n] != dead))* \
                ( -buffer_1/((buffer_3)*(buffer_4)) \
                - buffer_2/((buffer_5)*(buffer_4)) );
        // derivata funzione metabolica rispetto mGin calcolata in mGinOld
        DM_Gnow = ((int)(phase[n] != dead))* \
                ( -2*vmax2*mGinOld[n]/((buffer_3)*(buffer_4)) \
                + buffer_1/(pow(buffer_3,2)*(buffer_4)) \
                + buffer_1/((buffer_3)*pow(buffer_4,2)) \
                - 2*vmax22*mGinOld[n]/((buffer_5)*(buffer_4)) \
                + buffer_2/(pow(buffer_5,2)*(buffer_4)) \
                + buffer_2/((buffer_5)*pow(buffer_4,2)) );
            
        M_Gnow_0 = 0.;	// valore della funzione metabolica in mGin = 0
                
        // funzione di trasporto
        double buffer_6 = volume_extraOld[n]*Km1+mGextOld[n];
        double buffer_7 = volumeOld[n]*Km1+mGinOld[n];
        T_Gnow = a2c*v1max*mGextOld[n]/(buffer_6) - c2a*v1max*mGinOld[n]/(buffer_7);
        // derivata della funzione di trasporto rispetto mGin
        DT_Gnow_in = -c2a*v1max*volumeOld[n]*Km1/pow(buffer_7,2);
        // derivata della funzione di trasporto rispetto mGext
        DT_Gnow_ext = a2c*v1max*volume_extraOld[n]*Km1/pow(buffer_6,2);
        
        // valori speciali della funzione di trasporto
        T_Gnow_0in = a2c*v1max*mGextOld[n]/(buffer_6);	// mGin = 0
        T_Gnow_0ext = - c2a*v1max*mGinOld[n]/(buffer_7);			// mGext = 0
      }
                          
      // 1. calculation of the new value of mGin
      mGinNew[n] = mGinOld[n] - ( mGinOld[n] - ( G[n] + dt*( M_Gnow + T_Gnow ) ) )/( 1. - dt*( DM_Gnow + DT_Gnow_in ) );

      if( mGinNew[n] < 0 )	// if you get to this point it means that Newton's method has problems ...
      {
        double fhi = mGinOld[n] - ( G[n] + dt*( M_Gnow + T_Gnow ) );	// f valutata in mGinOld[n]
        double flo = - ( G[n] + dt*( M_Gnow_0 + T_Gnow_0in ) );		// f valutata in mGin = 0
        double xsol = mGinOld[n]*flo/(flo-fhi);
        mGinNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << endl;
        errorlog_file << scientific << "mGinOld["<< n << "] = " << mGinOld[n] << "; mGinNew["<< n << "] = " << mGinNew[n] << "\n" <<std::endl;
  #endif
        // TF. I think this is not necessary
//       }
//       if( mGinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mGinNew["<< n << "] = " << std::scientific << mGinNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mGinNew[n] = mGinOld[n] = " << mGinOld[n] << "\n" << std::endl;
        
        mGinNew[n] = mGinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }



      // 2. calculation of the new value of mGext *** MODIFIED FOR BLOOD VESSELS ***
      double cG = G_extra[n] + dt*Diff_ES_G*sG + dt*Diff_Env_G*(mG_envOld/volume_envOld)*g_env[n] + dt*Diff_BV_G*rhoG_bv*g_bv[n];	// calcolo della parte fissa
      //if( true )//local speed up variables
      {
        double buffer = 1.+Diff_ES_G*dt*s0/volume_extraOld[n]+Diff_Env_G*dt*g_env[n]/volume_extraOld[n]+Diff_BV_G*dt*g_bv[n]/volume_extraOld[n];
        mGextNew[n] = mGextOld[n] - ( mGextOld[n]*( buffer) + dt*T_Gnow - cG )\
            /( buffer + dt*DT_Gnow_ext );
      }
      if(mGextNew[n] < 0)	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = ( mGextOld[n]*(1+Diff_ES_G*dt*s0/volume_extraOld[n]+Diff_Env_G*dt*g_env[n]/volume_extraOld[n]+Diff_BV_G*dt*g_bv[n]/volume_extraOld[n]) + dt*T_Gnow - cG );
        double flo = ( dt*T_Gnow_0ext - cG );
        double xsol = mGextOld[n]*flo/(flo-fhi);
        mGextNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << endl;
        errorlog_file << scientific << "mGextOld["<< n << "] = " << mGextOld[n] << "; mGextNew["<< n << "] = " << mGextNew[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mGextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mGextNew["<< n << "] = " << std::scientific << mGextNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mGextNew[n] = mGextOld[n] = " << mGextOld[n] << "\n" << std::endl;
        mGextNew[n] = mGextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }
      
      // storage of consumption and transport
      GConsRate[n] = (double)M_Gnow;
      GAbsRate[n] = (double)T_Gnow;
      
      // *** END OF SECTION MODIFIED FOR BLOOD VESSELS ***



      
      // *** G6P (exact solution of the iterative equation derived from the Euler method)
      if( phase[n] != dead )
        mG6PNew[n] = (mG6POld[n] - dt*M_Gnow)/( 1. + dt*( tpH*coeffg1 + coeffg2*SensO2[n] + coeffg3) );
      else
        mG6PNew[n] = mG6POld[n];



      // *** O2      *** MODIFIED FOR BLOOD VESSELS *** 
          
      double M_O2now, DM_O2now;								// valori che servono al calcolo con il metodo di Newton
      double M_O2now_0;										// valore aggiuntivo per il calcolo con il metodo della secante
      
      if( phase[n] != dead )
      {
        //if( true )//local speed up environment
        {
          double fac1 = PM_O2/PM_G;
          double fac2 = PM_G/PM_ATP;
          double fac3 = ATP_St[n] - ATP_Ox[n];
          double buffer = StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n]);
          double buffer2 = mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]);
          M_O2now = -6.*(fac1) * SensO2[n] * ( coeffg2*mG6POld[n] \
            + coeffg3*mG6POld[n]*(buffer) \
            + SensATP*0.033333*(fac2)*( fac3 ) * (buffer) \
            + 0.1*SensATP*(fac2)*( fac3 ) * (buffer2)  );
          DM_O2now = -6.*(fac1) * dSensO2_O2 * ( coeffg2*mG6POld[n] \
            + coeffg3*mG6POld[n]*(buffer) \
            + SensATP*0.033333*(fac2)*( fac3 ) * (buffer) \
            + 0.1*SensATP*(fac2)*( fac3 ) * (buffer2)  );
            
          M_O2now_0 = 0.;	// valore della funzione metabolica in mO2 = 0
        }
      }
      else
      {
        M_O2now = DM_O2now = 0.;	// se la cellula e' morta non c'e' alcun consumo metabolico
        M_O2now_0 = 0.;
      }
        
      O2Flow += M_O2now;
        
      // calcolo del nuovo valore
      double cO2 = O2[n] + dt*Diff_ES_O2*sO2 + dt*Diff_Env_O2*(mO2_envOld/volume_envOld)*g_env[n] + dt*Diff_BV_O2*rhoO2_bv*g_bv[n];	// calcolo della parte fissa
      
      //if( true )// speed up environment
      {
        double buffer = 1+Diff_ES_O2*dt*s0/volumeOld[n]+Diff_Env_O2*dt*g_env[n]/volumeOld[n]+Diff_BV_O2*dt*g_bv[n]/volumeOld[n];
        mO2New[n] = mO2Old[n] - ( mO2Old[n] * ( buffer) - dt*M_O2now - cO2 )\
            /( buffer - dt*DM_O2now);
      }
      
      if( mO2New[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = mO2Old[n] * ( 1+Diff_ES_O2*dt*s0/volumeOld[n]+Diff_Env_O2*dt*g_env[n]/volumeOld[n]+Diff_BV_O2*dt*g_bv[n]/volumeOld[n]) - dt*M_O2now - cO2;
        double flo = - dt*M_O2now_0 - cO2;
        double xsol = mO2Old[n]*flo/(flo-fhi);
        mO2New[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << endl;
        errorlog_file << scientific << "mO2Old["<< n << "] = " << mO2Old[n] << "; mO2New["<< n << "] = " << mO2New[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mO2New[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mO2New["<< n << "] = " << std::scientific << mO2New[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mO2New[n] = mO2Old[n] = " << mO2Old[n] << "\n" << std::endl;
        mO2New[n] = mO2Old[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }
      
      // *** END OF SECTION MODIFIED FOR BLOOD VESSELS ***


      // *** glutammina
      
      double M_Anow, DM_Anow, T_Anow, DT_Anow_in, DT_Anow_ext;	// valori che servono al calcolo con il metodo di Newton
      double M_Anow_0, T_Anow_0in, T_Anow_0ext;					// valori aggiuntivi per il calcolo con il metodo della secante
      
      // metabolic function
      //if( true )//local speed up
      {
        double buffer = tp11*tpH*coeffr1*(1-SensATP);
        double buffer2 = 0.1*SensO2[n]*SensATP*(PM_G/PM_ATP)*(ATP_St[n]-ATP_Ox[n]);
        double buffer3 = volumeOld[n]*Kmd+mAinOld[n];
        double buffer4 = volume_extraOld[n]*KmA+mAextOld[n];
        double buffer5 = volumeOld[n]*KmA+mAinOld[n];
        M_Anow = ((int)(phase[n] != dead))* \
            ( -( buffer + buffer2 ) * (mAinOld[n]/(buffer3)) \
            - vmaxP_A*vp - vmaxDNA_A*vDNA - vmaxM_A*vM );
        // derivata funzione metabolica rispetto mAin calcolata in mAinOld
        DM_Anow = ((int)(phase[n] != dead))* \
            ( -( buffer + buffer2 ) * (volumeOld[n]*Kmd/pow(buffer3,2)) \
            - vmaxP_A*dvp_A - vmaxDNA_A*dvDNA_A - vmaxM_A*dvM_A ) ;
                
        M_Anow_0 = 0;	// valore della funzione metabolica in mAin = 0
        
        // funzione di trasporto
        T_Anow = a2cA*vmaxA * mAextOld[n]/(buffer4) - c2aA*vmaxA * mAinOld[n]/(buffer5);
        // derivata della funzione di trasporto rispetto mAin
        DT_Anow_in = - c2aA*vmaxA * volumeOld[n]*KmA/pow(buffer5,2);
        // derivata della funzione di trasporto rispetto mAext
        DT_Anow_ext = a2cA*vmaxA * volume_extraOld[n]*KmA/pow(buffer4,2);

        // valori speciali della funzione di trasporto
        T_Anow_0in = a2cA*vmaxA * mAextOld[n]/(buffer4);	// mAin = 0
        T_Anow_0ext = - c2aA*vmaxA * mAinOld[n]/(buffer5);		// mAext = 0
        
        // 1. calcolo del nuovo valore di mAin
        mAinNew[n] = mAinOld[n] - ( mAinOld[n] - ( A[n] + dt*( M_Anow + T_Anow ) ) )/( 1. - dt*( DM_Anow + DT_Anow_in ) );
      }
      
      if( mAinNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = mAinOld[n] - ( A[n] + dt*( M_Anow + T_Anow ) );
        double flo = - ( A[n] + dt*( M_Anow_0 + T_Anow_0in ) );
        double xsol = mAinOld[n]*flo/(flo-fhi);
        mAinNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << std::endl;
        errorlog_file << scientific << "mAinOld["<< n << "] = " << mAinOld[n] << "; mAinNew["<< n << "] = " << mAinNew[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mAinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mAinNew["<< n << "] = " << std::scientific << mAinNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mAinNew[n] = mAinOld[n] = " << mAinOld[n] << "\n" << std::endl;
        mAinNew[n] = mAinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }
      
      
      
      // 2. calcolo del nuovo valore di mAext      *** MODIFIED FOR BLOOD VESSELS ***
      double cA = A_extra[n] + dt*Diff_ES_A*sA + dt*Diff_Env_A*(mA_envOld/volume_envOld)*g_env[n] + dt*Diff_BV_A*rhoA_bv*g_bv[n];	// calcolo della parte fissa
      
      //if( true )// speed up
      {
        double buffer = (1.+Diff_ES_A*dt*s0/volume_extraOld[n]+Diff_Env_A*dt*g_env[n]/volume_extraOld[n]+Diff_BV_A*dt*g_bv[n]/volume_extraOld[n]);
        mAextNew[n] = mAextOld[n] - ( mAextOld[n]*buffer + dt*T_Anow - cA )\
            /(buffer + dt*DT_Anow_ext );
      }
      if( mAextNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = ( mAextOld[n]*(1+Diff_ES_A*dt*s0/volume_extraOld[n]+Diff_Env_A*dt*g_env[n]/volume_extraOld[n]+Diff_BV_A*dt*g_bv[n]/volume_extraOld[n]) + dt*T_Anow - cA );
        double flo = ( dt*T_Anow_0ext - cA );
        double xsol = mAextOld[n]*flo/(flo-fhi);
        mAextNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << std::endl;
        errorlog_file << scientific << "mAextOld["<< n << "] = " << mAextOld[n] << "; mAextNew["<< n << "] = " << mAextNew[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mAextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mAextNew["<< n << "] = " << std::scientific << mAextNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mAextNew[n] = mAextOld[n] = " << mAextOld[n] << "\n" << std::endl;
        mAextNew[n] = mAextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }

      // memorizzazione di consumo e trasporto
      AConsRate[n] = (double)M_Anow;
      AAbsRate[n] = (double)T_Anow;
      
      // *** END OF SECTION MODIFIED FOR BLOOD VESSELS ***


      
      
      // *** store (soluzione esatta dell'equazione iterativa derivata dal metodo di Eulero)

      if( phase[n] != dead )
      {
        //if( true )//local speed up
        {
          double buffer2 = (double) (coeffg3*mG6POld[n] + tp11*tpH*coeffr1*(mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]))*(1-SensATP));
          double buffer3 = (double) ( tpH*coeffr1 + coeffg3*mG6POld[n]*SensO2[n] + SensO2[n]*SensATP*0.033333*(PM_G/PM_ATP)*( ATP_St[n] - ATP_Ox[n] ) );
          double Ast = store[n] + dt * ( buffer2 );
          double Bst = -dt * ( buffer3 );
          double Kst = volumeOld[n]*Kmc;

        
          double buffer = Ast+Bst-Kst;
          StoreNew[n] = 0.5*( (buffer) + sqrt( pow(buffer,2) + 4*Ast*Kst ) );
          
          // memorizzazione di riempimento e consumo
          StoreFillRate[n] = buffer2;
          StoreConsRate[n] = (double)(-buffer3 * (StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])));
        }
      }
      else
      {
        StoreNew[n] = StoreOld[n];
        StoreFillRate[n] = StoreConsRate[n] = 0.;
      }

      // *** AcL 
      
      double M_AcLnow, DM_AcLnow, T_AcLnow, DT_AcLnow_in, DT_AcLnow_ext;
      double M_AcLnow_0, T_AcLnow_0in, T_AcLnow_0ext;
      
      //if( true )//speed up buffer
      {
        double buffer = volume_extraOld[n]*KmAL+mAcLextOld[n];
        double buffer2 = volumeOld[n]*KmAL+mAcLinOld[n];
        // funzione metabolica
        M_AcLnow = ((int)(phase[n] != dead))* \
            ( 2.*tpH*( coeffg1*mG6POld[n] + coeffr1*(StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) ) );
            
        // PAY LARGE ATTENTION TO THE DERIVATIVES COMPARED WITH THE COEFFICIENTS THAT ADJUST TRANSPORT ...
        // below ASSUME THAT THE DERIVATIVES RESPECT pHi ARE NULLS because at the moment pHi is kept constant and does not change as a function of mAcLin.
        
        // derivata funzione metabolica rispetto mAcLin calcolata in mAcLinOld			
        DM_AcLnow = 0.;
        
        M_AcLnow_0 = M_AcLnow;	// valore della funzione metabolica in mAcLin = 0: nella formulazione attuale la produzione di AcL non dipende dal contenuto di AcL della cellula
        
        // funzione di trasporto
        T_AcLnow = a2cAcL*vmaxAcL * mAcLextOld[n]/(buffer) - c2aAcL*vmaxAcL * mAcLinOld[n]/(buffer2);
        // derivata della funzione di trasporto rispetto mAcLin
        DT_AcLnow_in = - c2aAcL*vmaxAcL * volumeOld[n]*KmAL/pow(buffer2,2);
        // derivata della funzione di trasporto rispetto mAcLext			
        DT_AcLnow_ext = da2cAcL*vmaxAcL * mAcLextOld[n]/(buffer) + \
                        a2cAcL*vmaxAcL * volume_extraOld[n]*KmAL/pow(buffer,2);
        
        // valori speciali della funzione di trasporto
        T_AcLnow_0in = a2cAcL*vmaxAcL * mAcLextOld[n]/(buffer);	// mAcLin = 0
        T_AcLnow_0ext = - c2aAcL*vmaxAcL * mAcLinOld[n]/(buffer2);		// mAcLext = 0
          

        // 1. calcolo del nuovo valore di mAcLin
        mAcLinNew[n] = mAcLinOld[n] - ( mAcLinOld[n] - ( AcL[n] + dt*( M_AcLnow + T_AcLnow ) ) )/( 1 - dt*( DM_AcLnow + DT_AcLnow_in ) );
      }
      
      if( mAcLinNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = mAcLinOld[n] - ( AcL[n] + dt*( M_AcLnow + T_AcLnow ) );	// f valutata in mAcLinOld[n]
        double flo = - ( AcL[n] + dt*( M_AcLnow_0 + T_AcLnow_0in ) );			// f valutata in mAcLin = 0
        double xsol = mAcLinOld[n]*flo/(flo-fhi);
        mAcLinNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << std::endl;
        errorlog_file << std::scientific << "mAcLinOld["<< n << "] = " << mAcLinOld[n] << "; mAcLinNew["<< n << "] = " << mAcLinNew[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mAcLinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mAcLinNew["<< n << "] = " << std::scientific << mAcLinNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mAcLinNew[n] = mAcLinOld[n] = " << mAcLinOld[n] << "\n" << std::endl;
        mAcLinNew[n] = mAcLinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }
      
      // 2. calcolo del nuovo valore di mAcLext      *** MODIFIED FOR BLOOD VESSELS ***
      double cAcL = AcL_extra[n] + dt*Diff_ES_AcL*sAcL + dt*Diff_Env_AcL*(mAcL_envOld/volume_envOld)*g_env[n] + dt*Diff_BV_AcL*rhoAcL_bv*g_bv[n];	// calcolo della parte fissa
      
      if( true ) // local speed up
      {
        double buffer = Diff_Env_AcL*dt*g_env[n]/volume_extraOld[n];
        double buffer2 = 1.+Diff_ES_AcL*dt*s0/volume_extraOld[n];
        double buffer3 = Diff_BV_AcL*dt*g_bv[n]/volume_extraOld[n];
        mAcLextNew[n] = mAcLextOld[n] - ( mAcLextOld[n]*(buffer2+buffer+buffer3) + dt*T_AcLnow - cAcL )\
          /( buffer2 + buffer + buffer3 + dt*DT_AcLnow_ext );
      }
      if( mAcLextNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
      {
        double fhi = mAcLextOld[n]*(1+Diff_ES_AcL*dt*s0/volume_extraOld[n]+Diff_Env_AcL*dt*g_env[n]/volume_extraOld[n]+Diff_BV_AcL*dt*g_bv[n]/volume_extraOld[n]) + dt*T_AcLnow - cAcL;
        double flo = dt*T_AcLnow_0ext - cAcL;
        double xsol = mAcLextOld[n]*flo/(flo-fhi);
        mAcLextNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
        errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
        errorlog_file << " e' necessario passare al metodo della secante" << std::endl;
        errorlog_file << std::scientific << "mAcLextOld["<< n << "] = " << mAcLextOld[n] << "; mAcLextNew["<< n << "] = " << mAcLextNew[n] << "\n" <<std::endl;
  #endif
//       }
//       if( mAcLextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//       {
        errorlog_file << "ATTENZIONE: al passo " << nstep << " mAcLextNew["<< n << "] = " << std::scientific << mAcLextNew[n] << " < 0" << std::endl;
        errorlog_file << "in questa iterazione si mantiene mAcLextNew[n] = mAcLextOld[n] = " << mAcLextOld[n] << "\n" << std::endl;
        mAcLextNew[n] = mAcLextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
      }

      // memorizzazione di consumo e trasporto
      AcLRate[n] = (double)M_AcLnow;
      AcLOutRate[n] = (double)T_AcLnow;
      
      // *** END OF SECTION MODIFIED FOR BLOOD VESSELS ***




      // *** ATP
      
      if( phase[n] != dead )
      {
        
        double M_ATPnow, DM_ATPnow;								// valori che servono al calcolo con il metodo di Newton
        
        // memorizzazione di produzione e consumo
        ConsATP_1[n] = v_WORK*C1*ATPpOld[n];
        ConsATP_2[n] = vmaxP_ATP*vp;
        ConsATP_3[n] = vmaxDNA_ATP*vDNA;
        ConsATP_5[n] = vmaxM_ATP*vM;
        
        rATPprod = ATP_Ox[n] + ATP_NOx[n] + ATP2[n] + ATP3[n];
        rATPcons = ConsATP[n] + ConsATP_1[n] + ConsATP_2[n] + ConsATP_3[n] + ConsATP_5[n];

        // funzione metabolica
        M_ATPnow = rATPprod - rATPcons;
        // derivata funzione metabolica rispetto ATPp calcolata in ATPpOld
        DM_ATPnow = - v_WORK*C1 - vmaxP_ATP*dvp_ATPp - vmaxDNA_ATP*dvDNA_ATPp - vmaxM_ATP*dvM_ATPp; 

        // calcolo del nuovo valore di ATPp
        ATPpNew[n] = ATPpOld[n] - ( ATPpOld[n] - ATPp[n] - dt*M_ATPnow )/(1.-dt*DM_ATPnow);
        if( ATPpNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
        {
          double fhi = ATPpOld[n] - ATPp[n] - dt*M_ATPnow;
          double flo = - ( ATPp[n] + dt*(ATP_NOx[n]) );
          double xsol = ATPpOld[n]*flo/(flo-fhi);
          ATPpNew[n] = xsol;
  #ifdef SECANT_IN_ERRORLOG
          errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
          errorlog_file << " e' necessario passare al metodo della secante" << endl;
          errorlog_file << std::scientific << "ATPpOld["<< n << "] = " << ATPpOld[n] << "; ATPpNew["<< n << "] = " << ATPpNew[n] << "\n" <<std::endl;
  #endif
//         }
//         if( ATPpNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
//         {
          errorlog_file << "ATTENZIONE: al passo " << nstep << " ATPpNew["<< n << "] = " << std::scientific << ATPpNew[n] << " < 0" << std::endl;
          errorlog_file << "in questa iterazione si mantiene ATPpNew[n] = ATPpOld[n] = " << ATPpOld[n] << "\n" << std::endl;
          errorlog_file << "Altre informazioni sullo stato cellulare attuale: " << std::endl;
          errorlog_file << "eta' cellulare " << age[n] << std::endl;
          errorlog_file << "eta' di fase cellulare " << phase_age[n] << std::endl;
          errorlog_file << "ATP_St " << ATP_St[n] << std::endl;
          errorlog_file << "ATP_Ox " << ATP_Ox[n] << std::endl;
          errorlog_file << "ATP_NOx " << ATP_NOx[n] << std::endl;
          errorlog_file << "ATP2 " << ATP2[n] << std::endl;
          errorlog_file << "ATP3 " << ATP3[n] << std::endl;
          errorlog_file << "ConsATP " << ConsATP[n] << std::endl;
          errorlog_file << "ConsATP_1 " << ConsATP_1[n] << std::endl;
          errorlog_file << "ConsATP_2 " << ConsATP_2[n] << std::endl;
          errorlog_file << "ConsATP_3 " << ConsATP_3[n] << std::endl;
          errorlog_file << "ConsATP_4 " << ConsATP_4[n] << std::endl;
          errorlog_file << "ConsATP_5 " << ConsATP_5[n] << std::endl;
          errorlog_file << "ATPtot " << ATPtot[n] << std::endl;
          errorlog_file << "ATPp " << ATPp[n] << std::endl;
          errorlog_file << "ATPmin " << ATPmin[n] << std::endl;
          errorlog_file << std::endl;
          ATPpNew[n] = ATPpOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
        }
        
        // proliferazione dei mitocondri (la proliferazione dei mitocondri e' bloccata in fase G0)
        if(phase[n] != G0_phase)
          MitNew[n] = M[n] + vmaxM*vM*dt; 
        else
          MitNew[n] = M[n];
      }
      else
      {
        ConsATP_1[n] = 0.;
        ConsATP_2[n] = 0.;
        ConsATP_3[n] = 0.;
        ConsATP_5[n] = 0.;
        ATPpNew[n] = ATPpOld[n];
        MitNew[n] = M[n];
      }
        
      
      // *** proteine
      
      if( phase[n] != dead )
      {
        delta_protein[n] = vmaxP*vp*dt;
        proteinNew[n] = protein[n] + delta_protein[n];
        
        if( phase[n] == G2_phase || phase[n] == M_phase)
          pRbNew[n] = pRb[n] + type[n]->Get_pRb_fraction()*delta_protein[n];
        else
          pRbNew[n] = pRb[n];
          
          
        // calculation of the number of active pRb molecules
        // The approximate calculation of p, of Pk_pRb and of the number of active molecules NpRbk
        //
        // Here we assume a process of phosphorylation of this type:
        // 1. Cyclin D binds to ATP
        // 2. The complex is linked to kinasi
        // 3. all this complexity is to phosphorylate the pRb on the individual sites
        //
        // Let us assume that steps 1 and 2 can proceed (in other words we assume that the necessary ATP is)
        // much higher than the death threshold)
        // What we call ConcCyclinD is actually the concentration of the phosphorylating complex

        double ConcpRb = 1000.*pRbNew[n]/(PM_pRb * volume[n]);				// concentrazione MOLARE (!!!!) della pRb nella cellula
        double ConcCyclinD = 1000.*cyclinD[n]/(PM_cyclinD * volume[n] );	// concentrazione MOLARE (!!!!) della ciclina D nella cellula
        double ConcCyclinE = 1000.*cyclinE[n]/(PM_cyclinE * volume[n] );	// concentrazione MOLARE (!!!!) della ciclina E nella cellula
              // il coeff di conversione 1000 serve per passare da pg e micron alle unità per il calcolo
              // della conc. molare ... 
                  
        int N_pRb = type[n]->Get_N_pRb();		// number of phosphorylation sites of the pRb
        int k_pRb = type[n]->Get_k_pRb();		// number of threshold sites
        double pRb_ONOFFratio = type[n]->Get_pRb_ONOFFratio();
        double buffer = ConcCyclinD+ConcCyclinE;
        double buffer2 = N_pRb*ConcpRb;
        double p_pRb = ( (buffer2 + (buffer) + pRb_ONOFFratio) - \
            sqrt( pow(buffer2 + (buffer) + pRb_ONOFFratio,2) - 4.*buffer2*(buffer) ) )/(2. * buffer2);
        
        double Pk_pRb = 0.;
        // WARNING nested loop --> how much N_pRb are there?
        // NOTE bound by maximal number of neighbours!

        for(int l=k_pRb; l<=N_pRb; l++)
        {
          Pk_pRb += bico(N_pRb,k_pRb) * pow(p_pRb,l)*pow(1.-p_pRb,N_pRb-l);
        }
        
        double ConcE = ConcpRb*Pk_pRb;	// concentrazione MOLARE dell'enzima

        // reazione downstream (MM)
        double buffer3 = ConcS[n] - dt*ConcE*type[n]->Get_k3MM() - type[n]->Get_KmMM();
        ConcSNew[n] = 0.5*( (buffer3) \
          + sqrt( pow(buffer3,2) \
              + 4.*ConcS[n]*type[n]->Get_KmMM() ) );
      }
      else
      { //dead phase
        delta_protein[n] = 0.;
        proteinNew[n] = protein[n];
        pRbNew[n] = pRb[n];
        ConcSNew[n] = ConcS[n];
      }
      
      // *** DNA
      if(phase[n] == S_phase)
        DNANew[n] = DNA[n] + DNA_rate[n]*dt;
      else
        DNANew[n] = DNA[n];
      
    }	// end for(unsigned long n=0; n<ncells; n++)// loop on the cells
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff_loop = end-begin;
    myTiming.diff_loop_2 = myTiming.diff_loop_2 + myTiming.time_diff_loop.count();
    myTiming.bico_call = myTiming.bico_call + time_run_bico;
#endif
    
  // control over the number of iterations
  if( nrepeats < min_nrepeats ) isOK = false;

	// convergence control for environmental variables
	if( fabs(mG_envOld - mG_envNew) > eps*0.5*(fabs(mG_envOld)+fabs(mG_envNew))+TOL ) 
  {
		isOK = false;
		convergence_fail[0]++;
		double newprec = 2.*fabs(mG_envOld - mG_envNew)/(fabs(mG_envOld)+fabs(mG_envNew));
		prec = newprec > prec ? newprec : prec;
  }
//	if( fabs(mO2_envOld - mO2_envNew) > eps*0.5*(fabs(mO2_envOld)+fabs(mO2_envNew))+TOL ) 
//		{
//		isOK = false;
//		convergence_fail[1]++;
//		double newprec = 2.*fabs(mO2_envOld - mO2_envNew)/(fabs(mO2_envOld)+fabs(mO2_envNew));
//		prec = newprec > prec ? newprec : prec;
//		}
	if( fabs(mA_envOld - mA_envNew) > eps*0.5*(fabs(mA_envOld)+fabs(mA_envNew))+TOL ) 
  {
		isOK = false;
		convergence_fail[2]++;
		double newprec = 2.*fabs(mA_envOld - mA_envNew)/(fabs(mA_envOld)+fabs(mA_envNew));
		prec = newprec > prec ? newprec : prec;
	}

	if( fabs(mAcL_envOld - mAcL_envNew) > eps*0.5*(fabs(mAcL_envOld)+fabs(mAcL_envNew))+TOL )
  {
		isOK = false;
		convergence_fail[3]++;
		double newprec = 2*fabs(mAcL_envOld - mAcL_envNew)/(fabs(mAcL_envOld)+fabs(mAcL_envNew));
		prec = newprec > prec ? newprec : prec;
	}

	// assignment of environmental variables
	volume_envOld = volume_envNew;
	mG_envOld = mG_envNew;
	mO2_envOld = mO2_envNew;
	mA_envOld = mA_envNew;
	mAcL_envOld = mAcL_envNew;

	int ncell_fails = 0;										// variable that counts the number of cells for which convergence fails

#ifdef W_timing
    begin = std::chrono::steady_clock::now();
#endif
//#pragma omp parallel for ordered schedule(dynamic)
#pragma omp parallel for
	for(unsigned long n=0; n<ncells; n++)  // assignment of new values in preparation for the next iteration
  {
		bool cellisOK = true;
		
		// First check if the algorithm has reached convergence ...
		if( fabs(mGinOld[n] - mGinNew[n]) > eps*0.5*(fabs(mGinOld[n])+fabs(mGinNew[n]))+TOL ) 
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[4]++;
			double newprec = 2.*fabs(mGinOld[n] - mGinNew[n])/(fabs(mGinOld[n])+fabs(mGinNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mGinOld[n] - mGinNew[n]) > 0.5*(fabs(mGinOld[n])+fabs(mGinNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mGin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mGinNew[n]-mGinOld[n])/mGinOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mGextOld[n] - mGextNew[n]) > eps*0.5*(fabs(mGextOld[n])+fabs(mGextNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[5]++;
			double newprec = 2.*fabs(mGextOld[n] - mGextNew[n])/(fabs(mGextOld[n])+fabs(mGextNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mGextOld[n] - mGextNew[n]) > 0.5*(fabs(mGextOld[n])+fabs(mGextNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mGext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mGextNew[n]-mGextOld[n])/mGextOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mO2Old[n] - mO2New[n]) > eps*0.5*(fabs(mO2Old[n])+fabs(mO2New[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[6]++;
			double newprec = 2.*fabs(mO2Old[n] - mO2New[n])/(fabs(mO2Old[n])+fabs(mO2New[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mO2Old[n] - mO2New[n]) > 0.5*(fabs(mO2Old[n])+fabs(mO2New[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mO2: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mO2New[n]-mO2Old[n])/mO2Old[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mAinOld[n] - mAinNew[n]) > eps*0.5*(fabs(mAinOld[n])+fabs(mAinNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[7]++;
			double newprec = 2.*fabs(mAinOld[n] - mAinNew[n])/(fabs(mAinOld[n])+fabs(mAinNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mAinOld[n] - mAinNew[n]) > 0.5*(fabs(mAinOld[n])+fabs(mAinNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mAin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mAinNew[n]-mAinOld[n])/mAinOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mAextOld[n] - mAextNew[n]) > eps*0.5*(fabs(mAextOld[n])+fabs(mAextNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[8]++;
			double newprec = 2.*fabs(mAextOld[n] - mAextNew[n])/(fabs(mAextOld[n])+fabs(mAextNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mAextOld[n] - mAextNew[n]) > 0.5*(fabs(mAextOld[n])+fabs(mAextNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mAext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mAextNew[n]-mAextOld[n])/mAextOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mAcLinOld[n] - mAcLinNew[n]) > eps*0.5*(fabs(mAcLinOld[n])+fabs(mAcLinNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[9]++;
			double newprec = 2.*fabs(mAcLinOld[n] - mAcLinNew[n])/(fabs(mAcLinOld[n])+fabs(mAcLinNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mAcLinOld[n] - mAcLinNew[n]) > 0.5*fabs(fabs(mAcLinOld[n])+fabs(mAcLinNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mAcLin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mAcLinNew[n]-mAcLinOld[n])/mAcLinOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(mAcLextOld[n] - mAcLextNew[n]) > eps*0.5*(fabs(mAcLextOld[n])+fabs(mAcLextNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[10]++;
			double newprec = 2.*fabs(mAcLextOld[n] - mAcLextNew[n])/(fabs(mAcLextOld[n])+fabs(mAcLextNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(mAcLextOld[n] - mAcLextNew[n]) > 0.5*fabs(fabs(mAcLextOld[n])+fabs(mAcLextNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "mAcLext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(mAcLextNew[n]-mAcLextOld[n])/mAcLextOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		if( fabs(ATPpOld[n] - ATPpNew[n]) > eps*0.5*(fabs(ATPpOld[n])+fabs(ATPpNew[n]))+TOL )
    {
			isOK = false;
			cellisOK = false;
			convergence_fail[11]++;
			double newprec = 2.*fabs(ATPpOld[n] - ATPpNew[n])/(fabs(ATPpOld[n])+fabs(ATPpNew[n]));
			prec = newprec > prec ? newprec : prec;
			if ( fabs(ATPpOld[n] - ATPpNew[n]) > 0.5*(fabs(ATPpOld[n])+fabs(ATPpNew[n])) && EXTENDED_ERRORLOG ) 
      {
				errorlog_file << "CellsSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << std::endl;
				errorlog_file << "ATPp: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
				errorlog_file << "\tdifferenza : " << 100.*(ATPpNew[n]-ATPpOld[n])/ATPpOld[n] << "%\n";
				errorlog_file << "\tfase : " << phase[n] << std::endl;
				errorlog_file << "\teta' di fase: " << phase_age[n] << " s\n" << std::endl;
      }
    }
		
		if( !cellisOK ) ncell_fails++;
		
		// qui si rimpiazzano le vecchie approssimazioni con le nuove
		volumeOld[n] = volumeNew[n];
		volume_extraOld[n] = volume_extraNew[n];
		MitOld[n] = MitNew[n];
		pHiOld[n] = pHiNew[n];
		pHOld[n] = pHNew[n];
		mGinOld[n] = mGinNew[n];
		mGextOld[n] = mGextNew[n];
		mG6POld[n] = mG6PNew[n];
		mO2Old[n] = mO2New[n];
		StoreOld[n] = StoreNew[n];
		mAinOld[n] = mAinNew[n];
		mAextOld[n] = mAextNew[n];
		mAcLinOld[n] = mAcLinNew[n];
		mAcLextOld[n] = mAcLextNew[n];
		ATPpOld[n] = ATPpNew[n];
  }//end for(unsigned long n=0; n<ncells; n++)  // assignment of new values in preparation for the next iteration
#ifdef W_timing
    end= std::chrono::steady_clock::now();
    myTiming.time_diff_loop = end-begin;
    myTiming.diff_loop_3 = myTiming.diff_loop_3 + myTiming.time_diff_loop.count();
#endif
	nrepeats++;
	if( nrepeats >= MAXREPEATS )
  {
		errorlog_file << "\n*** ATTENZIONE: al passo " << nstep << " potenziale problema di convergenza dell'algoritmo in CellsSystem::Diff ***" << std::endl;
		errorlog_file << "precisione limitata a " << std::scientific << prec << " per " << ncell_fails << " cellule" << "\n" << std::endl;
	}
	
	} while( !isOK && nrepeats < MAXREPEATS );
	// end of the loop for the solution of the equations



	// calcolo del flusso di AcL nell'ambiente
	AcLFlow = ( mAcL_envOld - Env.GetEnvironmentAcL() )/dt;

	// assegnazione delle variabili ambientali
	Env.SetEnvironmentvolume(volume_envOld);
	Env.SetEnvironmentG(mG_envOld);
	Env.SetEnvironmentO2(mO2_envOld);
	Env.SetEnvironmentA(mA_envOld);
	Env.SetEnvironmentAcL(mAcL_envOld);
	Env.SetEnvironmentpH();
	

	// assegnazioni finali

	// assegnazioni comuni a cellule vive e a cellule morte: si tratta delle variabili relative alla diffusione che continua anche 
	// nel caso delle cellule morte
	pH = pHOld;
	G_extra = mGextOld;
	O2 = mO2Old;
	A_extra = mAextOld;
	AcL_extra = mAcLextOld;

#pragma omp parallel for
	for(unsigned long n=0; n<ncells; n++)
  {
		// assegnazioni che si fanno solo nel caso in cui la cellula e' viva
		if( phase[n] != dead )					
    {
			pHi[n] = pHiOld[n];
			G[n] = mGinOld[n];
			G6P[n] = mG6POld[n];
			store[n] = StoreOld[n];
			A[n] = mAinOld[n];
			AcL[n] = mAcLinOld[n];

			
			// ******* CAUTION !!! *******
			// this instruction defines and stores also radius, surface, volume, vol_extra, mass and ATPmin
			// WARNING: the DNA setting instruction must be set BEFORE the volume instruction, etc. because the value of the
			// Calculated DNA is used for volume setting, otherwise there is no consistency ...
			DNA[n] = DNANew[n];

			M[n] = MitOld[n];
			ATPp[n] = ATPpOld[n]; 
			
			ATPmin[n] = (type[n]->fATPmin)*(type[n]->C2 * M[n])/type[n]->C1;
			
			volume[n] = type[n]->Vmin * (1.+DNA[n]) + type[n]->C2 * M[n] + ATPp[n] * type[n]->C1;
			r[n] = pow(3.*volume[n]/(4.*PI), (double)1./3.); 
			surface[n] = 4.*PI*r[n]*r[n]; 	
			mass[n] = type[n]->density * volume[n];	
			
			volume_extra[n] = surface[n]*(type[n]->extvolume_thickness)*(type[n]->extvolume_fraction);

			// altre variabili e quantita' collegate
			protein[n] = proteinNew[n];
			prot_rate[n] = delta_protein[n]/dt;
						
			pRb[n] = pRbNew[n];
			
			if( phase[n] == G1m_phase)
      {
				cyclinD[n] += type[n]->Get_cyclinD_fraction()*delta_protein[n];
      }
			if( phase[n] == G1p_phase)
      {
				cyclinE[n] += type[n]->Get_cyclinE_fraction()*delta_protein[n];
      }
			if( phase[n] == G2_phase)
      {
				cyclinX[n] += type[n]->Get_cyclinX_fraction()*delta_protein[n];
      }
			
			ConcS[n] = ConcSNew[n];
						
			// variabili ausiliarie 
			
			rATPprod = ATP_Ox[n] + ATP_NOx[n] + ATP2[n] + ATP3[n];
			rATPcons = ConsATP[n] + ConsATP_1[n] + ConsATP_2[n] + ConsATP_3[n] + ConsATP_5[n];
			
			ATPprod[n] += rATPprod*dt;						// ATP prodotto in questo timestep
			ATPcons[n] += rATPcons*dt;						// ATP consumato in questo timestep
			ATPtot[n] = rATPprod - rATPcons;				// rate di variazione totale dell'ATPp

    }
		else	
		// the cell is dead and there are only residual processes 												
    {			

      // here the old concentrations are stored
			double concG6P = G6P[n]/volume[n];
      double concstore = store[n]/volume[n];
      double concATPp = ATPp[n]/volume[n];
      double concprotein = protein[n]/volume[n];
      double concpRb = pRb[n]/volume[n];
      double conccyclD = cyclinD[n]/volume[n];
      double conccyclE = cyclinE[n]/volume[n];
      double conccyclX = cyclinX[n]/volume[n];
      double concM = M[n]/volume[n];
            
            
      // *********
      // this instruction corresponds to the exact solution of the equation for the volume
      volume[n] = volumeOld[n]*exp(-dt*type[n]->Get_DVap()); 
      
      r[n] = pow(3.*volume[n]/(4.*PI), (double)1./3.); 
      surface[n] = 4.*PI*r[n]*r[n]; 
      mass[n] = type[n]->density * volume[n];
      
      volume_extra[n] = surface[n]*(type[n]->extvolume_thickness)*(type[n]->extvolume_fraction);
      // *********
            
            
      // here we redefine the masses of substances with constant concentration in dead cells
      G6P[n] = concG6P * volume[n];
      store[n] = concstore * volume[n];
      ATPp[n] = concATPp * volume[n];
 			protein[n] = concprotein * volume[n];
			pRb[n] = concpRb * volume[n];
      cyclinD[n] = conccyclD * volume[n];
      cyclinE[n] = conccyclE * volume[n];
      cyclinX[n] = conccyclX * volume[n];
      M[n] = concM * volume[n];
      
      // the DNA of dead cells is zeroed
      DNA[n] = 0;
    }

/* TF
 * see above
 */
		// controllo finale di consistenza
// 		if( phase[n] != dead ) 
//     {
// 			int code = CheckMVA(n);
// 			if(code < 0) errorlog_file << "Errore " << code << " alla fine di CellsSystem::Diff nel controllo di consistenza per la cellula " << n << "\n" << std::endl;
//     }
  }
}



