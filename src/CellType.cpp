// 
// metodi della classe CellType che contiene le informazioni sul fenotipo cellulare
//
// EM 19/1/2008
//
// **********************************************************************************

#include "sim.h"
#include "InputFromFile.h"

#include "CellType.h"
#include <cstdlib>
using namespace vbl;


// costruttore default
CellType::CellType()
{
	
	name = 0;					// nome del fenotipo standard 
	n_instances = 0;			// inizialmente non ci sono cellule con questo fenotipo
	
    VMAX_1 = 4e-6;  				// pg/(s·micron^2)
    VMAX_2 = 2.4e-4; 				// pg/s
    VMAX_22 = 2.4e-3;				// pg/s
    VMAX_A = 1.e-6;  				// pg/(s·micron^2) 
    VMAX_P = 3.e-4;  				// pg/s  
    VMAX_P_A = 1.018e-5;  			// pg/s
    VMAX_P_ATP = 1.81e-3;  			// pg/s 
    VMAX_DNA = 5.e-5;  				// molecole/s 
    VMAX_DNA_A = 5.847e-5;  		// pg/s 
    VMAX_DNA_ATP = 9.01e-5;  		// pg/s
    VMAX_M = 2e-3;  				// mitocondri/s 
    VMAX_M_A = 7.3e-8;  			// pg/s 
    VMAX_M_ATP = 1.125e-7;  		// pg/s 
    Km1 = 0.27024e-3;  				// pg/micron^3
    Km2 = 1.80e-3;  				// pg/micron^3
    Km22 = 1.80e-5;  				// pg/micron^3
    KmA = 0.023798e-3;  			// pg/micron^3
    Ka = 0.054e-3;  				// pg/micron^3
    Kmc = 0.096e-3;  				// pg/micron^3
    Kmd = 1.8e-5;  					// pg/micron^3
    KmO2 = 7.e-7;  					// pg/micron^3
    Kmp = 0.0067e-6;  				// (pg/micron^3)^2  
    KmDNA = 0.00045e-6;  			// (pg/micron^3)^2  
    KmM = 1.46e-8;  				// (pg/micron^3)^2  
    coeffg1 = 0.0075;  				// s^-1   
    coeffg2 = 0.00108;  			// s^-1   
    coeffg3 = 0.00067;  			// s^-1  
    coeffr1 = 3.e-5;  				// pg/s
    ATPSt = 2.3e-3;  				// pg/s ATP standard
    Vmin = 90;  					// micron^3
    DVap = 3e-6;  					// s^-1 DVap coefficiente di variazione di volume (shrinkage) per cellule apoptotiche
    fATPmin = 1.2;  				// frazione del volume cellulare e mitocondriale accessibile all'ATP
    pHimin = 6.8;  					// pH intracellulare minimo tollerabile
    VmaxAL0 = 9.583585e-5;  		// pg/(s·micron^2) 
    KmAL = 0.40536e-3;  			// pg/micron^3
    M_T_MEAN = 1800.;  				// s durata media della fase M
    DNA_MAX_SPREAD = 0.1;  			// fluttuazione massima di DNA_FRACTION nella singola cellula
    v_WORK = 1.5e-5;  				// pg/(s·micron^3) 
    PHASE_SPREAD = 0.5;  			// 
    k_pRb = 10;  					// 
    N_pRb = 16;  					// 
    pRb_ONOFFratio = 1.e-6;  		// k_ON/k_OFF in unita' di concentrazione molare per la fosforilazione della pRb
    pRb_fraction = 1.5e-2;  		// frazione delle proteine che e' pRb
    cyclinD_fraction = 10.e-3;  	// frazione delle proteine che e' ciclina D
    cyclinE_fraction = 3.e-3;  		// frazione delle proteine che e' ciclina E
    cyclinX_fraction = 14.e-3;  	// frazione delle proteine che e' ciclina A + B
    ConcS_0 = 0.001;  				// concentrazione MOLARE iniziale della sostanza S
    Thresh_S_start = 0.8;  			// frazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1m-G1p
    Thresh_S_stop = 0.05;  			// frazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1p-S
    k3MM = 1.e4;  					// rate per la reazione MM downstream
    KmMM = 1.e-3;  					// concentrazione MOLARE della Km per la MM downstream
    NUCLEAR_OBJ = 46;  				// numero di pseudooggetti legati alla matrice nucleare che tengono legata la pRb e che si separano a meta' circa al momento della mitosi (NB. la fluttuazione relativa e' uguale a 1/sqrt(NUCLEAR_OBJ) )
    ClusteringFactor = 15;  		// numero di compartimenti in cui si aggregano i mitocondri (costruzione fatta per adattare la stocasticita' a quella osservata)
    CycXThr = 0.08;  				// pg quantita' di ciclina A+B che definisce la soglia del checkpoint G2-M
    Vrif = -0.021;  				// V potenziale di Nernst in condizioni normali (interno a pH 7.2, esterno a pH 7.54)
    HPumpEff = 0.1;  				// inverso dell'efficienza del meccanismo di espulsione degli H+ (10 = efficienza 0.1)
    DiffH = 1.40e2;  				// micron^2/s Coefficiente di diffusione per gli H+
    C1 = 100;  						// micron^3/pg
    C2 = 0.2;  						// micron^3
    a_R = 5.e-2;  					// micron^3/(pg s) coefficiente che collega la concentrazione dell'AcL al danno endogeno
    alpha_R[G0_phase] = 0.;  		// 
    alpha_R[G1m_phase] = 0.43;  	// 
    alpha_R[G1p_phase] = 0.43;  	// 
    alpha_R[S_phase] = 0.17;  		// 
    alpha_R[G2_phase] = 0.77;  		// 
    alpha_R[M_phase] = 0.77;  		// 
    beta_R[G0_phase] = 0.;  		// 
    beta_R[G1m_phase] = 0.017;  	// 
    beta_R[G1p_phase] = 0.017;  	// 
    beta_R[S_phase] = 0.02;  		// 
    beta_R[G2_phase] = 0.;  		// 
    beta_R[M_phase] = 0.;  			// 
    YoungMod = 1e12;  				// pg/(micron·s^2) modulo di Young cellulare
    PoissonRatio = 0.5;  			// rapporto di Poisson cellulare
    density = 1.070;  				// pg/micron^3 densita' cellulare
    viscosity = 2e11;  				// pg/(micron·s) viscosita' citoplasmatica
    Mphase_correction = 1.1;  		// correzione al range della funzione di viscosita' durante la fase M e G1m
    adhesion_range = -0.5;  		// (frazione di raggio cellulare) range della forza di adesione
    adhesion_decay = 2.;  			// rate di decadimento della forza di adesione
    packing_factor = 0.9047;  		// fattore di impacchettamento
    extension_coeff = 1.1;  		// coefficiente di estensione
    extvolume_thickness = 0.1;  	// micron semispessore della matrice extracellulare
    extvolume_compression = 0.25;  	// micron-1 fattore di compressibilita' della matrice extracellulare
    extvolume_fraction = 0.3;  		// frazione di volume della matrice extracellulare utilizzabile per la diffusione
    tph_slope = 2.727273;  			// tuning dei rates in funzione del pH interno
    tph_thr = 6.55;  				// 
    tp11_slope = 10.90909;  		// tuning di p11 in funzione del pH interno
    tp11_thr = 6.9625;  			// 
    a2c_slope = 2.419355;  			// correzione al trasporto di glucosio da esterno a interno
    a2c_thr = 6.92;  				// 
    c2a_slope = 2.419355;  			// correzione al trasporto di glucosio da interno ad esterno
    c2a_thr = 6.92;  				// 
    a2cA_slope = 2.419355;  		// correzione al trasporto di glutammina da esterno a interno
    a2cA_thr = 6.92;  				// 
    c2aA_slope = 2.419355;  		// correzione al trasporto di glutammina da interno ad esterno
    c2aA_thr = 6.92;  				// 
    a2cAcL_slope = 1.5;  			// correzione al trasporto di AcL da esterno a interno
    a2cAcL_thr = 6.8;  				// 
    c2aAcL_slope = 1.5;  			// correzione al trasporto di AcL da interno a esterno
    c2aAcL_thr = 6.8;  				// 
}

// costruttore che prende i valori da un file

CellType::CellType(const std::string filename)
{

	std::ifstream ParameterFile( filename.c_str() );
	if( !ParameterFile ) 
		{
		std::cout << "CellType::CellType: errore nell'apertura del file dei parametri " << filename << std::endl;
		std::exit(-1);
		}

	name = 0;			// questo parametro viene assegnato dal programma
	n_instances = 0;	// questo parametro e' comunque nullo all'inizio
	
	VMAX_1 = InputRealPar(ParameterFile);
	VMAX_2 = InputRealPar(ParameterFile);
	VMAX_22 = InputRealPar(ParameterFile);
	VMAX_A = InputRealPar(ParameterFile);
	VMAX_P = InputRealPar(ParameterFile);
	VMAX_P_A = InputRealPar(ParameterFile);
	VMAX_P_ATP = InputRealPar(ParameterFile);
	VMAX_DNA = InputRealPar(ParameterFile);
	VMAX_DNA_A = InputRealPar(ParameterFile);
	VMAX_DNA_ATP = InputRealPar(ParameterFile);
	VMAX_M = InputRealPar(ParameterFile);
	VMAX_M_A = InputRealPar(ParameterFile);
	VMAX_M_ATP = InputRealPar(ParameterFile);
	Km1 = InputRealPar(ParameterFile);
	Km2 = InputRealPar(ParameterFile);
	Km22 = InputRealPar(ParameterFile);
	KmA = InputRealPar(ParameterFile);
	Ka = InputRealPar(ParameterFile);
	Kmc = InputRealPar(ParameterFile);
	Kmd = InputRealPar(ParameterFile);
	KmO2 = InputRealPar(ParameterFile);
	Kmp = InputRealPar(ParameterFile);
	KmDNA = InputRealPar(ParameterFile);
	KmM = InputRealPar(ParameterFile);
	coeffg1 = InputRealPar(ParameterFile);
	coeffg2 = InputRealPar(ParameterFile);
	coeffg3 = InputRealPar(ParameterFile);
	coeffr1 = InputRealPar(ParameterFile);
	ATPSt = InputRealPar(ParameterFile);
	Vmin = InputRealPar(ParameterFile);
	DVap = InputRealPar(ParameterFile);
	fATPmin = InputRealPar(ParameterFile);
	pHimin = InputRealPar(ParameterFile);
	VmaxAL0 = InputRealPar(ParameterFile);
	KmAL = InputRealPar(ParameterFile);
	M_T_MEAN = InputRealPar(ParameterFile);
	DNA_MAX_SPREAD = InputRealPar(ParameterFile);
	v_WORK = InputRealPar(ParameterFile);
	PHASE_SPREAD = InputRealPar(ParameterFile);
	k_pRb = InputIntPar(ParameterFile);
	N_pRb = InputIntPar(ParameterFile);
	pRb_ONOFFratio = InputRealPar(ParameterFile);
	pRb_fraction = InputRealPar(ParameterFile);
	cyclinD_fraction = InputRealPar(ParameterFile);
	cyclinE_fraction = InputRealPar(ParameterFile);
	cyclinX_fraction = InputRealPar(ParameterFile);
	ConcS_0 = InputRealPar(ParameterFile);
	Thresh_S_start = InputRealPar(ParameterFile);
	Thresh_S_stop = InputRealPar(ParameterFile);
	k3MM = InputRealPar(ParameterFile);
	KmMM = InputRealPar(ParameterFile);
	NUCLEAR_OBJ = InputIntPar(ParameterFile);
	ClusteringFactor = InputIntPar(ParameterFile);
	CycXThr = InputRealPar(ParameterFile);
	Vrif = InputRealPar(ParameterFile);
	HPumpEff = InputRealPar(ParameterFile);
	DiffH = InputRealPar(ParameterFile);
	C1 = InputRealPar(ParameterFile);
	C2 = InputRealPar(ParameterFile);
	a_R = InputRealPar(ParameterFile);
	alpha_R[G0_phase] = InputRealPar(ParameterFile);
	alpha_R[G1m_phase] = InputRealPar(ParameterFile);
	alpha_R[G1p_phase] = InputRealPar(ParameterFile);
	alpha_R[S_phase] = InputRealPar(ParameterFile);
	alpha_R[G2_phase] = InputRealPar(ParameterFile);
	alpha_R[M_phase] = InputRealPar(ParameterFile);
	alpha_R[dead] = 0.;	// valore default
	beta_R[G0_phase] = InputRealPar(ParameterFile);
	beta_R[G1m_phase] = InputRealPar(ParameterFile);
	beta_R[G1p_phase] = InputRealPar(ParameterFile);
	beta_R[S_phase] = InputRealPar(ParameterFile);
	beta_R[G2_phase] = InputRealPar(ParameterFile);
	beta_R[M_phase] = InputRealPar(ParameterFile);
	beta_R[dead] = 0.;	// valore default
	YoungMod = InputRealPar(ParameterFile);
	PoissonRatio = InputRealPar(ParameterFile);
	density = InputRealPar(ParameterFile);
	viscosity = InputRealPar(ParameterFile);
	Mphase_correction = InputRealPar(ParameterFile);
	adhesion_range = InputRealPar(ParameterFile);
	adhesion_decay = InputRealPar(ParameterFile);
	packing_factor = InputRealPar(ParameterFile);
	extension_coeff = InputRealPar(ParameterFile);
	extvolume_thickness = InputRealPar(ParameterFile);
	extvolume_compression = InputRealPar(ParameterFile);
	extvolume_fraction = InputRealPar(ParameterFile);

	tph_slope = InputRealPar(ParameterFile); 
	tph_thr = InputRealPar(ParameterFile);
	tp11_slope = InputRealPar(ParameterFile);
	tp11_thr = InputRealPar(ParameterFile);
	a2c_slope = InputRealPar(ParameterFile);
	a2c_thr = InputRealPar(ParameterFile);
	c2a_slope = InputRealPar(ParameterFile);
	c2a_thr = InputRealPar(ParameterFile);
	a2cA_slope = InputRealPar(ParameterFile);
	a2cA_thr = InputRealPar(ParameterFile);
	c2aA_slope = InputRealPar(ParameterFile);
	c2aA_thr = InputRealPar(ParameterFile);
	a2cAcL_slope = InputRealPar(ParameterFile);
	a2cAcL_thr = InputRealPar(ParameterFile);
	c2aAcL_slope = InputRealPar(ParameterFile);
	c2aAcL_thr = InputRealPar(ParameterFile);
	
}

// costruttore copia
CellType::CellType(const CellType& cCellType)
{

	int k;
	
	name = cCellType.name;
	n_instances = 0;	// ATTENZIONE: quando un nuovo fenotipo viene creato, il numero di instance e' sempre nullo ... 
	
	VMAX_1 = cCellType.VMAX_1;
	VMAX_2 = cCellType.VMAX_2;
	VMAX_22 = cCellType.VMAX_22;
	VMAX_A = cCellType.VMAX_A;
	VMAX_P = cCellType.VMAX_P;
	VMAX_P_A = cCellType.VMAX_P_A;
	VMAX_P_ATP = cCellType.VMAX_P_ATP;
	VMAX_DNA = cCellType.VMAX_DNA;
	VMAX_DNA_A = cCellType.VMAX_DNA_A;
	VMAX_DNA_ATP = cCellType.VMAX_DNA_ATP;
	VMAX_M = cCellType.VMAX_M;
	VMAX_M_A = cCellType.VMAX_M_A;
	VMAX_M_ATP = cCellType.VMAX_M_ATP;
	Km1 = cCellType.Km1;
	Km2 = cCellType.Km2;
	Km22 = cCellType.Km22;
	KmA = cCellType.KmA;
	Ka = cCellType.Ka;
	Kmc = cCellType.Kmc;
	Kmd = cCellType.Kmd;
	KmO2 = cCellType.KmO2;
	Kmp = cCellType.Kmp;
	KmDNA = cCellType.KmDNA;
	KmM = cCellType.KmM;
	coeffg1 = cCellType.coeffg1;
	coeffg2 = cCellType.coeffg2;
	coeffg3 = cCellType.coeffg3;
	coeffr1 = cCellType.coeffr1;
	ATPSt = cCellType.ATPSt;
	Vmin = cCellType.Vmin;
	DVap = cCellType.DVap;
	fATPmin = cCellType.fATPmin;
	pHimin = cCellType.pHimin;
	VmaxAL0 = cCellType.VmaxAL0;
	KmAL = cCellType.KmAL;
	M_T_MEAN = cCellType.M_T_MEAN;
	DNA_MAX_SPREAD = cCellType.DNA_MAX_SPREAD;
	v_WORK = cCellType.v_WORK;
	PHASE_SPREAD = cCellType.PHASE_SPREAD;
	k_pRb = cCellType.k_pRb;
	N_pRb = cCellType.N_pRb;
	pRb_ONOFFratio = cCellType.pRb_ONOFFratio;
	pRb_fraction = cCellType.pRb_fraction;
	cyclinD_fraction = cCellType.cyclinD_fraction;
	cyclinE_fraction = cCellType.cyclinE_fraction;
	cyclinX_fraction = cCellType.cyclinX_fraction;
	ConcS_0 = cCellType.ConcS_0;
	Thresh_S_start = cCellType.Thresh_S_start;
	Thresh_S_stop = cCellType.Thresh_S_stop;
	k3MM = cCellType.k3MM;
	KmMM = cCellType.KmMM;
	NUCLEAR_OBJ = cCellType.NUCLEAR_OBJ;
	ClusteringFactor = cCellType.ClusteringFactor;
	CycXThr = cCellType.CycXThr;
	Vrif = cCellType.Vrif;
	HPumpEff = cCellType.HPumpEff;
	DiffH = cCellType.DiffH;
	C1 = cCellType.C1;
	C2 = cCellType.C2;
	a_R = cCellType.a_R;
	for(k=0; k< Nphase; k++)
		{
		alpha_R[k] = cCellType.alpha_R[k];
		beta_R[k] = cCellType.beta_R[k];
		}
	YoungMod = cCellType.YoungMod;
	PoissonRatio = cCellType.PoissonRatio;
	density = cCellType.density;
	viscosity = cCellType.viscosity;
	Mphase_correction = cCellType.Mphase_correction;
	adhesion_range = cCellType.adhesion_range;
	adhesion_decay = cCellType.adhesion_decay;
	packing_factor = cCellType.packing_factor;
	extension_coeff = cCellType.extension_coeff;
	extvolume_thickness = cCellType.extvolume_thickness;
	extvolume_compression = cCellType.extvolume_compression;
	extvolume_fraction = cCellType.extvolume_fraction;

	tph_slope = cCellType.tph_slope; 
	tph_thr = cCellType.tph_thr;
	tp11_slope = cCellType.tp11_slope;
	tp11_thr = cCellType.tp11_thr;
	a2c_slope = cCellType.a2c_slope;
	a2c_thr = cCellType.a2c_thr;
	c2a_slope = cCellType.c2a_slope;
	c2a_thr = cCellType.c2a_thr;
	a2cA_slope = cCellType.a2cA_slope;
	a2cA_thr = cCellType.a2cA_thr;
	c2aA_slope = cCellType.c2aA_slope;
	c2aA_thr = cCellType.c2aA_thr;
	a2cAcL_slope = cCellType.a2cAcL_slope;
	a2cAcL_thr = cCellType.a2cAcL_thr;
	c2aAcL_slope = cCellType.c2aAcL_slope;
	c2aAcL_thr = cCellType.c2aAcL_thr;


}

CellType& CellType::operator=(const CellType& cCellType)
{
	int k;

	name = cCellType.name;
	n_instances = 0;	// ATTENZIONE: quando un nuovo fenotipo viene creato, il numero di instance e' sempre nullo ... 

	VMAX_1 = cCellType.VMAX_1;
	VMAX_2 = cCellType.VMAX_2;
	VMAX_22 = cCellType.VMAX_22;
	VMAX_A = cCellType.VMAX_A;
	VMAX_P = cCellType.VMAX_P;
	VMAX_P_A = cCellType.VMAX_P_A;
	VMAX_P_ATP = cCellType.VMAX_P_ATP;
	VMAX_DNA = cCellType.VMAX_DNA;
	VMAX_DNA_A = cCellType.VMAX_DNA_A;
	VMAX_DNA_ATP = cCellType.VMAX_DNA_ATP;
	VMAX_M = cCellType.VMAX_M;
	VMAX_M_A = cCellType.VMAX_M_A;
	VMAX_M_ATP = cCellType.VMAX_M_ATP;
	Km1 = cCellType.Km1;
	Km2 = cCellType.Km2;
	Km22 = cCellType.Km22;
	KmA = cCellType.KmA;
	Ka = cCellType.Ka;
	Kmc = cCellType.Kmc;
	Kmd = cCellType.Kmd;
	KmO2 = cCellType.KmO2;
	Kmp = cCellType.Kmp;
	KmDNA = cCellType.KmDNA;
	KmM = cCellType.KmM;
	coeffg1 = cCellType.coeffg1;
	coeffg2 = cCellType.coeffg2;
	coeffg3 = cCellType.coeffg3;
	coeffr1 = cCellType.coeffr1;
	ATPSt = cCellType.ATPSt;
	Vmin = cCellType.Vmin;
	DVap = cCellType.DVap;
	fATPmin = cCellType.fATPmin;
	pHimin = cCellType.pHimin;
	VmaxAL0 = cCellType.VmaxAL0;
	KmAL = cCellType.KmAL;
	M_T_MEAN = cCellType.M_T_MEAN;
	DNA_MAX_SPREAD = cCellType.DNA_MAX_SPREAD;
	v_WORK = cCellType.v_WORK;
	PHASE_SPREAD = cCellType.PHASE_SPREAD;
	k_pRb = cCellType.k_pRb;
	N_pRb = cCellType.N_pRb;
	pRb_ONOFFratio = cCellType.pRb_ONOFFratio;
	pRb_fraction = cCellType.pRb_fraction;
	cyclinD_fraction = cCellType.cyclinD_fraction;
	cyclinE_fraction = cCellType.cyclinE_fraction;
	cyclinX_fraction = cCellType.cyclinX_fraction;
	ConcS_0 = cCellType.ConcS_0;
	Thresh_S_start = cCellType.Thresh_S_start;
	Thresh_S_stop = cCellType.Thresh_S_stop;
	k3MM = cCellType.k3MM;
	KmMM = cCellType.KmMM;
	NUCLEAR_OBJ = cCellType.NUCLEAR_OBJ;
	ClusteringFactor = cCellType.ClusteringFactor;
	CycXThr = cCellType.CycXThr;
	Vrif = cCellType.Vrif;
	HPumpEff = cCellType.HPumpEff;
	DiffH = cCellType.DiffH;
	C1 = cCellType.C1;
	C2 = cCellType.C2;
	a_R = cCellType.a_R;
	for(k=0; k< Nphase; k++)
		{
		alpha_R[k] = cCellType.alpha_R[k];
		beta_R[k] = cCellType.beta_R[k];
		}
	YoungMod = cCellType.YoungMod;
	PoissonRatio = cCellType.PoissonRatio;
	density = cCellType.density;
	viscosity = cCellType.viscosity;
	Mphase_correction = cCellType.Mphase_correction;
	adhesion_range = cCellType.adhesion_range;
	adhesion_decay = cCellType.adhesion_decay;
	packing_factor = cCellType.packing_factor;
	extension_coeff = cCellType.extension_coeff;
	extvolume_thickness = cCellType.extvolume_thickness;
	extvolume_compression = cCellType.extvolume_compression;
	extvolume_fraction = cCellType.extvolume_fraction;

	tph_slope = cCellType.tph_slope; 
	tph_thr = cCellType.tph_thr;
	tp11_slope = cCellType.tp11_slope;
	tp11_thr = cCellType.tp11_thr;
	a2c_slope = cCellType.a2c_slope;
	a2c_thr = cCellType.a2c_thr;
	c2a_slope = cCellType.c2a_slope;
	c2a_thr = cCellType.c2a_thr;
	a2cA_slope = cCellType.a2cA_slope;
	a2cA_thr = cCellType.a2cA_thr;
	c2aA_slope = cCellType.c2aA_slope;
	c2aA_thr = cCellType.c2aA_thr;
	a2cAcL_slope = cCellType.a2cAcL_slope;
	a2cAcL_thr = cCellType.a2cAcL_thr;
	c2aAcL_slope = cCellType.c2aAcL_slope;
	c2aAcL_thr = cCellType.c2aAcL_thr;
	
	return *this;

}

// funzione che stampa su cout i valori 
void PrintCellType(CellType &cCellType)
{

	std::cout << "Size of CellType: " << sizeof(cCellType) << std::endl;
	
	std::cout << cCellType << std::endl;
	
}

// // funzione di stampa semplice (overloaded <<)
namespace vbl{ // see https://stackoverflow.com/questions/10712380/undefined-reference-to-operator
std::ostream& operator<<(std::ostream& stream, CellType const & cCellType)
{
	int k;
	
	stream << " cell type: " << cCellType.name << "\n" << std::endl;
	stream << " instances: " << cCellType.n_instances << "\n" << std::endl;	

	stream << " VMAX_1: " << cCellType.VMAX_1 << " \tpg*s^-1*micron^-2" << std::endl;
	stream << " VMAX_2: " << cCellType.VMAX_2 << " \tpg*s^-1" << std::endl;
	stream << " VMAX_22: " << cCellType.VMAX_22 << " \tpg*s^-1" << std::endl;
	stream << " VMAX_A: " << cCellType.VMAX_A << " \tpg*s^-1*micron^-2" << std::endl;
	stream << " VMAX_P: " << cCellType.VMAX_P << " \tpg*s^-1" << std::endl;
	stream << " VMAX_P_A: " << cCellType.VMAX_P_A << " \tpg*s^-1" << std::endl;
	stream << " VMAX_P_ATP: " << cCellType.VMAX_P_ATP << " \tpg*s^-1" << std::endl;
	stream << " VMAX_DNA: " << cCellType.VMAX_DNA << " \tmol*s^-1" << std::endl;
	stream << " VMAX_DNA_A: " << cCellType.VMAX_DNA_A << " \tpg*s^-1" << std::endl;
	stream << " VMAX_DNA_ATP: " << cCellType.VMAX_DNA_ATP << " \tpg*s^-1" << std::endl;
	stream << " VMAX_M: " << cCellType.VMAX_M << " \tmitocondri*s^-1" << std::endl;
	stream << " VMAX_M_A: " << cCellType.VMAX_M_A << " \tpg*s^-1" << std::endl;
	stream << " VMAX_M_ATP: " << cCellType.VMAX_M_ATP << " \tpg*s^-1" << std::endl;
	stream << " Km1: " << cCellType.Km1 << " \tpg*micron^-3" << std::endl;
	stream << " Km2: " << cCellType.Km2 << " \tpg*micron^-3" << std::endl;
	stream << " Km22: " << cCellType.Km22 << " \tpg*micron^-3" << std::endl;
	stream << " KmA: " << cCellType.KmA << " \tpg*micron^-3" << std::endl;
	stream << " Ka: " << cCellType.Ka << " \tpg*micron^-3" << std::endl;
	stream << " Kmc: " << cCellType.Kmc << " \tpg*micron^-3" << std::endl;
	stream << " Kmd: " << cCellType.Kmd << " \tpg*micron^-3" << std::endl;
	stream << " KmO2: " << cCellType.KmO2 << " \tpg*micron^-3" << std::endl;
	stream << " Kmp: " << cCellType.Kmp << " \t(pg*micron^-3)^2" << std::endl;
	stream << " KmDNA: " << cCellType.KmDNA << " \t(pg*micron^-3)^2" << std::endl;
	stream << " KmM: " << cCellType.KmM << " \t(pg*micron^-3)^2" << std::endl;
	stream << " coeffg1: " << cCellType.coeffg1 << " \ts^-1" << std::endl;
	stream << " coeffg2: " << cCellType.coeffg2 << " \ts^-1" << std::endl;
	stream << " coeffg3: " << cCellType.coeffg3 << " \ts^-1" << std::endl;
	stream << " coeffr1: " << cCellType.coeffr1 << " \ts^-1" << std::endl;
	stream << " ATPSt: " << cCellType.ATPSt << " \tpg*s-1 ATP standard" << std::endl;
	stream << " Vmin: " << cCellType.Vmin << " \tmicron^3" << std::endl;
	stream << " DVap: " << cCellType.DVap << " \ts^-1 DVap coefficiente di variazione di volume (shrinkage) per cellule apoptotiche" << std::endl;
	stream << " fATPmin: " << cCellType.fATPmin << " \t frazione del volume cellulare e mitocondriale accessibile all'ATP" << std::endl;
	stream << " pHimin: " << cCellType.pHimin << " \tpH intracellulare minimo tollerabile" << std::endl;
	stream << " VmaxAL0: " << cCellType.VmaxAL0 << " \tpg*s^-1*micron^-2" << std::endl;
	stream << " KmAL: " << cCellType.KmAL << " \tpg*micron^-3" << std::endl;
	stream << " M_T_MEAN: " << cCellType.M_T_MEAN << " \ts durata media della fase M" << std::endl;
	stream << " DNA_MAX_SPREAD: " << cCellType.DNA_MAX_SPREAD << " \tfluttuazione massima di DNA_FRACTION nella singola cellula" << std::endl;
	stream << " v_WORK: " << cCellType.v_WORK << " \tpg/s/micron^3 consumo di ATP per sostenere il funzionamento generale della cellula" << std::endl;
	stream << " PHASE_SPREAD: " << cCellType.PHASE_SPREAD << " \ts fluttuazione della durata di G2 e M dovuta ad effetti meccanici" << std::endl;
	stream << " k_pRb: " << cCellType.k_pRb << " \tk soglia pRB" << std::endl;
	stream << " N_pRb: " << cCellType.N_pRb << " \tN siti pRb" << std::endl;
	stream << " pRb_ONOFFratio: " << cCellType.pRb_ONOFFratio << " \tk_ON/k_OFF in unita' di concentrazione molare per la fosforilazione della pRb" << std::endl;
	stream << " pRb_fraction: " << cCellType.pRb_fraction << " \tfrazione delle proteine che e' pRb" << std::endl;
	stream << " cyclinD_fraction: " << cCellType.cyclinD_fraction << " \tfrazione delle proteine che e' ciclina D" << std::endl;
	stream << " cyclinE_fraction: " << cCellType.cyclinE_fraction << " \tfrazione delle proteine che e' ciclina E" << std::endl;
	stream << " cyclinX_fraction: " << cCellType.cyclinX_fraction << " \tfrazione delle proteine che e' ciclina A + B" << std::endl;
	stream << " ConcS_0: " << cCellType.ConcS_0 << " \tconcentrazione MOLARE iniziale della sostanza S" << std::endl;
	stream << " Thresh_S_start: " << cCellType.Thresh_S_start << " \tfrazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1m-G1p" << std::endl;
	stream << " Thresh_S_stop: " << cCellType.Thresh_S_stop << " \tfrazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1p-S" << std::endl;
	stream << " k3MM: " << cCellType.k3MM << " \ts^-1 rate per la reazione MM downstream" << std::endl;
	stream << " KmMM: " << cCellType.KmMM << " \tconcentrazione MOLARE della Km per la MM downstream" << std::endl;
	stream << " NUCLEAR_OBJ: " << cCellType.NUCLEAR_OBJ << " \tnumero di pseudooggetti legati alla matrice nucleare che tengono legata la pRb e che si separano a meta' circa al momento della mitosi (NB. la fluttuazione relativa e' uguale a 1/sqrt(NUCLEAR_OBJ) )" << std::endl;
	stream << " ClusteringFactor: " << cCellType.ClusteringFactor << " \tnumero di compartimenti in cui si aggregano i mitocondri (costruzione fatta per adattare la stocasticita' a quella osservata)" << std::endl;
	stream << " CycXThr: " << cCellType.CycXThr << " \tpg quantita' di ciclina A+B che definisce la soglia del checkpoint G2-M" << std::endl;
	stream << " Vrif: " << cCellType.Vrif << " \tV potenziale di Nernst in condizioni normali (interno a pH 7.2, esterno a pH 7.54)" << std::endl;
	stream << " HPumpEff: " << cCellType.HPumpEff << " \tinverso dell'efficienza del meccanismo di espulsione degli H+ (10 = efficienza 0.1)" << std::endl;
	stream << " DiffH: " << cCellType.DiffH << " \tmicron^2/s Coefficiente di diffusione per gli H+" << std::endl;
	stream << " C1: " << cCellType.C1 << "  \tparametri che collegano i mitocondri al volume" << std::endl;
	stream << " C2: " << cCellType.C2 << "  \tparametri che collegano i mitocondri al volume" << std::endl;
	stream << " a_R: " << cCellType.a_R << " \t(pg/micron^3)^-1 coefficiente che collega la concentrazione dell'AcL al danno endogeno" << std::endl;
	stream << " alpha_R: ";
		for(k=0; k<= Nphase-2; k++)
			stream << cCellType.alpha_R[k] << ", ";
		stream << cCellType.alpha_R[Nphase-1] << " \tGy^-1" << std::endl;
	stream << " beta_R: ";
		for(k=0; k<= Nphase-2; k++)
			stream << cCellType.beta_R[k] << ", ";
		stream << cCellType.beta_R[Nphase-1] << " \tGy^-2" << std::endl;
	stream << " YoungModulus: " << cCellType.YoungMod << " \tpg/(micron·s^2) modulo di Young cellulare" << std::endl;		
	stream << " PoissonRatio: " << cCellType.PoissonRatio << " \trapporto di Poisson cellulare" << std::endl;		
	stream << " density: " << cCellType.density << " \tpg/micron^3 densita' cellulare" << std::endl;
	stream << " viscosity: " << cCellType.viscosity << " \tpg/(micron·s) viscosita' citoplasmatica" << std::endl;
	stream << " Mphase_correction: " << cCellType.Mphase_correction << " \tcorrezione al range della funzione di viscosita' durante la fase M e G1m" << std::endl;	
	stream << " adhesion_range: " << cCellType.adhesion_range << " \tmicron range della forza di adesione" << std::endl;
	stream << " adhesion_decay: " << cCellType.adhesion_decay << " \t rate di decadimento della forza di adesione" << std::endl;
	stream << " packing_factor: " << cCellType.packing_factor << " \tfattore di impacchettamento" << std::endl;
	stream << " extension_coeff: " << cCellType.extension_coeff << " \tcoefficiente di estensione" << std::endl;
	stream << " extvolume_thickness: " << cCellType.extvolume_thickness << " \tmicron semispessore della matrice extracellulare" << std::endl;
	stream << " extvolume_compression: " << cCellType.extvolume_compression << " \tmicron^-1 fattore di compressibilita' della matrice extracellulare" << std::endl;
	stream << " extvolume_fraction: " << cCellType.extvolume_fraction << " \tfrazione di volume della matrice extracellulare utilizzabile per la diffusione" << std::endl;

	stream << " tph_slope: " << cCellType.tph_slope << " \t tuning dei rates in funzione del pH interno" << std::endl;
	stream << " tph_thr: " << cCellType.tph_thr << " \t" << std::endl;
	stream << " tp11_slope: " << cCellType.tp11_slope << " \t tuning di p11 in funzione del pH interno" << std::endl;
	stream << " tp11_thr: " << cCellType.tp11_thr << " \t" << std::endl;
	stream << " a2c_slope: " << cCellType.a2c_slope << " \t correzione al trasporto di glucosio da esterno a interno" << std::endl;
	stream << " a2c_thr: " << cCellType.a2c_thr << " \t" << std::endl;
	stream << " c2a_slope: " << cCellType.c2a_slope << " \t correzione al trasporto di glucosio da interno ad esterno" << std::endl;
	stream << " c2a_thr: " << cCellType.c2a_thr << " \t" << std::endl;
	stream << " a2cA_slope: " << cCellType.a2cA_slope << " \t correzione al trasporto di glutammina da esterno a interno" << std::endl;
	stream << " a2cA_thr: " << cCellType.a2cA_thr << " \t" << std::endl;
	stream << " c2aA_slope: " << cCellType.c2aA_slope << " \t correzione al trasporto di glutammina da interno ad esterno" << std::endl;
	stream << " c2aA_thr: " << cCellType.c2aA_thr << " \t" << std::endl;
	stream << " a2cAcL_slope: " << cCellType.a2cAcL_slope << " \t correzione al trasporto di AcL da esterno a interno" << std::endl;
	stream << " a2cAcL_thr: " << cCellType.a2cAcL_thr << " \t" << std::endl;
	stream << " c2aAcL_slope: " << cCellType.c2aAcL_slope << " \t correzione al trasporto di AcL da interno a esterno" << std::endl;
	stream << " c2aAcL_thr: " << cCellType.c2aAcL_thr << " \t" << std::endl;

	return stream;

}
}//namespace vbl{

// write binario su file del cell type
void CellType::WriteCellType( std::ofstream& stream )
{

	stream.write( (char*)(&name), sizeof(int) );
	stream.write( (char*)(&n_instances), sizeof(unsigned long) );	

	stream.write( (char*)(&VMAX_1), sizeof(double) );
	stream.write( (char*)(&VMAX_2), sizeof(double) );
	stream.write( (char*)(&VMAX_22), sizeof(double) );
	stream.write( (char*)(&VMAX_A), sizeof(double) );
	stream.write( (char*)(&VMAX_P), sizeof(double) );
	stream.write( (char*)(&VMAX_P_A), sizeof(double) );
	stream.write( (char*)(&VMAX_P_ATP), sizeof(double) );
	stream.write( (char*)(&VMAX_DNA), sizeof(double) );
	stream.write( (char*)(&VMAX_DNA_A), sizeof(double) );
	stream.write( (char*)(&VMAX_DNA_ATP), sizeof(double) );
	stream.write( (char*)(&VMAX_M), sizeof(double) );
	stream.write( (char*)(&VMAX_M_A), sizeof(double) );
	stream.write( (char*)(&VMAX_M_ATP), sizeof(double) );
	stream.write( (char*)(&Km1), sizeof(double) );
	stream.write( (char*)(&Km2), sizeof(double) );
	stream.write( (char*)(&Km22), sizeof(double) );
	stream.write( (char*)(&KmA), sizeof(double) );
	stream.write( (char*)(&Ka), sizeof(double) );
	stream.write( (char*)(&Kmc), sizeof(double) );
	stream.write( (char*)(&Kmd), sizeof(double) );
	stream.write( (char*)(&KmO2), sizeof(double) );
	stream.write( (char*)(&Kmp), sizeof(double) );
	stream.write( (char*)(&KmDNA), sizeof(double) );
	stream.write( (char*)(&KmM), sizeof(double) );
	stream.write( (char*)(&coeffg1), sizeof(double) );
	stream.write( (char*)(&coeffg2), sizeof(double) );
	stream.write( (char*)(&coeffg3), sizeof(double) );
	stream.write( (char*)(&coeffr1), sizeof(double) );
	stream.write( (char*)(&ATPSt), sizeof(double) );
	stream.write( (char*)(&Vmin), sizeof(double) );
	stream.write( (char*)(&DVap), sizeof(double) );
	stream.write( (char*)(&fATPmin), sizeof(double) );
	stream.write( (char*)(&pHimin), sizeof(double) );
	stream.write( (char*)(&VmaxAL0), sizeof(double) );
	stream.write( (char*)(&KmAL), sizeof(double) );
	stream.write( (char*)(&M_T_MEAN), sizeof(double) );
	stream.write( (char*)(&DNA_MAX_SPREAD), sizeof(double) );
	stream.write( (char*)(&v_WORK), sizeof(double) );
	stream.write( (char*)(&PHASE_SPREAD), sizeof(double) );
	stream.write( (char*)(&k_pRb), sizeof(int) );
	stream.write( (char*)(&N_pRb), sizeof(int) );
	stream.write( (char*)(&pRb_ONOFFratio), sizeof(double) );
	stream.write( (char*)(&pRb_fraction), sizeof(double) );
	stream.write( (char*)(&cyclinD_fraction), sizeof(double) );
	stream.write( (char*)(&cyclinE_fraction), sizeof(double) );
	stream.write( (char*)(&cyclinX_fraction), sizeof(double) );
	stream.write( (char*)(&ConcS_0), sizeof(double) );
	stream.write( (char*)(&Thresh_S_start), sizeof(double) );
	stream.write( (char*)(&Thresh_S_stop), sizeof(double) );
	stream.write( (char*)(&k3MM), sizeof(double) );
	stream.write( (char*)(&KmMM), sizeof(double) );
	stream.write( (char*)(&NUCLEAR_OBJ), sizeof(int) );
	stream.write( (char*)(&ClusteringFactor), sizeof(int) );
	stream.write( (char*)(&CycXThr), sizeof(double) );
	stream.write( (char*)(&Vrif), sizeof(double) );
	stream.write( (char*)(&HPumpEff), sizeof(double) );
	stream.write( (char*)(&DiffH), sizeof(double) );
	stream.write( (char*)(&C1), sizeof(double) );
	stream.write( (char*)(&C2), sizeof(double) );
	stream.write( (char*)(&a_R), sizeof(double) );
	for(int k=0; k< Nphase-1; k++)
		stream.write( (char*)(&alpha_R[k]), sizeof(double) );
	for(int k=0; k< Nphase-1; k++)
		stream.write( (char*)(&beta_R[k]), sizeof(double) );
	stream.write( (char*)(&YoungMod), sizeof(double) );	
	stream.write( (char*)(&PoissonRatio), sizeof(double) );		
	stream.write( (char*)(&density), sizeof(double) );
	stream.write( (char*)(&viscosity), sizeof(double) );
	stream.write( (char*)(&Mphase_correction), sizeof(double) );	
	stream.write( (char*)(&adhesion_range), sizeof(double) );
	stream.write( (char*)(&adhesion_decay), sizeof(double) );
	stream.write( (char*)(&packing_factor), sizeof(double) );
	stream.write( (char*)(&extension_coeff), sizeof(double) );
	stream.write( (char*)(&extvolume_thickness), sizeof(double) );
	stream.write( (char*)(&extvolume_compression), sizeof(double) );
	stream.write( (char*)(&extvolume_fraction), sizeof(double) );

	stream.write( (char*)(&tph_slope), sizeof(double) );
	stream.write( (char*)(&tph_thr), sizeof(double) );
	stream.write( (char*)(&tp11_slope), sizeof(double) );
	stream.write( (char*)(&tp11_thr), sizeof(double) );
	stream.write( (char*)(&a2c_slope), sizeof(double) );
	stream.write( (char*)(&a2c_thr), sizeof(double) );
	stream.write( (char*)(&c2a_slope), sizeof(double) );
	stream.write( (char*)(&c2a_thr), sizeof(double) );
	stream.write( (char*)(&a2cA_slope), sizeof(double) );
	stream.write( (char*)(&a2cA_thr), sizeof(double) );
	stream.write( (char*)(&c2aA_slope), sizeof(double) );
	stream.write( (char*)(&c2aA_thr), sizeof(double) );
	stream.write( (char*)(&a2cAcL_slope), sizeof(double) );
	stream.write( (char*)(&a2cAcL_thr), sizeof(double) );
	stream.write( (char*)(&c2aAcL_slope), sizeof(double) );
	stream.write( (char*)(&c2aAcL_thr), sizeof(double) );

	
}

// read binario da file del cell type
void CellType::ReadCellType( std::ifstream& stream )
{
	stream.read( (char*)(&name), sizeof(int) );
	stream.read( (char*)(&n_instances), sizeof(unsigned long) );	

	stream.read( (char*)(&VMAX_1), sizeof(double) );
	stream.read( (char*)(&VMAX_2), sizeof(double) );
	stream.read( (char*)(&VMAX_22), sizeof(double) );
	stream.read( (char*)(&VMAX_A), sizeof(double) );
	stream.read( (char*)(&VMAX_P), sizeof(double) );
	stream.read( (char*)(&VMAX_P_A), sizeof(double) );
	stream.read( (char*)(&VMAX_P_ATP), sizeof(double) );
	stream.read( (char*)(&VMAX_DNA), sizeof(double) );
	stream.read( (char*)(&VMAX_DNA_A), sizeof(double) );
	stream.read( (char*)(&VMAX_DNA_ATP), sizeof(double) );
	stream.read( (char*)(&VMAX_M), sizeof(double) );
	stream.read( (char*)(&VMAX_M_A), sizeof(double) );
	stream.read( (char*)(&VMAX_M_ATP), sizeof(double) );
	stream.read( (char*)(&Km1), sizeof(double) );
	stream.read( (char*)(&Km2), sizeof(double) );
	stream.read( (char*)(&Km22), sizeof(double) );
	stream.read( (char*)(&KmA), sizeof(double) );
	stream.read( (char*)(&Ka), sizeof(double) );
	stream.read( (char*)(&Kmc), sizeof(double) );
	stream.read( (char*)(&Kmd), sizeof(double) );
	stream.read( (char*)(&KmO2), sizeof(double) );
	stream.read( (char*)(&Kmp), sizeof(double) );
	stream.read( (char*)(&KmDNA), sizeof(double) );
	stream.read( (char*)(&KmM), sizeof(double) );
	stream.read( (char*)(&coeffg1), sizeof(double) );
	stream.read( (char*)(&coeffg2), sizeof(double) );
	stream.read( (char*)(&coeffg3), sizeof(double) );
	stream.read( (char*)(&coeffr1), sizeof(double) );
	stream.read( (char*)(&ATPSt), sizeof(double) );
	stream.read( (char*)(&Vmin), sizeof(double) );
	stream.read( (char*)(&DVap), sizeof(double) );
	stream.read( (char*)(&fATPmin), sizeof(double) );
	stream.read( (char*)(&pHimin), sizeof(double) );
	stream.read( (char*)(&VmaxAL0), sizeof(double) );
	stream.read( (char*)(&KmAL), sizeof(double) );
	stream.read( (char*)(&M_T_MEAN), sizeof(double) );
	stream.read( (char*)(&DNA_MAX_SPREAD), sizeof(double) );
	stream.read( (char*)(&v_WORK), sizeof(double) );
	stream.read( (char*)(&PHASE_SPREAD), sizeof(double) );
	stream.read( (char*)(&k_pRb), sizeof(int) );
	stream.read( (char*)(&N_pRb), sizeof(int) );
	stream.read( (char*)(&pRb_ONOFFratio), sizeof(double) );
	stream.read( (char*)(&pRb_fraction), sizeof(double) );
	stream.read( (char*)(&cyclinD_fraction), sizeof(double) );
	stream.read( (char*)(&cyclinE_fraction), sizeof(double) );
	stream.read( (char*)(&cyclinX_fraction), sizeof(double) );
	stream.read( (char*)(&ConcS_0), sizeof(double) );
	stream.read( (char*)(&Thresh_S_start), sizeof(double) );
	stream.read( (char*)(&Thresh_S_stop), sizeof(double) );
	stream.read( (char*)(&k3MM), sizeof(double) );
	stream.read( (char*)(&KmMM), sizeof(double) );
	stream.read( (char*)(&NUCLEAR_OBJ), sizeof(int) );
	stream.read( (char*)(&ClusteringFactor), sizeof(int) );
	stream.read( (char*)(&CycXThr), sizeof(double) );
	stream.read( (char*)(&Vrif), sizeof(double) );
	stream.read( (char*)(&HPumpEff), sizeof(double) );
	stream.read( (char*)(&DiffH), sizeof(double) );
	stream.read( (char*)(&C1), sizeof(double) );
	stream.read( (char*)(&C2), sizeof(double) );
	stream.read( (char*)(&a_R), sizeof(double) );
	for(int k=0; k< Nphase-1; k++)
		stream.read( (char*)(&alpha_R[k]), sizeof(double) );
	for(int k=0; k< Nphase-1; k++)
		stream.read( (char*)(&beta_R[k]), sizeof(double) );
	stream.read( (char*)(&YoungMod), sizeof(double) );	
	stream.read( (char*)(&PoissonRatio), sizeof(double) );		
	stream.read( (char*)(&density), sizeof(double) );
	stream.read( (char*)(&viscosity), sizeof(double) );
	stream.read( (char*)(&Mphase_correction), sizeof(double) );	
	stream.read( (char*)(&adhesion_range), sizeof(double) );
	stream.read( (char*)(&adhesion_decay), sizeof(double) );
	stream.read( (char*)(&packing_factor), sizeof(double) );
	stream.read( (char*)(&extension_coeff), sizeof(double) );
	stream.read( (char*)(&extvolume_thickness), sizeof(double) );
	stream.read( (char*)(&extvolume_compression), sizeof(double) );
	stream.read( (char*)(&extvolume_fraction), sizeof(double) );

	stream.read( (char*)(&tph_slope), sizeof(double) );
	stream.read( (char*)(&tph_thr), sizeof(double) );
	stream.read( (char*)(&tp11_slope), sizeof(double) );
	stream.read( (char*)(&tp11_thr), sizeof(double) );
	stream.read( (char*)(&a2c_slope), sizeof(double) );
	stream.read( (char*)(&a2c_thr), sizeof(double) );
	stream.read( (char*)(&c2a_slope), sizeof(double) );
	stream.read( (char*)(&c2a_thr), sizeof(double) );
	stream.read( (char*)(&a2cA_slope), sizeof(double) );
	stream.read( (char*)(&a2cA_thr), sizeof(double) );
	stream.read( (char*)(&c2aA_slope), sizeof(double) );
	stream.read( (char*)(&c2aA_thr), sizeof(double) );
	stream.read( (char*)(&a2cAcL_slope), sizeof(double) );
	stream.read( (char*)(&a2cAcL_thr), sizeof(double) );
	stream.read( (char*)(&c2aAcL_slope), sizeof(double) );
	stream.read( (char*)(&c2aAcL_thr), sizeof(double) );

	
}
