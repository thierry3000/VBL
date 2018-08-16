// 
// definizione della classe CellType che contiene le informazioni sul fenotipo cellulare
//
// EM 19/1/2008
//
// **********************************************************************************
#ifndef CELLTYPE_H
#define CELLTYPE_H  // header guard
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>
#include <vector>


namespace vbl{
// class CellsSystem; //forward declaration
class CellType 
{
// friends

		// funzione di stampa di tutti i dati della classe
	friend void PrintCellType( CellType &cCellType);		
	friend std::ostream& operator<<(std::ostream& s, CellType const & cCellType);
  void PrintCellType( CellType &cCellType);		
// 	std::ostream& operator<<(std::ostream& s, CellType& cCellType);
		
	private: 
//     std::vector<std::string> vector_of_cell_type_parameters = { 
//     "name","n_instances","VMAX_1","VMAX_2","VMAX_22","VMAX_A","VMAX_P","VMAX_P_A","VMAX_P_ATP","VMAX_DNA","VMAX_DNA_A","VMAX_DNA_ATP","VMAX_M","VMAX_M_A","VMAX_M_ATP","Km1","Km2","Km22","KmA","Ka","Kmc","Kmd","KmO2","Kmp","KmDNA","KmM","coeffg1","coeffg2","coeffg3","coeffr1","ATPSt","Vmin","DVap","fATPmin","pHimin","VmaxAL0","KmAL","M_T_MEAN","DNA_MAX_SPREAD","v_WORK","PHASE_SPREAD","k_pRb","N_pRb","pRb_ONOFFratio","pRb_fraction","cyclinD_fraction","cyclinE_fraction","cyclinX_fraction","ConcS_0","Thresh_S_start","Thresh_S_stop","k3MM","KmMM","NUCLEAR_OBJ","ClusteringFactor","CycXThr","Vrif","HPumpEff","DiffH","C1","C2","a_R","YoungMod","PoissonRatio","density","viscosity","Mphase_correction","adhesion_range","adhesion_decay","packing_factor","extension_coeff","extvolume_thickness","extvolume_compression","extvolume_fraction","tph_slope","tph_thr","tp11_slope","tp11_thr","a2c_slope","a2c_thr","c2a_slope","c2a_thr","a2cA_slope","a2cA_thr","c2aA_slope","c2aA_thr","a2cAcL_slope","a2cAcL_thr","c2aAcL_slope","c2aAcL_thr", "alpha_R", "beta_R"};
		//double alpha_R[Nphase];		// coefficiente alpha della legge lineare-quadratica
		//double beta_R[Nphase];			// coefficiente beta della legge lineare-quadratica
		
	public:
		// costruttori
		CellType();
		CellType(const std::string filename);
		// copia
		CellType(const CellType& cCellType);
    void assign( const boost::property_tree::ptree &pt);
    boost::property_tree::ptree as_ptree() const;
		// nessun distruttore, utilizzo il default del compilatore
		
		CellType& operator=(const CellType& ct);
				
		int name;							// nome (codice) del fenotipo (standard = 0)
		unsigned long n_instances;			// numero di cellule con questo fenotipo
		
		double VMAX_1;					// pg/(s·micron^2) VMAX_1
		double VMAX_2;					// pg/s VMAX_2
		double VMAX_22;				// pg/s VMAX_22
		double VMAX_A;					// pg/(s·micron^2) VMAX_A
		double VMAX_P;					// pg/s VMAX_P
		double VMAX_P_A;				// pg/s VMAX_P_A
		double VMAX_P_ATP;				// pg/s VMAX_P_ATP
		double VMAX_DNA;				// pg/s VMAX_DNA
		double VMAX_DNA_A;				// pg/s VMAX_DNA_A
		double VMAX_DNA_ATP;			// pg/s VMAX_DNA_ATP
		double VMAX_M;					// mitocondri/s VMAX_M
		double VMAX_M_A;				// pg/s VMAX_M_A
		double VMAX_M_ATP;				// pg/s VMAX_M_ATP
		double Km1;					// pg/micron^3 Km1
		double Km2;					// pg/micron^3 Km2
		double Km22;					// pg/micron^3 Km22
		double KmA;					// pg/micron^3 KmA
		double Ka;						// pg/micron^3 Ka
		double Kmc;					// pg/micron^3 Kmc
		double Kmd;					// pg/micron^3 Kmd
		double KmO2;					// pg/micron^3 KmO2
		double Kmp;					// (pg/micron^3)^2 Kmp
		double KmDNA;					// (pg/micron^3)^2 KmDNA
		double KmM;					// (pg/micron^3)^2 KmM
		double coeffg1;				// s^-1 coeff. g1
		double coeffg2;				// s^-1 coeff. g2
		double coeffg3;				// s^-1 coeff. g3
		double coeffr1;				// pg/s coeff. r1
		double ATPSt;					// pg/s ATP standard
		double Vmin;					// micron^3 Vmin  volume residuo minimo delle cellule apoptotiche
		double DVap;					// s^-1 DVap	coefficiente di variazione di volume (shrinkage) per cellule apoptotiche
		double fATPmin;				// moltiplicatore dell'ATPmin calcolato sulla base dell'occupazione del volume dei mitocondri
		double pHimin;					// pH intracellulare minimo tollerabile
		double VmaxAL0;				// pg/(s·micron^2) VmaxAL0
		double KmAL;					// pg/micron^3 KmAL
		double M_T_MEAN;				// durata media fase M
		double DNA_MAX_SPREAD;			// fluttuazione massima di DNA_FRACTION nella singola cellula
		double v_WORK;					// pg/(s·micron^3) consumo di ATP per sostenere il funzionamento generale della cellula
		double PHASE_SPREAD;			// fluttuazione della durata di G2 e M dovuta ad effetti "meccanici"
		int k_pRb;							// k threshold pRB
		int N_pRb;							// N sites pRb
		double pRb_ONOFFratio;			// k_ON/k_OFF in unità di concentrazione molare per la fosforilazione della pRb
		double pRb_fraction;			// frazione delle proteine che e' pRb
		double cyclinD_fraction;		// frazione delle proteine che e' ciclina D
		double cyclinE_fraction;		// frazione delle proteine che e' ciclina E
		double cyclinX_fraction;		// frazione delle proteine che e' ciclina A + B
		 
// i parametri di soglia seguenti servono a decidere l'attraversamento dei checkpoint G1m-G1p e G1p-S, a partire dal contenuto di 
// S che viene gradualmente consumato dalla reazione MM downstream nel meccanismo della pRb
// S puo' essere identificato genericamente con cdc6, e la soglia start definisce l'inizio del consumo di cdc6, e stop la fine, 
// quando puo' finalmente iniziare la duplicazione del DNA (che si e' liberato dalla guaina di cdc6)
// visto che si parla di consumo deve essere sempre Thresh_S_start > Thresh_S_stop

		double ConcS_0;				// concentrazione MOLARE iniziale della sostanza S
		double Thresh_S_start;			// frazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1m-G1p
		double Thresh_S_stop;			// frazione di molecola S (nella reazione MM downstream) che fissa il passaggio del checkpoint G1p-S
		double k3MM;					// s^-1 rate per la reazione MM downstream
		double KmMM;					// concentrazione MOLARE della Km per la MM downstream
		int NUCLEAR_OBJ;					// numero di pseudooggetti legati alla matrice nucleare che tengono legata la pRb e che si separano a meta' circa al 
											// momento della mitosi (NB. la fluttuazione relativa e' uguale a 1/sqrt(NUCLEAR_OBJ) )
		int ClusteringFactor;				// fattore di aggregazione dei mitocondri (serve per calcolare la stocasticita' al momento della divisione cellulare)
		double CycXThr;				// quantita' di ciclina A+B che definisce la soglia del checkpoint G2-M
		double Vrif;					// V potenziale di Nernst in condizioni normali (interno a pH 7.2, esterno a pH 7.54)
		double HPumpEff;				// Inverso dell'efficienza del meccanismo di espulsione degli H+ (10 = efficienza 0.1)
		double DiffH;					// micron^2/s Coefficiente di diffusione per gli H+
		double C1;						// micron^3/pg 
		double C2;                     // micron^3
		double a_R;					// micron^3/(pg s) coefficiente che collega la concentrazione dell'AcL al danno endogeno
		double alpha_R[Nphase];		// coefficiente alpha della legge lineare-quadratica
		double beta_R[Nphase];			// coefficiente beta della legge lineare-quadratica
		
		double YoungMod;				// pg/(micron·s^2) modulo di Young della cellula
		double PoissonRatio;			// rapporto di Poisson
		double density;				// pg/micron^3 densita' citoplasmatica di questo tipo cellulare
		double viscosity;				// pg/(micron·s) viscosita' citoplasmatica effettiva di questo tipo cellulare
		double Mphase_correction;		// correzione al range della funzione di viscosita' durante la fase M e G1m
		double adhesion_range;			// range dell'adesione cellulare in unità di raggio cellulare (in realta' questo e' un indicatore dell'estensibilita' cellulare)
		double adhesion_decay;			// decay rate della forza di adesione
		double packing_factor;			// fattore di impacchettamento: il fattore ottimale di impacchettamento di sfere rigide identiche e' 
											// 0.7405 in volume, e questo significa che con un raggio ridotto di un fattore 
											// 0.7405^(1/3) = 0.9047 si ottiene lo stesso volume totale nell'ipotesi di impacchettamento perfetto
											// il packing factor e' proprio il fattore di riduzione del raggio (0.9047 per sfere rigide identiche)
		double extension_coeff;		// il coeff. di estensione da' una misura di quanto una cellula si estende per restare a contatto 
											// con le vicine; in pratica questo corrisponde a calcolare la sup. di contatto tra cellule sferiche con
											// raggio equivalente diverso da quello della cellula. Si potrebbe stimare questo valore minimizzando la
											// somma di energia potenziale elastica e di energia di contatto dovuta all'adesione, pero' mentre 
											// scrivo non so come stimare l'energia di contatto di adesione, e quindi nel file che definisce il 
											// fenotipo metto un valore "ragionevole" uguale a 1.1
		double extvolume_thickness;	// micron, semi-spessore della regione extracellulare
		double extvolume_compression;	// micron^-1, fattore di compressibilita' della regione extracellulare
		double extvolume_fraction;		// frazione del volume extracellulare utilizzabile per la diffusione

// parametri che regolano le funzioni a soglia per l'acidita' 
// slope e' la pendenza della tanh, thr la posizione della soglia, le funzioni sono tutte del tipo 
// 0.5*(1 + tanh( slope*(pH-thr) ) )	
		double tph_slope;				// tuning dei rates in funzione del pH interno 
		double tph_thr;
		double tp11_slope;				// tuning di p11 in funzione del pH interno
		double tp11_thr;
		double a2c_slope;				// correzione al trasporto di glucosio da esterno a interno
		double a2c_thr;
		double c2a_slope;				// correzione al trasporto di glucosio da interno ad esterno
		double c2a_thr;
		double a2cA_slope;				// correzione al trasporto di glutammina da esterno a interno
		double a2cA_thr;
		double c2aA_slope;				// correzione al trasporto di glutammina da interno ad esterno
		double c2aA_thr;
		double a2cAcL_slope;			// correzione al trasporto di AcL da esterno a interno
		double a2cAcL_thr;
		double c2aAcL_slope;			// correzione al trasporto di AcL da interno a esterno
		double c2aAcL_thr;

		
// getters

	int Get_name() { return name; };
	unsigned long Get_n_instances() { return n_instances; };
	
	double Get_VMAX_1() { return VMAX_1; };
	double Get_VMAX_2() { return VMAX_2; };
	double Get_VMAX_22() { return VMAX_22; };
	double Get_VMAX_A() { return VMAX_A; };
	double Get_VMAX_P() { return VMAX_P; };
	double Get_VMAX_P_A() { return VMAX_P_A; };
	double Get_VMAX_P_ATP() { return VMAX_P_ATP; };
	double Get_VMAX_DNA() { return VMAX_DNA; };
	double Get_VMAX_DNA_A() { return VMAX_DNA_A; };
	double Get_VMAX_DNA_ATP() { return VMAX_DNA_ATP; };
	double Get_VMAX_M() { return VMAX_M; };
	double Get_VMAX_M_A() { return VMAX_M_A; };
	double Get_VMAX_M_ATP() { return VMAX_M_ATP; };
	double Get_Km1() { return Km1; };
	double Get_Km2() { return Km2; };
	double Get_Km22() { return Km22; };
	double Get_KmA() { return KmA; };
	double Get_Ka() { return Ka; };
	double Get_Kmc() { return Kmc; };
	double Get_Kmd() { return Kmd; };
	double Get_KmO2() { return KmO2; };
	double Get_Kmp() { return Kmp; };
	double Get_KmDNA() { return KmDNA; };
	double Get_KmM() { return KmM; };
	double Get_coeffg1() { return coeffg1; };
	double Get_coeffg2() { return coeffg2; };
	double Get_coeffg3() { return coeffg3; };
	double Get_coeffr1() { return coeffr1; };
	double Get_ATPSt() { return ATPSt; };
	double Get_Vmin() { return Vmin; };
	double Get_DVap() { return DVap; };
	double Get_fATPmin() { return fATPmin; };
	double Get_pHimin() { return pHimin; };
	double Get_VmaxAL0() { return VmaxAL0; };
	double Get_KmAL() { return KmAL; };
	double Get_M_T_MEAN() { return M_T_MEAN; };
	double Get_DNA_MAX_SPREAD() { return DNA_MAX_SPREAD; };
	double Get_v_WORK() { return v_WORK; };
	double Get_PHASE_SPREAD() { return PHASE_SPREAD; };
	int Get_k_pRb() { return k_pRb; };
	int Get_N_pRb() { return N_pRb; };
	double Get_pRb_ONOFFratio() { return pRb_ONOFFratio; };
	double Get_pRb_fraction() { return pRb_fraction; };
	double Get_cyclinD_fraction() { return cyclinD_fraction; };
	double Get_cyclinE_fraction() { return cyclinE_fraction; };
	double Get_cyclinX_fraction() { return cyclinX_fraction; };
	double Get_ConcS_0() { return ConcS_0; };
	double Get_Thresh_S_start() { return Thresh_S_start; };
	double Get_Thresh_S_stop() { return Thresh_S_stop; };
	double Get_k3MM() { return k3MM; };
	double Get_KmMM() { return KmMM; };
	int Get_NUCLEAR_OBJ() { return NUCLEAR_OBJ; };
	int Get_ClusteringFactor() { return ClusteringFactor; };
	double Get_CycXThr() { return CycXThr; };
	double Get_Vrif() { return Vrif; };
	double Get_HPumpEff() { return HPumpEff; };
	double Get_DiffH() { return DiffH; };
	double Get_C1() { return C1; };
	double Get_C2() { return C2; };
	double Get_a_R() { return a_R; };
	double Get_alpha_R( CellPhase phaseindex ) { return alpha_R[phaseindex]; };
	double Get_beta_R( CellPhase phaseindex ) { return beta_R[phaseindex]; };
	double Get_YoungMod() { return YoungMod; };
	double Get_PoissonRatio() { return PoissonRatio; };
	double Get_density() { return density; };
	double Get_viscosity() { return viscosity; };
	double Get_Mphase_correction()  { return Mphase_correction; };
	double Get_adhesion_range() { return adhesion_range; };
	double Get_adhesion_decay() { return adhesion_decay; };
	double Get_packing_factor() { return packing_factor; };
	double Get_extension_coeff() { return extension_coeff; };
	double Get_extvolume_thickness() { return extvolume_thickness; };
	double Get_extvolume_compression() { return extvolume_compression; };
	double Get_extvolume_fraction() {return extvolume_fraction; };

	double Get_tph_slope() {return tph_slope; };
	double Get_tph_thr() {return tph_thr; };
	double Get_tp11_slope() {return tp11_slope; };
	double Get_tp11_thr() {return tp11_thr; };
	double Get_a2c_slope() {return a2c_slope; };
	double Get_a2c_thr() {return a2c_thr; };
	double Get_c2a_slope() {return c2a_slope; };
	double Get_c2a_thr() {return c2a_thr; };
	double Get_a2cA_slope() {return a2cA_slope; };
	double Get_a2cA_thr() {return a2cA_thr; };
	double Get_c2aA_slope() {return c2aA_slope; };
	double Get_c2aA_thr() {return c2aA_thr; };
	double Get_a2cAcL_slope() {return a2cAcL_slope; };
	double Get_a2cAcL_thr() {return a2cAcL_thr; };
	double Get_c2aAcL_slope() {return c2aAcL_slope; };
	double Get_c2aAcL_thr() {return c2aAcL_thr; };


// setters

	void Set_name( const int cname ) { name = cname; };
	void Set_n_instances( const int cn_instances ) { n_instances = cn_instances; };
	void New_instance() { n_instances++; };
	void Delete_instance() { if(n_instances > 0) n_instances--; };
	
	void Set_VMAX_1(const double cVMAX_1) { VMAX_1 = cVMAX_1; };
	void Set_VMAX_2(const double cVMAX_2) { VMAX_2 = cVMAX_2; };
	void Set_VMAX_22(const double cVMAX_22) { VMAX_22 = cVMAX_22; };
	void Set_VMAX_A(const double cVMAX_A) { VMAX_A = cVMAX_A; };
	void Set_VMAX_P(const double cVMAX_P) { VMAX_P = cVMAX_P; };
	void Set_VMAX_P_A(const double cVMAX_P_A) { VMAX_P_A = cVMAX_P_A; };
	void Set_VMAX_P_ATP(const double cVMAX_P_ATP) { VMAX_P_ATP = cVMAX_P_ATP; };
	void Set_VMAX_DNA(const double cVMAX_DNA) { VMAX_DNA = cVMAX_DNA; };
	void Set_VMAX_DNA_A(const double cVMAX_DNA_A) { VMAX_DNA_A = cVMAX_DNA_A; };
	void Set_VMAX_DNA_ATP(const double cVMAX_DNA_ATP) { VMAX_DNA_ATP = cVMAX_DNA_ATP; };
	void Set_VMAX_M(const double cVMAX_M) { VMAX_M = cVMAX_M; };
	void Set_VMAX_M_A(const double cVMAX_M_A) { VMAX_M_A = cVMAX_M_A; };
	void Set_VMAX_M_ATP(const double cVMAX_M_ATP) { VMAX_M_ATP = cVMAX_M_ATP; };
	void Set_Km1(const double cKm1) { Km1 = cKm1; };
	void Set_Km2(const double cKm2) { Km2 = cKm2; };
	void Set_Km22(const double cKm22) { Km22 = cKm22; };
	void Set_KmA(const double cKmA) { KmA = cKmA; };
	void Set_Ka(const double cKa) { Ka = cKa; };
	void Set_Kmc(const double cKmc) { Kmc = cKmc; };
	void Set_Kmd(const double cKmd) { Kmd = cKmd; };
	void Set_KmO2(const double cKmO2) { KmO2 = cKmO2; };
	void Set_Kmp(const double cKmp) { Kmp = cKmp; };
	void Set_KmDNA(const double cKmDNA) { KmDNA = cKmDNA; };
	void Set_KmM(const double cKmM) { KmM = cKmM; };
	void Set_coeffg1(const double ccoeffg1) { coeffg1 = ccoeffg1; };
	void Set_coeffg2(const double ccoeffg2) { coeffg2 = ccoeffg2; };
	void Set_coeffg3(const double ccoeffg3) { coeffg3 = ccoeffg3; };
	void Set_coeffr1(const double ccoeffr1) { coeffr1 = ccoeffr1; };
	void Set_ATPSt(const double cATPSt) { ATPSt = cATPSt; };
	void Set_Vmin(const double cVmin) { Vmin = cVmin; };
	void Set_DVap(const double cDVap) { DVap = cDVap; };
	void Set_fATPmin(const double cfATPmin) { fATPmin = cfATPmin; };
	void Set_pHimin(const double cpHimin) { pHimin = cpHimin; };
	void Set_VmaxAL0(const double cVmaxAL0) { VmaxAL0 = cVmaxAL0; };
	void Set_KmAL(const double cKmAL) { KmAL = cKmAL; };
	void Set_M_T_MEAN(const double cM_T_MEAN) { M_T_MEAN = cM_T_MEAN; };
	void Set_DNA_MAX_SPREAD(const double cDNA_MAX_SPREAD) { DNA_MAX_SPREAD = cDNA_MAX_SPREAD; };
	void Set_v_WORK(const double cv_WORK) { v_WORK = cv_WORK; };
	void Set_PHASE_SPREAD(const double cPHASE_SPREAD) { PHASE_SPREAD = cPHASE_SPREAD; };
	void Set_k_pRb(const int ck_pRb) { k_pRb = ck_pRb; };
	void Set_N_pRb(const int cN_pRb) { N_pRb = cN_pRb; };
	void Set_pRb_ONOFFratio(const double cpRb_ONOFFratio) { pRb_ONOFFratio = cpRb_ONOFFratio; };
	void Set_pRb_fraction(const double cpRb_fraction) { pRb_fraction = cpRb_fraction; };
	void Set_cyclinD_fraction(const double ccyclinD_fraction) { cyclinD_fraction = ccyclinD_fraction; };
	void Set_cyclinE_fraction(const double ccyclinE_fraction) { cyclinE_fraction = ccyclinE_fraction; };
	void Set_cyclinX_fraction(const double ccyclinX_fraction) { cyclinX_fraction = ccyclinX_fraction; };
	void Set_ConcS_0(const double cConcS_0) { ConcS_0 = cConcS_0; };
	void Set_Thresh_S_start(const double cThresh_S_start) { Thresh_S_start = cThresh_S_start; };
	void Set_Thresh_S_stop(const double cThresh_S_stop) { Thresh_S_stop = cThresh_S_stop; };
	void Set_k3MM(const double ck3MM) { k3MM = ck3MM; };
	void Set_KmMM(const double cKmMM) { KmMM = cKmMM; };
	void Set_NUCLEAR_OBJ(const int cNUCLEAR_OBJ) { NUCLEAR_OBJ = cNUCLEAR_OBJ; };
	void Set_ClusteringFactor(const int cClusteringFactor) { ClusteringFactor = cClusteringFactor; };
	void Set_CycXThr(const double cCycXThr) { CycXThr = cCycXThr; };
	void Set_Vrif(const double cVrif) { Vrif = cVrif; };
	void Set_HPumpEff(const double cHPumpEff) { HPumpEff = cHPumpEff; };
	void Set_DiffH(const double cDiffH) { DiffH = cDiffH; };
	void Set_C1(const double cC1) { C1 = cC1; };
	void Set_C2(const double cC2) { C2 = cC2; };
	void Set_a_R(const double ca_R) { a_R = ca_R; };
	void Set_alpha_R( const CellPhase phaseindex, const double calpha_R ) { alpha_R[phaseindex] = calpha_R; };
	void Set_beta_R( const CellPhase phaseindex, const double cbeta_R ) { beta_R[phaseindex] = cbeta_R; };
	void Set_YoungMod( const double newYoungMod ) { YoungMod = newYoungMod; };
	void Set_PoissonRatio( const double newPoissonRatio) { PoissonRatio = newPoissonRatio; };
	void Set_density( const double newdensity ) { density = newdensity; }; 
	void Set_viscosity( const double newviscosity ) { viscosity = newviscosity; };
	void Set_Mphase_correction( const double newMphase_correction ) { Mphase_correction = newMphase_correction; };
	void Set_adhesion_range( const double newadhesion_range ) { adhesion_range = newadhesion_range; };
	void Set_adhesion_decay( const double newadhesion_decay ) { adhesion_decay = newadhesion_decay; };
	void Set_packing_factor( const double newpacking_factor ) { packing_factor = newpacking_factor; };
	void Set_extension_coeff( const double newextension_coeff ) { extension_coeff = newextension_coeff; };
	void Set_extvolume_thickness( const double newextvolume_thickness ) { extvolume_thickness = newextvolume_thickness; };
	void Set_extvolume_compression( const double newextvolume_compression ) { extvolume_compression = newextvolume_compression; };
	void Set_extvolume_fraction( const double newextvolume_fraction ) { extvolume_fraction = newextvolume_fraction; };

	void Set_tph_slope( const double newtph_slope ) { tph_slope = newtph_slope; };
	void Set_tph_thr( const double newtph_thr ) { tph_thr = newtph_thr; };
	void Set_tp11_slope( const double newtp11_slope ) { tp11_slope = newtp11_slope; };
	void Set_tp11_thr( const double newtp11_thr ) { tp11_thr = newtp11_thr; };
	void Set_a2c_slope( const double newa2c_slope ) { a2c_slope = newa2c_slope; };
	void Set_a2c_thr( const double newa2c_thr ) { a2c_thr = newa2c_thr; };
	void Set_c2a_slope( const double newc2a_slope ) { c2a_slope = newc2a_slope; };
	void Set_c2a_thr( const double newc2a_thr ) { c2a_thr = newc2a_thr; };
	void Set_a2cA_slope( const double newa2cA_slope ) { a2cA_slope = newa2cA_slope; };
	void Set_a2cA_thr( const double newa2cA_thr ) { a2cA_thr = newa2cA_thr; };
	void Set_c2aA_slope( const double newc2aA_slope ) { c2aA_slope = newc2aA_slope; };
	void Set_c2aA_thr( const double newc2aA_thr ) { c2aA_thr = newc2aA_thr; };
	void Set_a2cAcL_slope( const double newa2cAcL_slope ) { a2cAcL_slope = newa2cAcL_slope; };
	void Set_a2cAcL_thr( const double newa2cAcL_thr ) { a2cAcL_thr = newa2cAcL_thr; };
	void Set_c2aAcL_slope( const double newc2aAcL_slope ) { c2aAcL_slope = newc2aAcL_slope; };
	void Set_c2aAcL_thr( const double newc2aAcL_thr ) { c2aAcL_thr = newc2aAcL_thr; };
	
// read and write binario
	void WriteCellType( std::ofstream& stream );
	void ReadCellType( std::ifstream& stream );

};
//std::ostream& operator<<(std::ostream& s, CellType const & cCellType);
}//namespace vbl{
//std::ostream& vbl::operator<<(std::ostream& s, CellType const & cCellType);
#endif //#ifndef CELLTYPE_H
