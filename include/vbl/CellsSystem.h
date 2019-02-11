/*
 *  CellsSystem.h
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 20/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */

// Definition of the CellsSystem class that contains the cell structure
// and methods for operations on all cells
#ifndef CELLSSYSTEM_H
#define CELLSSYSTEM_H // header guard

#include <CGAL/Default.h>

#include <boost/optional.hpp>
#include <sys/utsname.h>	// header per i metodi di identificazione della macchina

#include "sim.h"
#include "InputFromFile.h"
#include "CellType.h"
#include "Environment.h"
#include "EnvironmentalSignals.h"
#include "BloodVessel.h"
#include "Utilities.h"
#include "MutEventCreator.h"
#include "geometry.h"

#ifdef _parallel
  #include <tbb/tbb.h>
  #include <tbb/concurrent_vector.h>
#endif
#include <mutex>

#include <limits>
#include <iterator>

#ifdef _OPENMP
  #include <omp.h>
#endif

#if VBL_USE_TUMORCODE
  //#include <mwlib/helpers-vec.h>
  //#include "common/hdfio.h"
  //#include "mwlib/lattice-data.h"
  #include "common/cell_based_oxygen_update_model.h"
#endif

#define W_timing
//#define serial

namespace vbl{
// #ifndef NDEBUG
//   // means debugging
//   #define W_timing
// #endif
//class CellsSystem;//forward declaration


struct ReadInParameters
{
  /* imulation with dispersed cells (0) or full 3D (1) */
  // questo serve a identificare il tipo di simulazione
  int sim_type;
  int run;						// numero del run
  double dt;					// timestep
  double dt_sm;               // timestep slow motion
  double t;					// tempo attuale di simulazione
  double t_ini;				// tempo di inizializzazione
  double treal;				// tempo di simulazione vero (dopo l'inizializzazione)
  double tmax;				// tempo max di simulazione
  double tsm_start;           // tempo di start dello slow motion
  double tsm_stop;            // tempo di stop dello slow motion
  bool slow_motion;           // indica lo stato dello slow motion
  double t_CPU_max;			// tempo di CPU massimo per una singola frazione di run (in s)
  unsigned long nstep;			// numero del passo
  unsigned long nstep_start;		// numero del passo alla partenza della simulazione (quando ready2start diventa true)
  unsigned long nmax;				// nmax = floor(tmax/dt) numero max di passi
  int idum;						// variabile globale utilizzata dalle routine di numeri casuali di Numerical Recipes 2
  unsigned long nprint;			// intervallo (in numero di passi) tra i passi di stampa su file
  unsigned long nscreen;			// intervallo (in numero di passi) tra i passi di stampa su schermo
  unsigned long nconfiguration;	// numero delle configurazioni scritte su file
  // precisione
  double eps;				// precisione del passo di diffusione
  double delta_vmax;			// precisione della determinazione di velocita'
  // numero di cellule 
  unsigned long ncells;
  
  // numero di cellule vive alla fine del passo di metabolismo
  unsigned long alive;

  // numero di tipi cellulari
  unsigned long ntypes;
  // flusso (flag che indica se il flusso e' non nullo, e classe che rappresenta il segnale non nullo)
  bool flowON;
  // dose rate (flag che indica se la dose e' non nulla, e classe che rappresenta il segnale non nullo)
  bool doseON;
  
  
  void assign(const boost::property_tree::ptree &pt);
  boost::property_tree::ptree as_ptree() const;
};

/** @brief 
 * Structure to log some timings.
 * This helps to analyze the runtimes
 */
template <class T>
struct Timing
{
  T diff = 0;
  T dynamics =0;
  T geometry = 0;
  T geometry_neighborhood = 0;
  T cellEvents = 0;
  T writeToFile = 0;
  T diff_loop_1 = 0;
  T diff_loop_2 = 0;
  T diff_loop_3 = 0;
  T bico_call = 0;
  std::chrono::duration<T> time_diff;
  std::chrono::duration<T> time_diff_loop;
  std::chrono::duration<T> time_run_bico;
  std::chrono::duration<T> time_geometry_neighborhood;
  void reset()
  {
    diff=0;
    diff_loop_1 = 0;
    diff_loop_2 = 0;
    diff_loop_3 = 0;
    bico_call = 0;
    dynamics=0;
    geometry=0;
    cellEvents=0;
    geometry_neighborhood = 0;
  };
};

class CellsSystem
{
private:
//prameters
ReadInParameters params;
// output file
std::string output_filename;
std::ofstream output_file;

// log file
std::string log_filename;
std::ofstream log_file;

// screen dump file
std::string screen_dump_filename;
std::ofstream screen_dump_file;

// configuration file
std::string configuration_filename;
std::ofstream configuration_file;

// configuration file (binary)
std::string configuration_b_filename;
std::ofstream configuration_b_file;

// error log file
std::string errorlog_filename;
std::ofstream errorlog_file;

// convergence log file
std::string convlog_filename;
std::ofstream convlog_file;

// cell log file
std::string cell_filename;
std::ofstream cell_file;

// environment log file
std::string env_filename;
std::ofstream env_file;

// flow file (binary only)
std::string flow_b_filename;
std::ofstream flow_b_file;

// CellsSystem files
std::string cellsys_in_filename;
std::ifstream cellsys_in_file;

std::string cellsys_out_filename;
std::ofstream cellsys_out_file;

// nome del file CellType
std::string CellTypeFile;

// nome del file CellTypeAlt
std::string CellTypeFileAlt;

// nome del file Environment
std::string EnvironmentFile;

// nome del file di comando
std::string commandFile;


// *** dati per la gestione del sistema di cellule ***


// questo identifica il tipo di distribuzione iniziale delle cellule (nel caso di simulazione Full3D)
int initial_cell_dist;

// timers
#ifdef _OPENMP
double CPU_time, CPU_time_0, CPU_time_1, CPU_time_2;	// CPU time
#else
clock_t CPU_time, CPU_time_0, CPU_time_1, CPU_time_2;	// CPU time
#endif
double t_CPU, t_CPU_int, t_CPU_0, delta_t_CPU;		// CPU time (real)
//double t_CPU_max;			// tempo di CPU massimo per una singola frazione di run (in s)
time_t	time_now, time_old;		// tempo reale trascorso dall'ultimo step
double time_from_CGAL;		// tempo simulato trascorso dall'ultimo update della triangolazione di Delaunay
double timing;

// double dt;					// timestep
// double dt_sm;               // timestep slow motion
// double t;					// tempo attuale di simulazione
// double t_ini;				// tempo di inizializzazione
// double treal;				// tempo di simulazione vero (dopo l'inizializzazione)
// double tmax;				// tempo max di simulazione
// double tsm_start;           // tempo di start dello slow motion
// double tsm_stop;            // tempo di stop dello slow motion
// bool slow_motion;           // indica lo stato dello slow motion


bool ready2start;				// Global variable by which the simulation can be performed on a regular basis
bool faketumAtCurrentTime;


// statistiche
int ncalls;						// numero di volte che si inizializza il loop dell'algoritmo di diffusione
int	nrepeats;					// ultimo valore del numero di ripetioni dell'algoritmo di diffusione, trasporto, metabolismo
long int min_nrepeats;          // numero minimo di ripetizioni dell'algoritmo di diffusione: per la convergenza ci deve essere almeno questo numero minimo che viene ricalcolato ad ogni inizio di loop
double nrepeats_average;
long int nrepeats_max, nrepeats_min;
unsigned long convergence_fail[NCONV_TEST];	// variabile che tiene conto della mancata convergenza nell'algoritmo di diffusione

unsigned long loop_count;		// numero di passaggi nel loop dell'algoritmo di calcolo della posizione
double loop_count_average;
unsigned long loop_count_max, loop_count_min;

long int n_mitoses;				// numero di mitosi avvenute in un singolo passo
double n_mitoses_average;
long int n_mitoses_max, n_mitoses_min;

double min_Gin;
double max_Gin;
double min_Gext;
double max_Gext;
double min_O2;
double max_O2;
double min_Ain;
double max_Ain;
double min_Aext;
double max_Aext;
double min_AcLin;
double max_AcLin;
double min_AcLext;
double max_AcLext;
double min_extvolume;
double max_extvolume;


// variabili globali che definiscono l'output
std::string machine;					// nome della macchina
//int run;						// numero del run
std::string dir;						// nome del directory di output
int part;						// per un run suddiviso in parti eseguite in tempi diversi, questa variabile indica la parte attuale

// numero iniziale di cellule
unsigned long nstart;



// ultimo nome assegnato (serve a gestire in modo univoco l'assegnazione del nome alle cellule)
unsigned long lastname;

// ambiente iniziale
Environment Env_0;

// ambiente attuale
Environment Env;

// SabryNew, gestione eventi di mutazione
MutEventCreator MutationEv;

// variazione dello stato dell'ambiente 
//Environment delta_Env;

EnvironmentalSignal flowSignal;
    
// flusso dell'ossigeno (flag che indica il tipo del funzionamento nel modo bioreattore; se la flag è accesa allora viene modulato anche
// l'ossigeno, che in caso contrario resta invece fisso al valore ambientale); questa flag non serve a nulla se flowON e' spenta.
bool oxygenflowON;

EnvironmentalSignal dose_rateSignal;

// definizione del vettore dei fenotipi cellulari
std::vector<CellType> CellTypeVector;
std::vector<unsigned long> CellTypeIndexVector;


// parameters related to geometry and dynamics of the cluster
double maxdr;				// max displacement in the cellular system at a dynamic pace

// parametri metabolici del cluster
double O2Flow;				// flusso di O2 dalla periferia alle regioni interne (in kg/s)

double AcLFlow;			// flusso di AcL nell'ambiente (in kg/s)

// *** fine dei dati per la gestione del sistema di cellule ***

//****************************************************************************************************

// *** vettori di dati associati alle singole cellule ***


	// informazioni fondamentali sulla cellula: nome ed eventuale marcatura (che viene ereditata), 
	// fenotipo e caratteristiche generali dell'ambiente
		
	std::vector<unsigned long> name;			// identificatore della cellula 	
	std::vector<int> mark;					// label applicato alla cellula
	std::vector<CellType*> type;				// puntatore a un CellType (fenotipo cellulare)
	
	// temperatura della cellula
	
	std::vector<double> Temperature;	
		
	// stato cellulare (associato al metabolismo)
	
	std::vector<CellPhase> phase;			// fase cellulare
	std::vector<int> death_condition;		// variabile che registra la ragione della morte cellulare
	std::vector<float> age;					// eta' della cellula (tempo in secondi dalla nascita)
	std::vector<float> phase_age;			// eta' dello stato cellulare attuale (tempo in secondi dall'inizio)
	std::vector<float> age_mother;			// eta' della cellula madre (tempo in secondi dalla nascita alla mitosi)
	std::vector<int> n_mitosis;				// numero di mitosi dall'inizio della simulazione
	
	// Geometric and topological variables
	
	// seeding spheroid postion
  double New_x_0 = 0;
  double New_y_0 = 0;
  double New_z_0 = 0;

	std::vector<double> x;					// posizione nello spazio del centro della cellula
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> vx;					// velocità della cellula
	std::vector<double> vy;
	std::vector<double> vz;
	std::vector<double> r;					// raggio cellulare	
	std::vector<double> surface;             // cell surface
	std::vector<double> volume;				// volume cellulare
	std::vector<double> mass;				// massa cellulare
	
	// copies for use in the geometry calculation method
	std::vector<double> vxnew;				// velocità della cellula
	std::vector<double> vynew;
	std::vector<double> vznew;

	std::vector<double> fx;					// forze
	std::vector<double> fy;
	std::vector<double> fz;
  
  std::vector< std::pair<Point,unsigned> > v;	// vector of points with info passed to CGAL

	std::vector<double> volume_extra;		// volume of the extracellular region surrounding the cell

	//in random experiment max was 82
	static const int MAX_N_NEIGH = 122;
	std::vector<unsigned long> neigh;					// number of neighbors
	std::vector< std::array<unsigned long, MAX_N_NEIGH> > vneigh;		// vector of neighbors
	std::vector< std::array<double, MAX_N_NEIGH> > vdist;     // vector of distances from neighbors
	std::vector< std::array<double, MAX_N_NEIGH> > vcsurf;	// vector of contact surfaces with neighbors (approx calculation)
	std::vector< std::array<double, MAX_N_NEIGH> > gnk;		// vector of geometric factors
	std::vector<double> contact_surf;		// The total area of ​​the contact surface with adjacent cells
	
	std::vector<bool> isonCH;				// label indicating if the cell is on convex hull
	std::vector<bool> isonAS;				// label indicating if the cell is on the alpha shape
	std::vector<int> isonBV;					// this label is true if the cell is in contact
										// with a blood vessel
										// if a cell is not in contact with blood vessel, 
										// then isonBV[n] = 0, else isonBV[n] = position 
										// of vessel in BloodVesselVector + 1
										// this shift is done so that isonBV can also be 
										// used as a boolean variable 
										
	std::vector<double> env_surf;			// contact surface with the environment
	std::vector<double> g_env;				// The geometric factor relating to contact with the environment

	std::vector<double> bv_surf;				// contact surface with blood vessels
	std::vector<double> g_bv;				// The geometric factor relating to contact with blood vessels

	std::vector<double> M;					// numero di mitocondri 

	// lista delle variabili metaboliche all'interno della cellula
	
	std::vector<double> G;					// quantita' (massa) interna di glucosio
	std::vector<double> G6P;                 // quantita' (massa) interna di G6P
	std::vector<double> O2;					// quantita' (massa) interna di ossigeno
	std::vector<double> store;				// riserva in unita' di massa di glucosio equivalente
	std::vector<double> A;					// quantita' (massa) interna di altri nutrienti (glutammina)
	std::vector<double> AcL;                 // quantita' (massa) interna di acido lattico
    
  std::vector<double> h;                   // parametrizzazione dell'attivita' di trasporto del glucosio in funzione di O2
	std::vector<double> pHi;                 // valore del pH interno alla cellula
	// double H;                        // quantita' (massa) interna di protoni liberi
	// double CO2;                      // quantita' (massa) interna di anidride carbonica

	std::vector<double> protein;             // quantita' totale di proteine nelle cellule (in kg)
	std::vector<double> prot_rate;			// rate di produzione delle proteine
	std::vector<double> DNA;                 // quantita' di DNA prodotto per la replicazione 
                                        // (in molecole di DNA. Una molecola = intero genoma)
	std::vector<double> DNA_rate;			// rate di produzione del DNA

	// rates associate al glucosio, allo store e all'acido lattico
	
	std::vector<double> GAbsRate;			// rate di assorbimento del glucosio
	std::vector<double> GConsRate;			// rate di consumo del glucosio
	std::vector<double> AAbsRate;			// rate di assorbimento della glutammina
	std::vector<double> AConsRate;			// rate di consumo del glutammmina
	std::vector<double> StoreFillRate;		// rate di riempiemento dello store
	std::vector<double> StoreConsRate;		// rate di consumo dello store
	std::vector<double> AcLRate;				// rate di produzione dell'acido lattico
	std::vector<double> AcLOutRate;			// rate di trasporto all'esterno dell'acido lattico
	std::vector<double> O2Rate;				// O2 concentration rate of change (latest computed value)  //*************************** new for O2 rate ******* april 2018

	// calcolo dell'ATP (richiede parecchie variabili ausiliarie che hanno significato biologico)
	
	std::vector<double> ATP_St;				// ATP standard
	std::vector<double> ATP_Ox;				// rate di produzione di ATP tramite fosforilazione ossidativa
	std::vector<double> ATP_NOx;             // rate di produzione di ATP tramite glicolisi anaerobica
	std::vector<double> ATP2;				// rate di produzione di ATP dallo store
	std::vector<double> ATP3;				// rate di produzione di ATP dagli altri nutrienti
	std::vector<double> ConsATP;             // rate di consumo di ATP associato al metabolismo
	std::vector<double> ConsATP_1;			// rate di consumo di ATP associato ai mitocondri
	std::vector<double> ConsATP_2;			// rate di consumo di ATP associato alla produzione di proteine
	std::vector<double> ConsATP_3;			// rate di consumo di ATP associato alla produzione di DNA
	std::vector<double> ConsATP_4;			// rate di consumo di ATP associato alla pompa protonica
	std::vector<double> ConsATP_5;			// rate di consumo di ATP associato alla proliferazione dei mitocondri
	std::vector<double> ATPtot;				// rate totale di variazione dell'ATP = somma dei rates
	std::vector<double> ATPp;				// massa totale di ATP all'interno della cellula
	std::vector<double> ATPmin;				// valore minimo dell'ATPp (al di sotto di questo la cellula muore)

	// altre variabili correlate all'ATP
	
	std::vector<double> ATPstart;			// ATP iniziale nella cellula (in kg)
	std::vector<double> ATPprod;             // ATP prodotto dalla cellula (in kg) nel passo metabolico
	std::vector<double> ATPcons;             // ATP consumato dalla cellula (in kg) nel passo metabolico


	// variabili metaboliche negli spazi extracellulari
	
	std::vector<double> G_extra;             // massa del glucosio nella matrice extracellulare (negli spazi intercellulari)
	std::vector<double> A_extra;             // massa degli altri nutrienti (glutammina) nella matrice extracellulare (negli spazi intercellulari)
	std::vector<double> AcL_extra;			// massa dell'AcL nella matrice extracellulare (negli spazi intercellulari)
	
	std::vector<double> pH;					// valore del pH extracellulare
	// double H_extra;					// massa dei protoni liberi nella matrice extracellulare (negli spazi intercellulari)
	// double CO2_extra;				// massa dell'anidride carbonica nella matrice extracellulare (negli spazi intercellulari)

	// altre quantita' derivate
	
	std::vector<double> SensO2;				// frazione di ossigeno disponibile rispetto la richiesta
	std::vector<double> ConsO;				// rate di consumo dell'ossigeno (kg/s)
	// double ProdCO2;					// rate di produzione di CO2 (kg/s)
	
	// proteine e DNA

	std::vector<double> DNA_spread;			// variazione individuale della velocita' di sintesi del DNA 
                                        // (parte relativa al consumo di ATP)
                                        // questa variazione e' una modellizzazione grezza della velocita' 
                                        // variabile di sintesi dovuta ai 
                                        // danni al DNA (da modificare in futuro)
	
	std::vector<double> M_T;                 // duration of the media M (in s)
	std::vector<double> pRb;                 // quantita' (massa) della pRb

	std::vector<double> ConcS;				// concentrazione MOLARE della sostanza S che fa da substrato alla reazione di MM per la soglia

	std::vector<double> cyclinD;             // quantita' totale di ciclina D (in kg)
	std::vector<double> cyclinE;             // quantita' totale di ciclina E (in kg)
	std::vector<double> cyclinX;             // quantita' totale di cicline (A e B) attive in fase G2 (in kg)
	
	std::vector<double> NpRbk;				// numero di molecole di pRb attive

// copie per i calcoli in Diff

	std::vector<double> volumeOld;			// vettore dei volumi cellulari (valore vecchio)
	std::vector<double> volumeNew;			// vettore dei volumi cellulari (valore nuovo)
	std::vector<double> volume_extraOld;     // volume dei volume extracellulari (valore vecchio)
	std::vector<double> volume_extraNew;     // volume dei volume extracellulari (valore nuovo)
                                        // anche se il volume non e' una variabile dinamica, va immagazzinato in vettori 
                                        // per poter effettuare le somme necessarie al calcolo della diffusione
																
	std::vector<double> MitOld;				// vettore dei mitocondri (valore vecchio)
	std::vector<double> MitNew;				// vettore dei mitocondri (valore nuovo)

	std::vector<double> pHiOld;				// vettore dei pH cellulari (valore vecchio)
	std::vector<double> pHiNew;				// vettore dei pH cellulari (valore nuovo)
	std::vector<double> pHOld;				// vettore dei pH extracellulari (valore vecchio)
	std::vector<double> pHNew;				// vettore dei pH extracellulari (valore nuovo)
                                        // anche se il pH non e' una variabile dinamica, va immagazzinato in vettori 
                                        // per poter effettuare le somme necessarie al calcolo della diffusione
	
	std::vector<double> mGinOld;             // vettore della massa di glucosio dentro le cellule (valore vecchio)
	std::vector<double> mGinNew;             // vettore della massa di glucosio dentro le cellule (valore nuovo)
	std::vector<double> mGextOld;			// vettore della massa di glucosio negli spazi extracellulari (valore vecchio)
	std::vector<double> mGextNew;			// vettore della massa di glucosio negli spazi extracellulari (valore nuovo)

	std::vector<double> mG6POld;             // vettore della massa di G6P nelle cellule (valore vecchio)
	std::vector<double> mG6PNew;             // vettore della massa di G6P nelle cellule (valore nuovo)

	std::vector<double> mO2Old;				// vettore della massa di O2 nelle cellule (valore vecchio)
	std::vector<double> mO2New;				// vettore della massa di O2 nelle cellule (valore nuovo)

	std::vector<double> StoreOld;			// vettore dello Store nelle cellule (valore vecchio)
	std::vector<double> StoreNew;			// vettore dello Store nelle cellule (valore nuovo)

	std::vector<double> mAinOld;             // vettore della massa di glutammina dentro le cellule (valore vecchio)
	std::vector<double> mAinNew;             // vettore della massa di glutammina dentro le cellule (valore nuovo)
	std::vector<double> mAextOld;			// vettore della massa di glutammina negli spazi extracellulari (valore vecchio)
	std::vector<double> mAextNew;			// vettore della massa di glutammina negli spazi extracellulari (valore nuovo)

	std::vector<double> mAcLinOld;			// vettore della massa di acido lattico dentro le cellule (valore vecchio)
	std::vector<double> mAcLinNew;			// vettore della massa di acido lattico dentro le cellule (valore nuovo)
	std::vector<double> mAcLextOld;			// vettore della massa di acido lattico negli spazi extracellulari (valore vecchio)
	std::vector<double> mAcLextNew;			// vettore della massa di acido lattico negli spazi extracellulari (valore nuovo)

	std::vector<double> ATPpOld;             // vettore dell'ATPp nelle cellule (valore vecchio)
	std::vector<double> ATPpNew;             // vettore dell'ATPp nelle cellule (valore nuovo)

	// altre allocazioni di memoria utili per il metodo di diffusione
	
	std::vector<double> proteinNew;			// vettore della massa di proteina nelle cellule
	std::vector<double> pRbNew;				// vettore della massa di pRb nelle cellule
	std::vector<double> delta_protein;		// vettore della variazione di massa totale di proteine nelle cellule
	std::vector<double> ConcSNew;			// concentrazione del substrato S per il calcolo delle soglie dei checkpoints

	std::vector<double> DNANew;				// vettore che contiene la lunghezza relativa della molecola di DNA (1 = intera molecola)

	// vector<double> DNA_rate;			// vettore che contiene il rate di sintesi del DNA
	
	// vector<double> GAbsRate;			// vettore che memorizza l'assorbimento di glucosio
	// vector<double> GConsRate;		// vettore che memorizza il consumo di glucosio

	// vector<double> AAbsRate;			// vettore che memorizza l'assorbimento di glutammina
	// vector<double> AConsRate;		// vettore che memorizza il consumo di glutammina
	
	// vector<double> StoreFillRate;	// vettore che memorizza il rate di riempiemento dello store
	// vector<double> StoreConsRate;	// vettore che memorizza il rate di consumo dello store
	
	// vector<double> AcLRate;			// vettore che memorizza il rate di produzione dell'acido lattico
	// vector<double> AcLOutRate;		// vettore che memorizza il rate di espulsione dell'acido lattico
	
	// vector<double> ATP_Ox;			// vettore che memorizza il rate ATP_Ox
	// vector<double> ATP_NOx;			// vettore che memorizza il rate ATP_NOx
	// vector<double> ATP2;				// vettore che memorizza il rate ATP2
	// vector<double> ATP3;				// vettore che memorizza il rate ATP3
	// vector<double> ConsATP;			// vettore che memorizza il rate ConsATP
	// vector<double> ConsATP_1;		// vettore che memorizza il rate ConsATP_1
	// vector<double> ConsATP_2;		// vettore che memorizza il rate ConsATP_2
	// vector<double> ConsATP_3;		// vettore che memorizza il rate ConsATP_3
	// vector<double> ConsATP_5;		// vettore che memorizza il rate ConsATP_5

// *** fine dei vettori di dati associati alle singole cellule ***

// *** dati associati ai vasi sanguigni
    
    int nbv;                                // numero di vasi sanguigni
    // the vector structure messes with the memory, in fact with the stack!
    std::vector<BloodVessel> BloodVesselVector;  // vettore dei vasi sanguigni
    //BloodVessel BloodVesselVector;
    //boost::unordered_map<uint, vbl::BloodVessel> bloodVesselMap;
#if VBL_USE_TUMORCODE
    CellBasedO2Uptake o2_uptake_model;
#endif
// *** fine dei dati associati a vasi sanguigni


//****************************************************************************************************


public:
friend class ApplyGeometricCalculation;
ReadInParameters* get_params_pointer()
{
  return &params;
};

std::array<float,3> get_seeding_position();
//****************************************************************************************************

// *** Methods for managing the system ***

// Redeployment of the dynamic reserve
void Set_reserve(const int reserve); // cell vectors
void Set_BV_reserve(const int reserve); // blood vessel vector

int checkNeighbourhood_consistency(std::string atPlace);

// builder that builds a cell array of no length, but assigns a dynamic reserve to carriers
// Note that the creation of the CellsSystem also resets cell counts, types, blood vessels
//CellsSystem(const int reserve, const int reserve_bv) { ncells = 0; ntypes = 0; nbv=0; Set_reserve(reserve); Set_BV_reserve(reserve_bv); };
CellsSystem(const int reserve, const int reserve_bv) { params.ncells = 0; params.ntypes = 0; nbv=0; Set_reserve(reserve); };
// default builder, builds a set of cells of zero length and assigns the standard dynamic reserve to the carriers
CellsSystem() { params.ncells = 0; params.ntypes = 0; nbv=0; Set_reserve(RESERVE); Set_BV_reserve(RESERVE_BV); };
//CellsSystem() { ncells = 0; ntypes = 0; nbv=0; Set_reserve(RESERVE); };
// aggiunta di cellule non inizializzate al sistema
void AddCells( const int newcells );
// Adding a single standardized standardized cell
void AddInitializedCell(int& idum, CellType* cType, Environment* cEnv);
// inizializzazione standard dell'intero sistema di cellule
void InitCells( int& idum, CellType* cType, Environment* cEnv ) { for(unsigned long int k=0; k<params.ncells; k++) AddInitializedCell(idum, cType, cEnv ); };
// copia di una cellula (k-esima) in un'altra sezione del sistema di cellule (da kstart incluso a kstop incluso)
int CopyCell( const unsigned long int k, const unsigned long int kstart, const unsigned long int kstop);
// metodo di copia: replica la cellula k-esima inserendone una copia alla fine (tutto tranne le caratteristiche geometriche-topologiche)
// inoltre ridefinisce il celltype
int CopyCell( const unsigned long int k, const unsigned long int kstart, const unsigned long int kstop, CellType* newtype);
// metodo di copia: replica la cellula k-esima inserendone una copia alla fine (tutto tranne le caratteristiche geometriche-topologiche),
int ReplicateCell( const unsigned long int k );
// metodo di copia: replica la cellula k-esima inserendo n copie alla fine (tutto tranne le caratteristiche geometriche-topologiche)
int ReplicateCell( const unsigned long int k, const unsigned long int n );
// metodo di copia: replica la cellula k-esima inserendone una copia alla fine (tutto tranne le caratteristiche geometriche-topologiche),
// inoltre ridefinisce il celltype
int ReplicateCell( const unsigned long int k, CellType* newtype );

// Method of printing a single cell on a file
void PrintCell( std::ostream& stream, const unsigned long int k );
// Printing data in the form of a single string readable by a spreadsheet program
void PrintCellData( const unsigned long int k, std::ofstream& stream, long int nrec );

// rimozione di una cellula (il metodo gestisce tutti i vettori)
void RemoveCell( const unsigned long  n );

// altri metodi

// calcolo della forza tra le cellule nella configurazione attuale
void GetForces();
// calcolo della posizione e della velocità di ciascuna cellula (si calcola solo se ci sono almeno due cellule)
void NewPositionsAndVelocities( );


// getters
std::vector<unsigned long int> get_CellTypeIndexVector();
void set_CellTypeFromIndexVector(std::vector<unsigned long> &cellIndexVector);
void set_CellPhaseFromIntVector(std::vector<int> &int_buffer);

std::string Get_Commands( ) { return commandFile; };
std::string Get_CellTypeFile( ) { return CellTypeFile; };
std::string Get_CellTypeFileAlt( ) { return CellTypeFileAlt; };
std::string Get_EnvironmentFile( ) { return EnvironmentFile; };
int Get_sim_type() { return params.sim_type; };
int Get_initial_cell_dist() { return initial_cell_dist; };
double Get_dt();
double Get_dt_sm();
bool Get_slow_motion();
double Get_t();
double Get_treal();
double Get_tmax();
double Get_time_from_CGAL();

unsigned long Get_ncalls() { return ncalls; };
unsigned long Get_nstep();
unsigned long Get_nstep_start();
unsigned long Get_nmax();

int Get_idum();

bool Get_ready2start() { return ready2start; };
unsigned long Get_nprint();
unsigned long Get_nscreen(); 
unsigned long Get_nconfiguration();

double Get_eps();
double Get_delta_vmax();
double Get_nrepeats() { return nrepeats; };
long int Get_min_nrepeats() { return min_nrepeats; };
double Get_nrepeats_average() { return nrepeats_average; };
double Get_nrepeats_max() { return nrepeats_max; };
double Get_nrepeats_min() { return nrepeats_min; };

double Get_max_x();
double Get_max_y();
double Get_max_z();
double Get_min_x();
double Get_min_y();
double Get_min_z();

unsigned long Get_nstart() { return nstart; };
unsigned long Get_ncells();
unsigned long Get_alive();
unsigned long Get_ntypes();
unsigned long Get_lastname() { return lastname; };
bool Get_flowON();
bool Get_doseON();
Environment Get_Env() { return Env; };
Environment Get_Env_0() { return Env_0; };
EnvironmentalSignal Get_flowSignal() { return flowSignal; };
EnvironmentalSignal Get_dose_rateSignal() { return dose_rateSignal; };
double Get_maxdr() { return maxdr; };
double Get_O2Flow() { return O2Flow; };
double Get_AcLFlow() { return AcLFlow; };

// setters
void Set_Commands( const std::string newcommandFile ) 
{ 
  commandFile = newcommandFile; 
};
void Set_CellTypeFile( std::string newCellTypeFile ) 
{ 
  CellTypeFile = newCellTypeFile; 
};
void Set_CellTypeFileAlt( std::string newCellTypeFile ) 
{ 
  CellTypeFileAlt = newCellTypeFile; 
};
void Set_EnvironmentFile( std::string newEnvironmentFile ) 
{ 
  EnvironmentFile = newEnvironmentFile; 
};
void Set_idum ( int newidum );
void Set_dt( double newdt );
void Set_time_from_CGAL( double newtime_from_CGAL ) { time_from_CGAL = newtime_from_CGAL; };
void Set_ready2start( bool newr2s ) { ready2start = newr2s; };
void Set_nconfiguration( unsigned long newnconfiguration );
void Step_nconfiguration();
void Set_eps( double neweps );
void Set_delta_vmax( double newdelta_vmax );



// inizializzazione del run da terminale
void InitializeCellsSystem( bool terminal );

// inizializzazione del run da file
void InitializeCellsSystem( const std::string filename );

// metodo per gli eventi cellulari (crescita, mitosi, etc.; restituisce una flag che indica se c'e' stata almeno una mitosi nel sistema di cellule)
//bool CellEvents( );
int CellEvents( );

// metodo per la pulizia della memoria (al momento usa operazioni non molto efficienti)
void CleanCellsSystem( );

// geometria
void Geometry();
#ifdef serial
void Geometry_serial(); //needed for the case of new cells
#endif
void NoGeometry();	// versione ridotta, calcoli minimi per cellule disperse

// meccanica del cluster
void Dynamics( );
void DummyDynamics( ); // dinamica dummy

// summary printout
void Printout();

// close files
void CloseOutputFiles() 
{ 
  std::cout << "\nFine del run " << params.run << " al passo " << params.nstep << std::endl; 
  output_file.close(); 
  log_file.close(); 
  screen_dump_file.close();
  errorlog_file.close();
  cell_file.close();
  env_file.close();
  convlog_file.close();
  std::cout << "\nclosed all vbl files "<< std::endl; 
};

// summary printout on file
void Print2file();

// singola stampa su logfile
void Print2logfile(std::string str);			// cellula ben formattata
void PrintAll2logfile(std::string str);		// stampa TUTTE le cellule su logfile (!!!)

// definizione del run
void RunDefinition( );
void RunDefinition( std::string run_name );

// scrittura e lettura del CellsSystem
void WriteCellsSystem( );
void ReadCellsSystem( );


// avanzamento dei timers
//bool TimersAdvance( );
// if no endtime is present, the behaviour of old TimersAdvance is restored

/* T.F. 06.02.2019 
 * obviously there are problems with different compilers and boost::optional,
 * therefore I changed the behaviour: endtime == 0.0 means unlimited,
 * otherwise run until endtime
 */
bool TimersAdvanceUntil( double &endtime );

// CPU timing
void CPU_timer( timer_button button );

// Timing
double Timing( bool reset );

// distanza tra cellule con indici j e k
inline double Distance(const int j, const int k)
{
  return sqrt( SQR(x[j]-x[k]) + SQR(y[j]-y[k]) + SQR(z[j]-z[k]) );
};

// printout dell'header di una singola configurazione
void PrintHeader(bool isBinary);

// printout dei punti su file 
void PrintPoints(bool isBinary);

// questo metodo  stampa solo i links su file
void PrintLinks(bool isBinary);

// questo metodo  stampa solo la flag di appartenenza al CH su file
void PrintCHFlag(bool isBinary);

// questo metodo  stampa solo la flag di appartenenza all'AS su file
void PrintASFlag(bool isBinary);

// this method prints on file the isonBV flag
void PrintBVFlag(bool isBinary);

// This method prints only the cellular rays on files
void PrintR(bool isBinary);

// questo metodo stampa il codice della fase cellulare
void PrintPhase(bool isBinary);

// questo metodo stampa il codice della fase cellulare
void PrintType(bool isBinary);
    
// questo metodo stampa parametri importanti dentro le cellule
void PrintVar(bool isBinary);

// questo metodo stampa i dati ambientali essenziali
void PrintEnv(bool isBinary);

// questo metodo scrive un record di configurazione
void PrintConfiguration(bool isBinary);

// questo metodo calcola e stampa i flussi extracellulari
void PrintFlows();

// metodo per il calcolo della diffusione e metodi collegati
void Diff();

// metodo per il calcolo delle statistiche, passo per passo
void StepStat( bool reset_stat );

// metodo per la stampa delle log files dettagliate
// void PrintLog(long int n) { PrintCellData(cell_file,nstep); Env.PrintEnvironmentData(env_file,nstep); };

// *** fine della parte dei metodi per la gestione del sistema ***

//****************************************************************************************************

unsigned int runMainLoop( double &endtime);
//unsigned int runMainLoop( );
#if VBL_USE_TUMORCODE
  //void Set_Tumorcode_Continuous_lattice(LatticeDataQuad3d &field_ld);
void Set_Tumorcode_O2_uptake_model(CellBasedO2Uptake  &p_o2_uptake_model)
{
  std::cout << "setting to " << &p_o2_uptake_model << std::endl;
  o2_uptake_model = p_o2_uptake_model;
  std::cout << "after setting to " << &o2_uptake_model << std::endl;
}

//void interpolate_O2_uptake_to_tumorcode_2(CellBasedO2Uptake &o2_uptake_model, std::vector<double> &O2Rates);
#endif
// *** metodi per la gestione della parte biofisica *** 

	

	// setters (nella forma di vettore e di singolo elemento)

	void Set_name(const std::vector<unsigned long>& newname) { name = newname; };
	void Set_name(const int k, const unsigned long newname) { name[k] = newname; };
	
	void Set_mark(const std::vector<int>& newmark) { mark = newmark; };
	void Set_mark(const int k, const int newmark) { mark[k] = newmark; };
	
	void Set_type(const std::vector<CellType*>& cType) { type = cType; };
	void Set_type(const int k, CellType* cType) { type[k] = cType; };
	
	void Set_Temperature(const std::vector<double>& newTemperature) { Temperature = newTemperature; };
	void Set_Temperature(const int k, const double newTemperature) { Temperature[k] = newTemperature; };
	
	void Set_phase(const std::vector<CellPhase>& newphase) { phase = newphase; };
	void Set_phase(const int k, const CellPhase newphase) { phase[k] = newphase; };
	
	void Set_death_condition( const std::vector<int>& newdeath_condition ) { death_condition = newdeath_condition; };
	void Set_death_condition( const int k, const int newdeath_condition ) { death_condition[k] = newdeath_condition; };

	void Set_age(const std::vector<float>& newage) { age = newage; };
	void Set_age(const int k, const float newage) { age[k] = newage; };

	void Set_phase_age(const std::vector<float>& newphase_age) { phase_age = newphase_age; };
	void Set_phase_age(const int k, const float newphase_age) { phase_age[k] = newphase_age; };

	void Set_age_mother(const std::vector<float>& newage_mother) { age_mother = newage_mother; };
	void Set_age_mother(const int k, const float newage_mother) { age_mother[k] = newage_mother; };

	void Set_n_mitosis(const std::vector<int>& newn_mitosis ) { n_mitosis = newn_mitosis; };
	void Set_n_mitosis(const int k, const int newn_mitosis ) { n_mitosis[k] = newn_mitosis; };

	void Add_BloodVesselVector( BloodVessel NewBV );
  //BloodVessel& Get_BloodVessel_at(int index){ return BloodVesselVector[index];}
  void clean_BloodVesselVector();
// 	{
// 	  nbv++;
//   //Vector()[index] = NewBV;
//   //BloodVesselVector.push_back(NewBV);
//   BloodVesselVector[index] = vbl::BloodVessel(*NewBV); /*cout << "New blood vessel in CellsSystem" << endl;*/
//   //bloodVesselMap[index] = NewBV;
// 	}
// 	; 
    
    // geometria, cinematica e dinamica

	void Set_x( const std::vector<double>& newx ) { x = newx; };
	void Set_x( const int k, const double newx ) { x[k] = newx; };

	void Set_y( const std::vector<double>& newy ) { y = newy; };
	void Set_y( const int k, const double newy ) { y[k] = newy; };

	void Set_z( const std::vector<double>& newz ) { z = newz; };
	void Set_z( const int k, const double newz ) { z[k] = newz; };

	void Set_vx( const std::vector<double>& newvx ) { vx = newvx; };
	void Set_vx( const int k, const double newvx ) { vx[k] = newvx; };

	void Set_vy( const std::vector<double>& newvy ) { vy = newvy; };
	void Set_vy( const int k, const double newvy ) { vy[k] = newvy; };

	void Set_vz( const std::vector<double>& newvz ) { vz = newvz; };
	void Set_vz( const int k, const double newvz ) { vz[k] = newvz; };
  
  void Set_surface( const std::vector<double>& newsurface ) { surface = newsurface; };
  void Set_mass( const std::vector<double>& newmass ) { mass = newmass; };
	
	void Set_r( const std::vector<double>& newr ) { r = newr; };
	void Set_volume( const std::vector<double>& newvolume ) { volume = newvolume; };
	void Set_M( const std::vector<double>& newM ) { M = newM; };
	void Set_ATPp( const std::vector<double>& newATPp ) { ATPp = newATPp; };
			
	void Init_fx( ) { fx.clear(); };
	void Init_fy( ) { fy.clear(); };
	void Init_fz( ) { fz.clear(); };
	
	
	// fEnd of useful setters to the geometric part

	void Set_volume_extra( const std::vector<double>& newvolume_extra ) { volume_extra = newvolume_extra; };
	void Set_volume_extra( const int k, const double newvolume_extra ) { volume_extra[k] = newvolume_extra; };


	//void Set_neigh( const std::vector<int>& neighin ) {    neigh = neighin; };
	void Set_neigh( const int k, const int neighin ) { neigh[k] = neighin; };

	// questi setters esistono solo nella forma per singole cellule (inseriscono vettori di lunghezza variabile)
// 	void Set_vneigh( const int k, int* vneighin ) { vneigh[k].clear(); vneigh[k].insert( vneigh[k].begin(), vneighin, vneighin+neigh[k]); };
// 	void Set_vdist( const int k, double* vdistin ) { vdist[k].clear(); vdist[k].insert( vdist[k].begin(), vdistin, vdistin+neigh[k]); };
// 	void Set_vcsurf( const int k, double* vcsurfin ) { vcsurf[k].clear(); vcsurf[k].insert( vcsurf[k].begin(), vcsurfin, vcsurfin+neigh[k]); };
// 	void Set_gnk( const int k, double* newgnk ) { gnk[k].clear(); gnk[k].insert( gnk[k].begin(), newgnk, newgnk+neigh[k]); };
	// fine della parte dei setters non standard
  
  // use this one if we already know, that it is safe!
  void Set_vneigh_quick( const int k, const int kk, unsigned long vneighin ) { vneigh[k][kk]= vneighin; };
	void Set_vdist_quick( const int k,const int kk, double vdistin ) { vdist[k][kk]= vdistin; };
	void Set_vcsurf_quick( const int k,const int kk, double vcsurfin ) { vcsurf[k][kk] = vcsurfin; };
	void Set_gnk_quick( const int k,const int kk, double newgnk ) { gnk[k][kk]= newgnk; };
	// fine della parte dei setters non standard

	
	
	void Set_contact_surf( const std::vector<double>& newcontact_surf ) { contact_surf = newcontact_surf; };
	void Set_contact_surf( const int k, const double newcontact_surf ) { contact_surf[k] = newcontact_surf; };

	void Set_isonCH( const std::vector<bool>& isonCHnow ) { isonCH = isonCHnow; };
	void Set_isonCH( const int k, const bool isonCHnow ) 
  { 
    isonCH[k] = isonCHnow; 
  };

	void Set_isonAS( const std::vector<bool>& isonASnow ) { isonAS = isonASnow; };
	void Set_isonAS( const int k, const bool isonASnow ) { isonAS[k] = isonASnow; };

	void Set_isonBV( const std::vector<int>& isonBVnow ) { isonBV = isonBVnow; };
	void Set_isonBV( const int k, const int isonBVnow ) { isonBV[k] = isonBVnow; };

	void Set_env_surf( const std::vector<double>& newenv_surf ) { env_surf = newenv_surf; };
	void Set_env_surf( const int k, const double newenv_surf ) { env_surf[k] = newenv_surf; };

	void Set_g_env( const std::vector<double>& newg_env ) { g_env = newg_env; };
	void Set_g_env( const int k, const double newg_env ) { g_env[k] = newg_env; };

	void Set_bv_surf( const std::vector<double>& newbv_surf ) { bv_surf = newbv_surf; };
	void Set_bv_surf( const int k, const double newbv_surf ) { bv_surf[k] = newbv_surf; };

	void Set_g_bv( const std::vector<double>& newg_bv ) { g_bv = newg_bv; };
	void Set_g_bv( const int k, const double newg_bv ) { g_bv[k] = newg_bv; };


	void Set_G( const std::vector<double>& newG ) { G = newG; };
	void Set_G( const int k, const double newG ) { G[k] = newG; };

	void Set_G6P( const std::vector<double>& newG6P ) { G6P = newG6P; };
	void Set_G6P( const int k, const double newG6P ) { G6P[k] = newG6P; };

	void Set_O2( const std::vector<double>& newO2 ) { O2 = newO2; };
	void Set_O2( const int k, const double newO2 ) { O2[k] = newO2; };
  
  void Set_SensO2( const std::vector<double>& newSensO2 ) { SensO2 = newSensO2; };
	void Set_SensO2( const int k, const double newSensO2 ) { SensO2[k] = newSensO2; };
  void Set_ConsO( const std::vector<double>& newConsO ) { ConsO = newConsO; };
	void Set_ConsO( const int k, const double newConsO ) { ConsO[k] = newConsO; };

	void Set_store( const std::vector<double>& newstore ) { store = newstore; };
	void Set_store( const int k, const double newstore ) { store[k] = newstore; };

	void Set_A( const std::vector<double>& newA ) { A = newA; };
	void Set_A( const int k, const double newA ) { A[k] = newA; };

	void Set_AcL( const std::vector<double>& newAcL ) { AcL = newAcL; };
	void Set_AcL( const int k, const double newAcL ) { AcL[k] = newAcL; };
	
	void Set_pHi( const std::vector<double>& newpHi ) 
  { pHi = newpHi; };
	void Set_pHi( const int k, const double newpHi ) 
  { pHi[k] = newpHi; };
  
  void Set_h( const std::vector<double>& newph ) 
  { h = newph; };

	// void Set_H( const double newH ) { H = newH; };
	// void Set_CO2( const double newCO2 ) { CO2 = newCO2; };
	
	
	void Set_G_extra( const std::vector<double>& newG_extra ) { G_extra = newG_extra; };
	void Set_G_extra( const int k, const double newG_extra ) { G_extra[k] = newG_extra; };

	void Set_A_extra( const std::vector<double>& newA_extra ) { A_extra = newA_extra; };
	void Set_A_extra( const int k, const double newA_extra ) { A_extra[k] = newA_extra; };

	void Set_AcL_extra( const std::vector<double>& newAcL_extra ) { AcL_extra = newAcL_extra; };
	void Set_AcL_extra( const int k, const double newAcL_extra ) { AcL_extra[k] = newAcL_extra; };


	void Set_pH( const std::vector<double>& newpH ) { pH = newpH; };
	void Set_pH( const int k, const double newpH ) { pH[k] = newpH; };
	
	
	void Set_protein( const std::vector<double>& newprot ) { protein = newprot; }; 
	void Set_protein( const int k, const double newprot ) { protein[k] = newprot; }; 

	void Set_prot_rate( const std::vector<double>& newprot_rate ) { prot_rate = newprot_rate; };
	void Set_prot_rate( const int k, const double newprot_rate ) { prot_rate[k] = newprot_rate; };

	void Set_DNA_rate( const std::vector<double>& newDNA_rate ) { DNA_rate = newDNA_rate; }; 
	void Set_DNA_rate( const int k, const double newDNA_rate ) { DNA_rate[k] = newDNA_rate; }; 

	void Set_pRb ( const std::vector<double>& newpRb ) { pRb = newpRb; };
	void Set_pRb ( const int k, const double newpRb ) { pRb[k] = newpRb; };

	void Set_cyclinD( const std::vector<double>& newcyclinD ) { cyclinD = newcyclinD; };
	void Set_cyclinD( const int k, const double newcyclinD ) { cyclinD[k] = newcyclinD; };

	void Set_cyclinE( const std::vector<double>& newcyclinE ) { cyclinE = newcyclinE; };
	void Set_cyclinE( const int k, const double newcyclinE ) { cyclinE[k] = newcyclinE; };

	void Set_cyclinX( const std::vector<double>& newcyclinX ) { cyclinX = newcyclinX; };
	void Set_cyclinX( const int k, const double newcyclinX ) { cyclinX[k] = newcyclinX; };

	void Set_ConcS( const std::vector<double>& newconcS ) { ConcS = newconcS; };
	void Set_ConcS( const int k, const double newconcS ) { ConcS[k] = newconcS; };
  void Set_NpRbk( const std::vector<double>& newNpRbk ) { NpRbk = newNpRbk; };
	void Set_NpRbk( const int k, const double newNpRbk ) { NpRbk[k] = newNpRbk; };

	
	void Set_DNA( const std::vector<double>& newDNA ) { DNA = newDNA; };
	void Set_DNA( const int k, const double newDNA ) { DNA[k] = newDNA; };

	void Set_DNA_spread( const std::vector<double>& newDNA_spread ) { DNA_spread = newDNA_spread; };
	void Set_DNA_spread( const int k, const double newDNA_spread ) { DNA_spread[k] = newDNA_spread; };

	
	void Set_M_T( const std::vector<double>& newM_T ) { M_T = newM_T; };
	void Set_M_T( const int k, const double newM_T ) { M_T[k] = newM_T; };
	

	void Set_GAbsRate( const std::vector<double>& newGAbsRate ) { GAbsRate = newGAbsRate; }; 
	void Set_GAbsRate( const int k, const double newGAbsRate ) { GAbsRate[k] = newGAbsRate; }; 

	void Set_GConsRate( const std::vector<double>& newGConsRate ) { GConsRate = newGConsRate; };
	void Set_GConsRate( const int k, const double newGConsRate ) { GConsRate[k] = newGConsRate; };

	void Set_AAbsRate( const std::vector<double>& newAAbsRate ) { AAbsRate = newAAbsRate; }; 
	void Set_AAbsRate( const int k, const double newAAbsRate ) { AAbsRate[k] = newAAbsRate; }; 

	void Set_AConsRate( const std::vector<double>& newAConsRate ) { AConsRate = newAConsRate; };
	void Set_AConsRate( const int k, const double newAConsRate ) { AConsRate[k] = newAConsRate; };

	void Set_StoreFillRate( const std::vector<double>& newStoreFillRate ) { StoreFillRate = newStoreFillRate; };
	void Set_StoreFillRate( const int k, const double newStoreFillRate ) { StoreFillRate[k] = newStoreFillRate; };

	void Set_StoreConsRate( const std::vector<double>& newStoreConsRate ) { StoreConsRate = newStoreConsRate; };
	void Set_StoreConsRate( const int k, const double newStoreConsRate ) { StoreConsRate[k] = newStoreConsRate; };

	void Set_AcLRate( const std::vector<double>& newAcLRate ) { AcLRate = newAcLRate; };
	void Set_AcLRate( const int k, const double newAcLRate ) { AcLRate[k] = newAcLRate; };

	void Set_AcLOutRate( const std::vector<double>& newAcLOutRate ) { AcLOutRate = newAcLOutRate; };
	void Set_AcLOutRate( const int k, const double newAcLOutRate ) { AcLOutRate[k] = newAcLOutRate; };
	void Set_O2Rate( const std::vector<double>& newO2Rate ) { O2Rate = newO2Rate; }; //*************************** new for O2 rate ******* april 2018


	void Set_ATP_Ox( const std::vector<double>& newATP_Ox ) { ATP_Ox = newATP_Ox; };
	void Set_ATP_Ox( const int k, const double newATP_Ox ) { ATP_Ox[k] = newATP_Ox; };

	void Set_ATP_NOx( const std::vector<double>& newATP_NOx ) { ATP_NOx = newATP_NOx; };
	void Set_ATP_NOx( const int k, const double newATP_NOx ) { ATP_NOx[k] = newATP_NOx; };

	void Set_ATP2( const std::vector<double>& newATP2 ) { ATP2 = newATP2; };
	void Set_ATP2( const int k, const double newATP2 ) { ATP2[k] = newATP2; };

	void Set_ATP3( const std::vector<double>& newATP3 ) { ATP3 = newATP3; };
	void Set_ATP3( const int k, const double newATP3 ) { ATP3[k] = newATP3; };

	void Set_ConsATP( const std::vector<double>& newConsATP ) { ConsATP = newConsATP; };
	void Set_ConsATP( const int k, const double newConsATP ) { ConsATP[k] = newConsATP; };

	void Set_ConsATP_1( const std::vector<double>& newConsATP_1 ) { ConsATP_1 = newConsATP_1; };
	void Set_ConsATP_1( const int k, const double newConsATP_1 ) { ConsATP_1[k] = newConsATP_1; };

	void Set_ConsATP_2( const std::vector<double>& newConsATP_2 ) { ConsATP_2 = newConsATP_2; };
	void Set_ConsATP_2( const int k, const double newConsATP_2 ) { ConsATP_2[k] = newConsATP_2; };

	void Set_ConsATP_3( const std::vector<double>& newConsATP_3 ) { ConsATP_3 = newConsATP_3; };
	void Set_ConsATP_3( const int k, const double newConsATP_3 ) { ConsATP_3[k] = newConsATP_3; };

	void Set_ConsATP_4( const std::vector<double>& newConsATP_4 ) { ConsATP_4 = newConsATP_4; };
	void Set_ConsATP_4( const int k, const double newConsATP_4 ) { ConsATP_4[k] = newConsATP_4; };

	void Set_ConsATP_5( const std::vector<double>& newConsATP_5 ) { ConsATP_5 = newConsATP_5; };
	void Set_ConsATP_5( const int k, const double newConsATP_5 ) { ConsATP_5[k] = newConsATP_5; };

	void Set_ATPtot( const std::vector<double>& newATPtot ) { ATPtot = newATPtot; };
	void Set_ATPtot( const int k, const double newATPtot ) { ATPtot[k] = newATPtot; };
  void Set_ATP_St( const std::vector<double>& newATP_St ) { ATP_St = newATP_St; };
	void Set_ATP_St( const int k, const double newATP_St ) { ATP_St[k] = newATP_St; };
  void Set_ATPmin( const std::vector<double>& newATPmin ) { ATPmin = newATPmin; };
	void Set_ATPmin( const int k, const double newATPmin ) { ATPmin[k] = newATPmin; };



	// setters che inizializzano i contatori di produzione dell'ATP (normalmente vanno chiamati solo al momento della mitosi)
	void Set_ATPstart(  ) { for(unsigned long int k=0; k<params.ncells; k++) ATPstart[k] = (double)ATPp[k]; };	
	void Set_ATPstart( const unsigned long int k ) { ATPstart[k] = (double)ATPp[k]; };
  void Set_ATPstart( const std::vector<double>& newATPstart ) { ATPstart = newATPstart; };

	void Set_ATPprod( const std::vector<double>& newATPprod ) { ATPprod = newATPprod; };
	void Set_ATPprod( const unsigned long int k, const double newATPprod ) { ATPprod[k] = newATPprod; };

	void Set_ATPcons( const std::vector<double>& newATPcons ) { ATPcons = newATPcons; };
	void Set_ATPcons( const unsigned long int k, const double newATPcons ) { ATPcons[k] = newATPcons; };
	



	// getters
	
	std::vector<unsigned long> Get_name() { return name; };
	unsigned long Get_name( const unsigned long int k ) { return name[k]; };

	std::vector<int> Get_mark() { return mark; };
	int Get_mark( const unsigned long int k ) { return mark[k]; };
	
	std::vector<CellType*> Get_type() { return type; };
	CellType* Get_type( const unsigned long int k ) { return type[k]; };
	
	std::vector<double> Get_Temperature() { return Temperature; };
	double Get_Temperature( const unsigned long int k ) { return Temperature[k]; };
	
	std::vector<CellPhase> Get_phase() { return phase; };	
  std::vector<int> Get_phase_int();
	CellPhase Get_phase( const unsigned long int k ) { return phase[k]; };	

	std::vector<int> Get_death_condition() { return death_condition; };
	int Get_death_condition( const unsigned long int k ) { return death_condition[k]; };

	std::vector<float> Get_age() { return age; };
	float Get_age( const unsigned long int k ) { return age[k]; };

	std::vector<float> Get_phase_age() { return phase_age; };
	float Get_phase_age( const unsigned long int k ) { return phase_age[k]; };

	std::vector<float> Get_age_mother() { return age_mother; };
	float Get_age_mother( const unsigned long int k ) { return age_mother[k]; };

	std::vector<int> Get_n_mitosis() { return n_mitosis; };
	int Get_n_mitosis( const unsigned long int k ) { return n_mitosis[k]; };

	std::vector<BloodVessel> Get_BloodVesselVector() { return BloodVesselVector; };
	BloodVessel* Get_Pointer_to_Vessel_at(const int k) { return &BloodVesselVector[k];}
	//boost::unordered_map<uint, vbl::BloodVessel> Get_BloodVesselMap() { return bloodVesselMap; };
	//vbl::BloodVessel* Get_BloodVessel(int k) { return &bloodVesselMap[k];}
	int Get_nbv() { return nbv; };

	std::vector<double> Get_x() { return x; };
	double Get_x( const unsigned long int k ) { return x[k]; };

	std::vector<double> Get_y() { return y; };
	double Get_y( const unsigned long int k ) { return y[k]; };

	std::vector<double> Get_z() { return z; };
	double Get_z( const unsigned long int k ) { return z[k]; };

	std::vector<double> Get_vx() { return vx; };
	double Get_vx( const unsigned long int k ) { return vx[k]; };

	std::vector<double> Get_vy() { return vy; };
	double Get_vy( const unsigned long int k ) { return vy[k]; };

	std::vector<double> Get_vz() { return vz; };
	double Get_vz( const unsigned long int k ) { return vz[k]; };

	std::vector<double> Get_r() { return r; };
	double Get_r( const unsigned long int k ) { return r[k]; };

	std::vector<double> Get_surface() { return surface; };
	double Get_surface( const unsigned long int k ) { return surface[k]; };

	std::vector<double> Get_volume() { return volume; };
	double Get_volume( const unsigned long int k ) { return volume[k]; };

	std::vector<double> Get_mass() { return mass; };
	double Get_mass( const unsigned long int k ) { return mass[k]; };


	std::vector<double> Get_volume_extra() { return volume_extra; };
	double Get_volume_extra( const unsigned long int k ) { return volume_extra[k]; };



	std::vector<unsigned long> Get_neigh() { return neigh; };
	int Get_neigh( const unsigned long int k ) { return neigh[k]; };

	// questi getters esistono solo nella forma per singole cellule (restituiscono vettori di lunghezza variabile)
// 	std::vector<int> Get_vneigh( const unsigned long int k ) { return vneigh[k]; };
// 	std::vector<double> Get_vdist( const unsigned long int k ) { return vdist[k]; };
// 	std::vector<double> Get_vcsurf( const unsigned long int k ) { return vcsurf[k]; };
// 	std::vector<double> Get_gnk( const unsigned long int k ) { return gnk[k]; };
	// fine della parte dei getters non standard

	std::vector<double> Get_contact_surf() { return contact_surf; };
	double Get_contact_surf( const unsigned long int k ) { return contact_surf[k]; };
	
	std::vector<bool> Get_isonCH() { return isonCH; };
	bool Get_isonCH( const unsigned long int k ) { return isonCH[k]; };

	std::vector<bool> Get_isonAS() { return isonAS; };
	bool Get_isonAS( const unsigned long int k ) { return isonAS[k]; };

	std::vector<int> Get_isonBV() { return isonBV; };
	int Get_isonBV( const unsigned long int k ) { return isonBV[k]; };

	std::vector<double> Get_env_surf() { return env_surf; }; 
	double Get_env_surf( const unsigned long int k ) { return env_surf[k]; }; 

	std::vector<double> Get_g_env() { return g_env; }; 
	double Get_g_env( const unsigned long int k ) { return g_env[k]; }; 

	std::vector<double> Get_bv_surf() { return bv_surf; }; 
	double Get_bv_surf( const unsigned long int k ) { return bv_surf[k]; }; 

	std::vector<double> Get_g_bv() { return g_bv; }; 
	double Get_g_bv( const unsigned long int k ) { return g_bv[k]; }; 


	std::vector<double> Get_M() { return M; }; 
	double Get_M( const unsigned long int k ) { return M[k]; }; 


	std::vector<double> Get_G() { return G; };
	double Get_G( const unsigned long int k ) { return G[k]; };

	std::vector<double> Get_G6P() { return G6P; };
	double Get_G6P( const unsigned long int k ) { return G6P[k]; };

	std::vector<double> Get_O2() { return O2; };
	double Get_O2( const unsigned long int k ) { return O2[k]; };

	std::vector<double> Get_store() { return store; };
	double Get_store( const unsigned long int k ) { return store[k]; };

	std::vector<double> Get_A() { return A; };
	double Get_A( const unsigned long int k ) { return A[k]; };

	std::vector<double> Get_AcL() { return AcL; };
	double Get_AcL( const unsigned long int k ) { return AcL[k]; };


	std::vector<double> Get_h() { return h; };
	double Get_h( const unsigned long int k ) { return h[k]; };
	
	std::vector<double> Get_pHi() { return pHi; };
	double Get_pHi( const unsigned long int k ) { return pHi[k]; };

	// double Get_H() { return H; };
	// double Get_CO2() { return CO2; };

	std::vector<double> Get_GAbsRate() { return GAbsRate; };
	double Get_GAbsRate( const unsigned long int k ) { return GAbsRate[k]; };

	std::vector<double> Get_GConsRate() { return GConsRate; };
	double Get_GConsRate( const unsigned long int k ) { return GConsRate[k]; };

	std::vector<double> Get_AAbsRate() { return AAbsRate; };
	double Get_AAbsRate( const unsigned long int k ) { return AAbsRate[k]; };

	std::vector<double> Get_AConsRate() { return AConsRate; };
	double Get_AConsRate( const unsigned long int k ) { return AConsRate[k]; };

	std::vector<double> Get_StoreFillRate() { return StoreFillRate; };
	double Get_StoreFillRate( const unsigned long int k ) { return StoreFillRate[k]; };

	std::vector<double> Get_StoreConsRate() { return StoreConsRate; };
	double Get_StoreConsRate( const unsigned long int k ) { return StoreConsRate[k]; };

	std::vector<double> Get_AcLRate() { return AcLRate; };
	double Get_AcLRate( const unsigned long int k ) { return AcLRate[k]; };

	std::vector<double> Get_AcLOutRate() { return AcLOutRate; };
	double Get_AcLOutRate( const unsigned long int k ) { return AcLOutRate[k]; };

	//*************************** new for O2 rate ******* april 2018
	std::vector<double> Get_O2Rate() { return O2Rate; }; 
	double Get_O2Rate( const unsigned long int k ) { return O2Rate[k]; };

	std::vector<double> Get_ATP_St() { return ATP_St; };
	double Get_ATP_St( const unsigned long int k ) { return ATP_St[k]; };

	std::vector<double> Get_ATP_Ox() { return ATP_Ox; };
	double Get_ATP_Ox( const unsigned long int k ) { return ATP_Ox[k]; };

	std::vector<double> Get_ATP_NOx() { return ATP_NOx; };
	double Get_ATP_NOx( const unsigned long int k ) { return ATP_NOx[k]; };

	std::vector<double> Get_ATP2() { return ATP2; };
	double Get_ATP2( const unsigned long int k ) { return ATP2[k]; };

	std::vector<double> Get_ATP3() { return ATP3; };
	double Get_ATP3( const unsigned long int k ) { return ATP3[k]; };

	std::vector<double> Get_ConsATP() { return ConsATP; };
	double Get_ConsATP( const unsigned long int k ) { return ConsATP[k]; };

	std::vector<double> Get_ConsATP_1() { return ConsATP_1; };
	double Get_ConsATP_1( const unsigned long int k ) { return ConsATP_1[k]; };

	std::vector<double> Get_ConsATP_2() { return ConsATP_2; };
	double Get_ConsATP_2( const unsigned long int k ) { return ConsATP_2[k]; };

	std::vector<double> Get_ConsATP_3() { return ConsATP_3; };
	double Get_ConsATP_3( const unsigned long int k ) { return ConsATP_3[k]; };

	std::vector<double> Get_ConsATP_4() { return ConsATP_4; };
	double Get_ConsATP_4( const unsigned long int k ) { return ConsATP_4[k]; };

	std::vector<double> Get_ConsATP_5() { return ConsATP_5; };
	double Get_ConsATP_5( const unsigned long int k ) { return ConsATP_5[k]; };

	std::vector<double> Get_ATPtot() { return ATPtot; };
	double Get_ATPtot( const unsigned long int k ) { return ATPtot[k]; };

	std::vector<double> Get_ATPp() { return ATPp; };	
	double Get_ATPp( const unsigned long int k ) { return ATPp[k]; };	

	std::vector<double> Get_ATPmin() { return ATPmin; };	
	double Get_ATPmin( const unsigned long int k ) { return ATPmin[k]; };	


	std::vector<double> Get_ATPstart() { return ATPstart; };
	double Get_ATPstart( const unsigned long int k ) { return ATPstart[k]; };

	std::vector<double> Get_ATPprod() { return ATPprod; };
	double Get_ATPprod( const unsigned long int k ) { return ATPprod[k]; };

	std::vector<double> Get_ATPcons() { return ATPcons; };
	double Get_ATPcons( const unsigned long int k ) { return ATPcons[k]; };


	std::vector<double> Get_G_extra() { return G_extra; };
	double Get_G_extra( const unsigned long int k ) { return G_extra[k]; };

	std::vector<double> Get_A_extra() { return A_extra; };
	double Get_A_extra( const unsigned long int k ) { return A_extra[k]; };

	std::vector<double> Get_AcL_extra() { return AcL_extra; };
	double Get_AcL_extra( const unsigned long int k ) { return AcL_extra[k]; };

	
	std::vector<double> Get_pH() { return pH; };
	double Get_pH( const unsigned long int k ) { return pH[k]; };

	// double Get_H_extra() { return H_extra; };
	// double Get_CO2_extra() { return CO2_extra; };
	
	std::vector<double> Get_SensO2() { return SensO2; };
	double Get_SensO2( const unsigned long int k ) { return SensO2[k]; };

	std::vector<double> Get_ConsO() { return ConsO; };
	double Get_ConsO( const unsigned long int k ) { return ConsO[k]; };

	// double Get_ProdCO2() { return ProdCO2; };

	std::vector<double> Get_DNA_spread() { return DNA_spread; };
	double Get_DNA_spread( const unsigned long int k ) { return DNA_spread[k]; };

	
	std::vector<double> Get_M_T() { return M_T; };
	double Get_M_T( const int k ) { return M_T[k]; };


	std::vector<double> Get_DNA() { return DNA; };
	double Get_DNA( const unsigned long int k ) { return DNA[k]; };

	std::vector<double> Get_ConcS() { return ConcS; };
	double Get_ConcS( const unsigned long int k ) { return ConcS[k]; };


	std::vector<double> Get_protein() { return protein; };
	double Get_protein( const unsigned long int k ) { return protein[k]; };

	std::vector<double> Get_prot_rate() { return prot_rate; };
	double Get_prot_rate( const unsigned long int k ) { return prot_rate[k]; };

	std::vector<double> Get_DNA_rate() { return DNA_rate; };
	double Get_DNA_rate( const unsigned long int k ) { return DNA_rate[k]; };

	std::vector<double> Get_pRb() { return pRb; };
	double Get_pRb( const unsigned long int k ) { return pRb[k]; };

	std::vector<double> Get_cyclinD() { return cyclinD; };
	double Get_cyclinD( const unsigned long int k ) { return cyclinD[k]; };

	std::vector<double> Get_cyclinE() { return cyclinE; };
	double Get_cyclinE( const unsigned long int k ) { return cyclinE[k]; };

	std::vector<double> Get_cyclinX() { return cyclinX; };
	double Get_cyclinX( const unsigned long int k ) { return cyclinX[k]; };



	std::vector<double> Get_NpRbk() { return NpRbk; };
	double Get_NpRbk( const unsigned long int k ) { return NpRbk[k]; };

  vbl::Timing<double> myTiming;

	// overloaded =
	// Cells& operator=(const Cells& newcell);
	
	// setters piu' complessi

	// setter del raggio: attenzione non e' protetto e se r e' troppo piccolo si puo' avere un valore negativo di ATPp
	void Set_r( const unsigned long int k, const double newr) 
	{ 
		r[k] = newr; 
		surface[k] = 4.*PI*newr*newr; 
		volume[k] = surface[k]*newr/3.; 
		mass[k] = type[k]->density * volume[k];
		
  		volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);
		
		// variabili interne che dipendono dal volume (si assume comunque type->C1 > 0 )


		ATPp[k] = (volume[k] - type[k]->C2 * M[k] - type[k]->Vmin * (1.+DNA[k]))/type[k]->C1;						// inizializzazione di ATP pool


		
	};

	// setter del volume: attenzione non e' protetto e se il volume e' troppo piccolo si puo' avere un valore negativo di ATPp
	void Set_volume( const unsigned long int k, const double newvolume ) 
	{ 
		volume[k] = newvolume; 
		r[k] = pow(3.*newvolume/(4.*PI), (double)1./3.); 
		surface[k] = 4.*PI*r[k]*r[k]; 
		mass[k] = type[k]->density * newvolume;
		
  		volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);
		
		// variabili interne che dipendono dal volume
		if( phase[k] != dead )
			{


			ATPp[k] = (newvolume - type[k]->C2 * M[k] - type[k]->Vmin * (1.+DNA[k]))/type[k]->C1;						// inizializzazione di ATP pool

			}
	};

	// questa funzione ha senso solo se la cellula e' viva
	void Set_ATPp( const unsigned long int k, const double newATPp ) 
	{ 
		ATPp[k] = newATPp; 
		volume[k] = type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1;
		r[k] = pow(3.*volume[k]/(4.*PI), (double)1./3.); 
		surface[k] = 4.*PI*r[k]*r[k]; 	
		mass[k] = type[k]->density * volume[k];	
		
		volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);
	
	};
	

	// questa funzione ha senso solo se la cellula e' viva
	void Set_M( const unsigned long int k, const double newM ) 
	{ 
	
	M[k] = newM;
	
	ATPmin[k] = (type[k]->fATPmin)*(type[k]->C2 * M[k])/type[k]->C1;
	
	volume[k] = type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1;
	r[k] = pow(3.*volume[k]/(4.*PI), (double)1./3.); 
	surface[k] = 4.*PI*r[k]*r[k]; 	
	mass[k] = type[k]->density * volume[k];	
	
	volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);
	
	};

	// questa funzione ha senso solo se la cellula e' viva e serve a settare contemporaneamente M e ATPp (e le variabili derivate)
	void Set_M_and_ATPp( const unsigned long int k, const double newM,  const double newATPp ) 
	{ 
	
	M[k] = newM;
	ATPp[k] = newATPp; 
	
	ATPmin[k] = (type[k]->fATPmin)*(type[k]->C2 * M[k])/type[k]->C1;
	
	volume[k] = type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1;
	r[k] = pow(3.*volume[k]/(4.*PI), (double)1./3.); 
	surface[k] = 4.*PI*r[k]*r[k]; 	
	mass[k] = type[k]->density * volume[k];	
	
	volume_extra[k] = surface[k]*(type[k]->extvolume_thickness)*(type[k]->extvolume_fraction);
	
	};
	
	// Function that controls the consistency of linked mitochondria, volume, ATPp values
	int CheckMVA( const unsigned long int k )
	{
	
	int return_code = 0;
	static int count = 0;
    
    if(phase[k] != dead)
    {
      if(M[k] < -(std::numeric_limits<double>::epsilon( )))
      {
        std::cout << "Inconsistent value of M in the cell " << name[k] << ": M=" << M[k] << std::endl;
        return_code = -1;
      }
      if(volume[k] < type[k]->Vmin)
      {
        std::cout << "Inconsistent volume value in the cell " << " k: " << k << name[k] 
        << ": Vmin=" << type[k]->Vmin << "volume: " << volume[k] << std::endl;
        return_code = -2;
      }
        if(ATPp[k] < -(std::numeric_limits<double>::epsilon( )))
            {
            std::cout << "Inconsistent value of ATPp in the cell " << name[k] << ": ATPp=" << ATPp[k] << std::endl;
            return_code = -3;
            }
        if( fabs( ATPmin[k] - ((type[k]->fATPmin)*((type[k]->C2) * M[k])/(type[k]->C1)) ) > (std::numeric_limits<double>::epsilon( )) )
            {
              /* T.F.
               * this often broke when running on the cluster
               */
            std::cout << "Inconsistent value of ATPmin in the cell " << name[k] << ": ATPmin=" << ATPmin[k] << std::endl;
            std::cout << std::scientific << "ATPmin[k] = " << ATPmin[k] << 
                                  "type[k]->fATPminVmin " << type[k]->fATPmin <<
                                  "type[k]->C2 = " <<     type[k]->C2<<
                                  "type[k]->C1 = " <<     type[k]->C1<<
                                  "M[k] = "        <<     M[k] << std::endl;
            return_code = -4;
            }
        if( fabs( volume[k] - (type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1) )/volume[k] > 1.e-6 )
            {
            std::cout << "Inconsistent values ​​of ATPp, M, volume, in the cell " << name[k] << std::endl;
            std::cout << std::scientific << "volume = " << volume[k] << " = type->Vmin * (1.+DNA) + type->C2 * M + ATPp * type->C1 = " << (type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1) << std::endl;
            std::cout << "DNA: " << DNA[k] << std::endl;
            std::cout << "M: " << M[k] << std::endl;
            std::cout << "ATPp: " << ATPp[k] << std::endl;
            return_code = -5;
            }
        }
		
	if(return_code < 0) 
		{
		std::cout << std::endl;
		count++;
		}

	
	// dopo 20 errori di questo tipo il programma si ferma
	if(count > 20) exit(-1);
	
	
	
	return return_code;

	}
	

// *** fine dei metodi per la gestione della parte biofisica *** 

  boost::property_tree::ptree as_ptree() const;
  void assign(const boost::property_tree::ptree &pt);
};


  //CellsSystem CellsSystem;
}//namespace vbl{

#endif //#ifndef CELLSSYSTEM_H

