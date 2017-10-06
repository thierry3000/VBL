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

#include "BloodVessel.h"
// to use new and old version of TimersAdvance and TimersAdvanceUntil simultaniously
#include <boost/optional.hpp>

//for timing output
#include <iostream>
#include <fstream>
#define W_timing
// #ifndef NDEBUG
//   // means debugging
//   #define W_timing
// #endif
struct ReadInParameters
{
  /* imulation with dispersed cells (0) or full 3D (1) */
  int sim_type;
};
class CellsSystem
{

private:

// output file
string output_filename;
ofstream output_file;

// log file
string log_filename;
ofstream log_file;

// screen dump file
string screen_dump_filename;
ofstream screen_dump_file;

// configuration file
string configuration_filename;
ofstream configuration_file;

// configuration file (binary)
string configuration_b_filename;
ofstream configuration_b_file;

// error log file
string errorlog_filename;
ofstream errorlog_file;

// convergence log file
string convlog_filename;
ofstream convlog_file;

// cell log file
string cell_filename;
ofstream cell_file;

// environment log file
string env_filename;
ofstream env_file;

// flow file (binary only)
string flow_b_filename;
ofstream flow_b_file;

// CellsSystem files
string cellsys_in_filename;
ifstream cellsys_in_file;

string cellsys_out_filename;
ofstream cellsys_out_file;

// nome del file CellType
string CellTypeFile;

// nome del file CellTypeAlt
string CellTypeFileAlt;

// nome del file Environment
string EnvironmentFile;

// nome del file di comando
string commandFile;


// *** dati per la gestione del sistema di cellule ***

// questo serve a identificare il tipo di simulazione
int sim_type;

// questo identifica il tipo di distribuzione iniziale delle cellule (nel caso di simulazione Full3D)
int initial_cell_dist;

// timers
#ifdef _OPENMP
double CPU_time, CPU_time_0, CPU_time_1, CPU_time_2;	// CPU time
#else
clock_t CPU_time, CPU_time_0, CPU_time_1, CPU_time_2;	// CPU time
#endif
double t_CPU, t_CPU_int, t_CPU_0, delta_t_CPU;		// CPU time (real)
double t_CPU_max;			// tempo di CPU massimo per una singola frazione di run (in s)
time_t	time_now, time_old;		// tempo reale trascorso dall'ultimo step
double time_from_CGAL;		// tempo simulato trascorso dall'ultimo update della triangolazione di Delaunay
double timing;

double dt;					// timestep
double dt_sm;               // timestep slow motion
double t;					// tempo attuale di simulazione
double t_ini;				// tempo di inizializzazione
double treal;				// tempo di simulazione vero (dopo l'inizializzazione)
double tmax;				// tempo max di simulazione
double tsm_start;           // tempo di start dello slow motion
double tsm_stop;            // tempo di stop dello slow motion
bool slow_motion;           // indica lo stato dello slow motion


unsigned long nstep;			// numero del passo
unsigned long nstep_start;		// numero del passo alla partenza della simulazione (quando ready2start diventa true)
unsigned long nmax;				// nmax = floor(tmax/dt) numero max di passi
int idum;						// variabile globale utilizzata dalle routine di numeri casuali di Numerical Recipes 2
bool ready2start;				// Global variable by which the simulation can be performed on a regular basis
bool faketumAtCurrentTime;
unsigned long nprint;			// intervallo (in numero di passi) tra i passi di stampa su file
unsigned long nscreen;			// intervallo (in numero di passi) tra i passi di stampa su schermo
unsigned long nconfiguration;	// numero delle configurazioni scritte su file

// precisione
double eps;				// precisione del passo di diffusione
double delta_vmax;			// precisione della determinazione di velocita'

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
string machine;					// nome della macchina
int run;						// numero del run
string dir;						// nome del directory di output
int part;						// per un run suddiviso in parti eseguite in tempi diversi, questa variabile indica la parte attuale

// numero iniziale di cellule
unsigned long nstart;

// numero di cellule 
unsigned long ncells;

// numero di cellule vive alla fine del passo di metabolismo
unsigned long alive;

// numero di tipi cellulari
unsigned long ntypes;		

// ultimo nome assegnato (serve a gestire in modo univoco l'assegnazione del nome alle cellule)
unsigned long lastname;

// ambiente iniziale
Environment Env_0;

// ambiente attuale
Environment Env;

// variazione dello stato dell'ambiente 
//Environment delta_Env;

// flusso (flag che indica se il flusso e' non nullo, e classe che rappresenta il segnale non nullo)
bool flowON;
EnvironmentalSignal flowSignal;
    
// flusso dell'ossigeno (flag che indica il tipo del funzionamento nel modo bioreattore; se la flag è accesa allora viene modulato anche
// l'ossigeno, che in caso contrario resta invece fisso al valore ambientale); questa flag non serve a nulla se flowON e' spenta.
bool oxygenflowON;

// dose rate (flag che indica se la dose e' non nulla, e classe che rappresenta il segnale non nullo)
bool doseON;
EnvironmentalSignal dose_rateSignal;

// definizione del vettore dei fenotipi cellulari
vector<CellType> CellTypeVector;

// parametri legati alla geometria e alla dinamica del cluster
double maxdr;				// spostamento max nel sistema di cellule in un passo di dinamica

// parametri metabolici del cluster
double O2Flow;				// flusso di O2 dalla periferia alle regioni interne (in kg/s)

double AcLFlow;			// flusso di AcL nell'ambiente (in kg/s)

// *** fine dei dati per la gestione del sistema di cellule ***

//****************************************************************************************************

// *** vettori di dati associati alle singole cellule ***


	// informazioni fondamentali sulla cellula: nome ed eventuale marcatura (che viene ereditata), 
	// fenotipo e caratteristiche generali dell'ambiente
		
	vector<unsigned long> name;			// identificatore della cellula 	
	vector<int> mark;					// label applicato alla cellula
	vector<CellType*> type;				// puntatore a un CellType (fenotipo cellulare)
	
	// temperatura della cellula
	
	vector<double> Temperature;	
		
	// stato cellulare (associato al metabolismo)
	
	vector<CellPhase> phase;			// fase cellulare
	vector<int> death_condition;		// variabile che registra la ragione della morte cellulare
	vector<float> age;					// eta' della cellula (tempo in secondi dalla nascita)
	vector<float> phase_age;			// eta' dello stato cellulare attuale (tempo in secondi dall'inizio)
	vector<float> age_mother;			// eta' della cellula madre (tempo in secondi dalla nascita alla mitosi)
	vector<int> n_mitosis;				// numero di mitosi dall'inizio della simulazione
	
	// Geometric and topological variables

	vector<double> x;					// posizione nello spazio del centro della cellula
	vector<double> y;
	vector<double> z;
	vector<double> vx;					// velocità della cellula
	vector<double> vy;
	vector<double> vz;
	vector<double> r;					// raggio cellulare	
	vector<double> surface;             // superficie cellulare
	vector<double> volume;				// volume cellulare
	vector<double> mass;				// massa cellulare
	
	// copie per utilizzo nel metodo di calcolo della geometria
	vector<double> vxnew;				// velocità della cellula
	vector<double> vynew;
	vector<double> vznew;

	vector<double> fx;					// forze
	vector<double> fy;
	vector<double> fz;
	
	vector<Point> v;					// vettore dei punti passati alla triangolazione di CGAL


	vector<double> volume_extra;		// volume della regione extracellulare che circonda la cellula
	
	vector<int> neigh;					// numero di vicini
	vector< vector<int> > vneigh;		// vettore dei vicini
	vector< vector<double> > vdist;     // vettore delle distanze dai vicini
	vector< vector<double> > vcsurf;	// vettore delle superfici di contatto con i vicini (calcolo approx)
	vector< vector<double> > gnk;		// vettore dei fattori geometrici
	vector<double> contact_surf;		// area totale della superficie di contatto con le cellule adiacenti
	
	vector<bool> isonCH;				// label che indica se la cellula si trova sul convex hull
	vector<bool> isonAS;				// label che indica se la cellula si trova sull'alpha shape
	vector<int> isonBV;					// this label is true if the cell is in contact
										// with a blood vessel
										// if a cell is not in contact with blood vessel, 
										// then isonBV[n] = 0, else isonBV[n] = position 
										// of vessel in BloodVesselVector + 1
										// this shift is done so that isonBV can also be 
										// used as a boolean variable 
										
	vector<double> env_surf;			// superficie di contatto con l'ambiente
	vector<double> g_env;				// fattore geometrico relativo al contatto con l'ambiente

	vector<double> bv_surf;				// superficie di contatto con vasi sanguigni
	vector<double> g_bv;				// fattore geometrico relativo al contatto con vasi sanguigni

	vector<double> M;					// numero di mitocondri 

	// lista delle variabili metaboliche all'interno della cellula
	
	vector<double> G;					// quantita' (massa) interna di glucosio
	vector<double> G6P;                 // quantita' (massa) interna di G6P
	vector<double> O2;					// quantita' (massa) interna di ossigeno
	vector<double> store;				// riserva in unita' di massa di glucosio equivalente
	vector<double> A;					// quantita' (massa) interna di altri nutrienti (glutammina)
	vector<double> AcL;                 // quantita' (massa) interna di acido lattico
    
    vector<double> h;                   // parametrizzazione dell'attivita' di trasporto del glucosio in funzione di O2
	vector<double> pHi;                 // valore del pH interno alla cellula
	// double H;                        // quantita' (massa) interna di protoni liberi
	// double CO2;                      // quantita' (massa) interna di anidride carbonica

	vector<double> protein;             // quantita' totale di proteine nelle cellule (in kg)
	vector<double> prot_rate;			// rate di produzione delle proteine
	vector<double> DNA;                 // quantita' di DNA prodotto per la replicazione 
                                        // (in molecole di DNA. Una molecola = intero genoma)
	vector<double> DNA_rate;			// rate di produzione del DNA

	// rates associate al glucosio, allo store e all'acido lattico
	
	vector<double> GAbsRate;			// rate di assorbimento del glucosio
	vector<double> GConsRate;			// rate di consumo del glucosio
	vector<double> AAbsRate;			// rate di assorbimento della glutammina
	vector<double> AConsRate;			// rate di consumo del glutammmina
	vector<double> StoreFillRate;		// rate di riempiemento dello store
	vector<double> StoreConsRate;		// rate di consumo dello store
	vector<double> AcLRate;				// rate di produzione dell'acido lattico
	vector<double> AcLOutRate;			// rate di trasporto all'esterno dell'acido lattico

	// calcolo dell'ATP (richiede parecchie variabili ausiliarie che hanno significato biologico)
	
	vector<double> ATP_St;				// ATP standard
	vector<double> ATP_Ox;				// rate di produzione di ATP tramite fosforilazione ossidativa
	vector<double> ATP_NOx;             // rate di produzione di ATP tramite glicolisi anaerobica
	vector<double> ATP2;				// rate di produzione di ATP dallo store
	vector<double> ATP3;				// rate di produzione di ATP dagli altri nutrienti
	vector<double> ConsATP;             // rate di consumo di ATP associato al metabolismo
	vector<double> ConsATP_1;			// rate di consumo di ATP associato ai mitocondri
	vector<double> ConsATP_2;			// rate di consumo di ATP associato alla produzione di proteine
	vector<double> ConsATP_3;			// rate di consumo di ATP associato alla produzione di DNA
	vector<double> ConsATP_4;			// rate di consumo di ATP associato alla pompa protonica
	vector<double> ConsATP_5;			// rate di consumo di ATP associato alla proliferazione dei mitocondri
	vector<double> ATPtot;				// rate totale di variazione dell'ATP = somma dei rates
	vector<double> ATPp;				// massa totale di ATP all'interno della cellula
	vector<double> ATPmin;				// valore minimo dell'ATPp (al di sotto di questo la cellula muore)

	// altre variabili correlate all'ATP
	
	vector<double> ATPstart;			// ATP iniziale nella cellula (in kg)
	vector<double> ATPprod;             // ATP prodotto dalla cellula (in kg) nel passo metabolico
	vector<double> ATPcons;             // ATP consumato dalla cellula (in kg) nel passo metabolico


	// variabili metaboliche negli spazi extracellulari
	
	vector<double> G_extra;             // massa del glucosio nella matrice extracellulare (negli spazi intercellulari)
	vector<double> A_extra;             // massa degli altri nutrienti (glutammina) nella matrice extracellulare (negli spazi intercellulari)
	vector<double> AcL_extra;			// massa dell'AcL nella matrice extracellulare (negli spazi intercellulari)
	
	vector<double> pH;					// valore del pH extracellulare
	// double H_extra;					// massa dei protoni liberi nella matrice extracellulare (negli spazi intercellulari)
	// double CO2_extra;				// massa dell'anidride carbonica nella matrice extracellulare (negli spazi intercellulari)

	// altre quantita' derivate
	
	vector<double> SensO2;				// frazione di ossigeno disponibile rispetto la richiesta
	vector<double> ConsO;				// rate di consumo dell'ossigeno (kg/s)
	// double ProdCO2;					// rate di produzione di CO2 (kg/s)
	
	// proteine e DNA

	vector<double> DNA_spread;			// variazione individuale della velocita' di sintesi del DNA 
                                        // (parte relativa al consumo di ATP)
                                        // questa variazione e' una modellizzazione grezza della velocita' 
                                        // variabile di sintesi dovuta ai 
                                        // danni al DNA (da modificare in futuro)
	
	vector<double> M_T;                 // durata media fase M (in s)
	vector<double> pRb;                 // quantita' (massa) della pRb

	vector<double> ConcS;				// concentrazione MOLARE della sostanza S che fa da substrato alla reazione di MM per la soglia

	vector<double> cyclinD;             // quantita' totale di ciclina D (in kg)
	vector<double> cyclinE;             // quantita' totale di ciclina E (in kg)
	vector<double> cyclinX;             // quantita' totale di cicline (A e B) attive in fase G2 (in kg)
	
	vector<double> NpRbk;				// numero di molecole di pRb attive

// copie per i calcoli in Diff

	vector<double> volumeOld;			// vettore dei volumi cellulari (valore vecchio)
	vector<double> volumeNew;			// vettore dei volumi cellulari (valore nuovo)
	vector<double> volume_extraOld;     // volume dei volume extracellulari (valore vecchio)
	vector<double> volume_extraNew;     // volume dei volume extracellulari (valore nuovo)
                                        // anche se il volume non e' una variabile dinamica, va immagazzinato in vettori 
                                        // per poter effettuare le somme necessarie al calcolo della diffusione
																
	vector<double> MitOld;				// vettore dei mitocondri (valore vecchio)
	vector<double> MitNew;				// vettore dei mitocondri (valore nuovo)

	vector<double> pHiOld;				// vettore dei pH cellulari (valore vecchio)
	vector<double> pHiNew;				// vettore dei pH cellulari (valore nuovo)
	vector<double> pHOld;				// vettore dei pH extracellulari (valore vecchio)
	vector<double> pHNew;				// vettore dei pH extracellulari (valore nuovo)
                                        // anche se il pH non e' una variabile dinamica, va immagazzinato in vettori 
                                        // per poter effettuare le somme necessarie al calcolo della diffusione
	
	vector<double> mGinOld;             // vettore della massa di glucosio dentro le cellule (valore vecchio)
	vector<double> mGinNew;             // vettore della massa di glucosio dentro le cellule (valore nuovo)
	vector<double> mGextOld;			// vettore della massa di glucosio negli spazi extracellulari (valore vecchio)
	vector<double> mGextNew;			// vettore della massa di glucosio negli spazi extracellulari (valore nuovo)

	vector<double> mG6POld;             // vettore della massa di G6P nelle cellule (valore vecchio)
	vector<double> mG6PNew;             // vettore della massa di G6P nelle cellule (valore nuovo)

	vector<double> mO2Old;				// vettore della massa di O2 nelle cellule (valore vecchio)
	vector<double> mO2New;				// vettore della massa di O2 nelle cellule (valore nuovo)

	vector<double> StoreOld;			// vettore dello Store nelle cellule (valore vecchio)
	vector<double> StoreNew;			// vettore dello Store nelle cellule (valore nuovo)

	vector<double> mAinOld;             // vettore della massa di glutammina dentro le cellule (valore vecchio)
	vector<double> mAinNew;             // vettore della massa di glutammina dentro le cellule (valore nuovo)
	vector<double> mAextOld;			// vettore della massa di glutammina negli spazi extracellulari (valore vecchio)
	vector<double> mAextNew;			// vettore della massa di glutammina negli spazi extracellulari (valore nuovo)

	vector<double> mAcLinOld;			// vettore della massa di acido lattico dentro le cellule (valore vecchio)
	vector<double> mAcLinNew;			// vettore della massa di acido lattico dentro le cellule (valore nuovo)
	vector<double> mAcLextOld;			// vettore della massa di acido lattico negli spazi extracellulari (valore vecchio)
	vector<double> mAcLextNew;			// vettore della massa di acido lattico negli spazi extracellulari (valore nuovo)

	vector<double> ATPpOld;             // vettore dell'ATPp nelle cellule (valore vecchio)
	vector<double> ATPpNew;             // vettore dell'ATPp nelle cellule (valore nuovo)

	// altre allocazioni di memoria utili per il metodo di diffusione
	
	vector<double> proteinNew;			// vettore della massa di proteina nelle cellule
	vector<double> pRbNew;				// vettore della massa di pRb nelle cellule
	vector<double> delta_protein;		// vettore della variazione di massa totale di proteine nelle cellule
	vector<double> ConcSNew;			// concentrazione del substrato S per il calcolo delle soglie dei checkpoints

	vector<double> DNANew;				// vettore che contiene la lunghezza relativa della molecola di DNA (1 = intera molecola)

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
    std::vector<vbl::BloodVessel> BloodVesselVector;  // vettore dei vasi sanguigni
    
// *** fine dei dati associati a vasi sanguigni


//****************************************************************************************************


public:

//****************************************************************************************************

// *** Methods for managing the system ***

// Redeployment of the dynamic reserve
void Set_reserve(const int reserve); // cell vectors
void Set_BV_reserve(const int reserve); // blood vessel vector

// builder that builds a cell array of no length, but assigns a dynamic reserve to carriers
// Note that the creation of the CellsSystem also resets cell counts, types, blood vessels
CellsSystem(const int reserve, const int reserve_bv) { ncells = 0; ntypes = 0; nbv=0; Set_reserve(reserve); Set_BV_reserve(reserve_bv); };
// default builder, builds a set of cells of zero length and assigns the standard dynamic reserve to the carriers
CellsSystem() { ncells = 0; ntypes = 0; nbv=0; Set_reserve(RESERVE); Set_BV_reserve(RESERVE_BV); };

// aggiunta di cellule non inizializzate al sistema
void AddCells( const int newcells );
// Adding a single standardized standardized cell
void AddInitializedCell(int& idum, CellType* cType, Environment* cEnv);
// inizializzazione standard dell'intero sistema di cellule
void InitCells( int& idum, CellType* cType, Environment* cEnv ) { for(unsigned long int k=0; k<ncells; k++) AddInitializedCell(idum, cType, cEnv ); };
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
void PrintCell( ostream& stream, const unsigned long int k );
// Printing data in the form of a single string readable by a spreadsheet program
void PrintCellData( const unsigned long int k, ofstream& stream, long int nrec );

// rimozione di una cellula (il metodo gestisce tutti i vettori)
void RemoveCell( const unsigned long  n );

// altri metodi

// calcolo della forza tra le cellule nella configurazione attuale
void GetForces();
// calcolo della posizione e della velocità di ciascuna cellula (si calcola solo se ci sono almeno due cellule)
void NewPositionsAndVelocities( );


// getters
string Get_Commands( ) { return commandFile; };
string Get_CellTypeFile( ) { return CellTypeFile; };
string Get_CellTypeFileAlt( ) { return CellTypeFileAlt; };
string Get_EnvironmentFile( ) { return EnvironmentFile; };
int Get_sim_type() { return sim_type; };
int Get_initial_cell_dist() { return initial_cell_dist; };
double Get_dt() { return dt; };
double Get_dt_sm() { return dt_sm; };
bool Get_slow_motion() { return slow_motion; };
double Get_t() { return t; };
double Get_treal() { return treal; };
double Get_tmax() { return tmax; };
double Get_time_from_CGAL() { return time_from_CGAL; };

unsigned long Get_ncalls() { return ncalls; };
unsigned long Get_nstep() { return nstep; };
unsigned long Get_nstep_start() { return nstep_start; };
unsigned long Get_nmax() { return nmax; };

int Get_idum() { return idum; };

bool Get_ready2start() { return ready2start; };
unsigned long Get_nprint() { return nprint; };
unsigned long Get_nscreen() { return nscreen; };
unsigned long Get_nconfiguration() { return nconfiguration; };

double Get_eps() { return eps; };
double Get_delta_vmax() { return delta_vmax; };
double Get_nrepeats() { return nrepeats; };
long int Get_min_nrepeats() { return min_nrepeats; };
double Get_nrepeats_average() { return nrepeats_average; };
double Get_nrepeats_max() { return nrepeats_max; };
double Get_nrepeats_min() { return nrepeats_min; };


unsigned long Get_nstart() { return nstart; };
unsigned long Get_ncells() { return ncells; };
unsigned long Get_alive() { return alive; };
unsigned long Get_ntypes() { return ntypes; };
unsigned long Get_lastname() { return lastname; };
bool Get_flowON() { return flowON; };
bool Get_doseON() { return doseON; };
Environment Get_Env() { return Env; };
Environment Get_Env_0() { return Env_0; };
EnvironmentalSignal Get_flowSignal() { return flowSignal; };
EnvironmentalSignal Get_dose_rateSignal() { return dose_rateSignal; };
double Get_maxdr() { return maxdr; };
double Get_O2Flow() { return O2Flow; };
double Get_AcLFlow() { return AcLFlow; };

// setters
void Set_Commands( const string newcommandFile ) { commandFile = newcommandFile; };
void Set_CellTypeFile( string newCellTypeFile ) { CellTypeFile = newCellTypeFile; };
void Set_CellTypeFileAlt( string newCellTypeFile ) { CellTypeFileAlt = newCellTypeFile; };
void Set_EnvironmentFile( string newEnvironmentFile ) { EnvironmentFile = newEnvironmentFile; };
void Set_idum ( int newidum ) { idum = newidum; };
void Set_dt( double newdt ) { dt = newdt; };
void Set_time_from_CGAL( double newtime_from_CGAL ) { time_from_CGAL = newtime_from_CGAL; };
void Set_ready2start( bool newr2s ) { ready2start = newr2s; };
void Set_nconfiguration( unsigned long newnconfiguration ) { nconfiguration = newnconfiguration; };
void Step_nconfiguration() { nconfiguration++; };
void Set_eps( double neweps ) { eps = neweps; };
void Set_delta_vmax( double newdelta_vmax ) { delta_vmax = newdelta_vmax; };

// inizializzazione del run da terminale
void InitializeCellsSystem( bool terminal );

// inizializzazione del run da file
void InitializeCellsSystem( const string filename );

// metodo per gli eventi cellulari (crescita, mitosi, etc.; restituisce una flag che indica se c'e' stata almeno una mitosi nel sistema di cellule)
bool CellEvents( );

// metodo per la pulizia della memoria (al momento usa operazioni non molto efficienti)
void CleanCellsSystem( );

// geometria
void Geometry();
void NoGeometry();	// versione ridotta, calcoli minimi per cellule disperse

// meccanica del cluster
void Dynamics( );
void DummyDynamics( ); // dinamica dummy

// summary printout
void Printout();

// close files
void CloseOutputFiles() { cout << "\nFine del run " << run << " al passo " << nstep << endl; 
						output_file.close(); 
						log_file.close(); 
						screen_dump_file.close();
						errorlog_file.close();
						cell_file.close();
						env_file.close();
						convlog_file.close(); };

// summary printout on file
void Print2file();

// singola stampa su logfile
void Print2logfile(string str);			// cellula ben formattata
void PrintAll2logfile(string str);		// stampa TUTTE le cellule su logfile (!!!)

// definizione del run
void RunDefinition( );
void RunDefinition( string run_name );

// scrittura e lettura del CellsSystem
void WriteCellsSystem( );
void ReadCellsSystem( );


// avanzamento dei timers
//bool TimersAdvance( );
// if no endtime is present, the behaviour of old TimersAdvance is restored
bool TimersAdvanceUntil( boost::optional<double> endtime );

// CPU timing
void CPU_timer( timer_button button );

// Timing
double Timing( bool reset );

// distanza tra cellule con indici j e k
double Distance(const int j, const int k);

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

unsigned int runMainLoop( boost::optional<double> endtime);
//unsigned int runMainLoop( );

// *** metodi per la gestione della parte biofisica *** 

	

	// setters (nella forma di vettore e di singolo elemento)

	void Set_name(const vector<unsigned long>& newname) { name = newname; };
	void Set_name(const int k, const unsigned long newname) { name[k] = newname; };
	
	void Set_mark(const vector<int>& newmark) { mark = newmark; };
	void Set_mark(const int k, const int newmark) { mark[k] = newmark; };
	
	void Set_type(const vector<CellType*>& cType) { type = cType; };
	void Set_type(const int k, CellType* cType) { type[k] = cType; };
	
	void Set_Temperature(const vector<double>& newTemperature) { Temperature = newTemperature; };
	void Set_Temperature(const int k, const double newTemperature) { Temperature[k] = newTemperature; };
	
	void Set_phase(const vector<CellPhase>& newphase) { phase = newphase; };
	void Set_phase(const int k, const CellPhase newphase) { phase[k] = newphase; };
	
	void Set_death_condition( const vector<int>& newdeath_condition ) { death_condition = newdeath_condition; };
	void Set_death_condition( const int k, const int newdeath_condition ) { death_condition[k] = newdeath_condition; };

	void Set_age(const vector<float>& newage) { age = newage; };
	void Set_age(const int k, const float newage) { age[k] = newage; };

	void Set_phase_age(const vector<float>& newphase_age) { phase_age = newphase_age; };
	void Set_phase_age(const int k, const float newphase_age) { phase_age[k] = newphase_age; };

	void Set_age_mother(const vector<float>& newage_mother) { age_mother = newage_mother; };
	void Set_age_mother(const int k, const float newage_mother) { age_mother[k] = newage_mother; };

	void Set_n_mitosis(const vector<int>& newn_mitosis ) { n_mitosis = newn_mitosis; };
	void Set_n_mitosis(const int k, const int newn_mitosis ) { n_mitosis[k] = newn_mitosis; };

	void Add_BloodVessel_at(uint index, vbl::BloodVessel NewBV ); 
    
    // geometria, cinematica e dinamica

	void Set_x( const vector<double>& newx ) { x = newx; };
	void Set_x( const int k, const double newx ) { x[k] = newx; };

	void Set_y( const vector<double>& newy ) { y = newy; };
	void Set_y( const int k, const double newy ) { y[k] = newy; };

	void Set_z( const vector<double>& newz ) { z = newz; };
	void Set_z( const int k, const double newz ) { z[k] = newz; };

	void Set_vx( const vector<double>& newvx ) { vx = newvx; };
	void Set_vx( const int k, const double newvx ) { vx[k] = newvx; };

	void Set_vy( const vector<double>& newvy ) { vy = newvy; };
	void Set_vy( const int k, const double newvy ) { vy[k] = newvy; };

	void Set_vz( const vector<double>& newvz ) { vz = newvz; };
	void Set_vz( const int k, const double newvz ) { vz[k] = newvz; };
	
			
	void Init_fx( ) { fx.clear(); };
	void Init_fy( ) { fy.clear(); };
	void Init_fz( ) { fz.clear(); };
	
	
	// fEnd of useful setters to the geometric part

	void Set_volume_extra( const vector<double>& newvolume_extra ) { volume_extra = newvolume_extra; };
	void Set_volume_extra( const int k, const double newvolume_extra ) { volume_extra[k] = newvolume_extra; };


	void Set_neigh( const vector<int>& neighin ) { neigh = neighin; };
	void Set_neigh( const int k, const int neighin ) { neigh[k] = neighin; };

	// questi setters esistono solo nella forma per singole cellule (inseriscono vettori di lunghezza variabile)
	void Set_vneigh( const int k, int* vneighin ) { vneigh[k].clear(); vneigh[k].insert( vneigh[k].begin(), vneighin, vneighin+neigh[k]); };
	void Set_vdist( const int k, double* vdistin ) { vdist[k].clear(); vdist[k].insert( vdist[k].begin(), vdistin, vdistin+neigh[k]); };
	void Set_vcsurf( const int k, double* vcsurfin ) { vcsurf[k].clear(); vcsurf[k].insert( vcsurf[k].begin(), vcsurfin, vcsurfin+neigh[k]); };
	void Set_gnk( const int k, double* newgnk ) { gnk[k].clear(); gnk[k].insert( gnk[k].begin(), newgnk, newgnk+neigh[k]); };
	// fine della parte dei setters non standard

	
	
	void Set_contact_surf( const vector<double>& newcontact_surf ) { contact_surf = newcontact_surf; };
	void Set_contact_surf( const int k, const double newcontact_surf ) { contact_surf[k] = newcontact_surf; };

	void Set_isonCH( const vector<bool>& isonCHnow ) { isonCH = isonCHnow; };
	void Set_isonCH( const int k, const bool isonCHnow ) { isonCH[k] = isonCHnow; };

	void Set_isonAS( const vector<bool>& isonASnow ) { isonAS = isonASnow; };
	void Set_isonAS( const int k, const bool isonASnow ) { isonAS[k] = isonASnow; };

	void Set_isonBV( const vector<int>& isonBVnow ) { isonBV = isonBVnow; };
	void Set_isonBV( const int k, const int isonBVnow ) { isonBV[k] = isonBVnow; };

	void Set_env_surf( const vector<double>& newenv_surf ) { env_surf = newenv_surf; };
	void Set_env_surf( const int k, const double newenv_surf ) { env_surf[k] = newenv_surf; };

	void Set_g_env( const vector<double>& newg_env ) { g_env = newg_env; };
	void Set_g_env( const int k, const double newg_env ) { g_env[k] = newg_env; };

	void Set_bv_surf( const vector<double>& newbv_surf ) { bv_surf = newbv_surf; };
	void Set_bv_surf( const int k, const double newbv_surf ) { bv_surf[k] = newbv_surf; };

	void Set_g_bv( const vector<double>& newg_bv ) { g_bv = newg_bv; };
	void Set_g_bv( const int k, const double newg_bv ) { g_bv[k] = newg_bv; };


	void Set_G( const vector<double>& newG ) { G = newG; };
	void Set_G( const int k, const double newG ) { G[k] = newG; };

	void Set_G6P( const vector<double>& newG6P ) { G6P = newG6P; };
	void Set_G6P( const int k, const double newG6P ) { G6P[k] = newG6P; };

	void Set_O2( const vector<double>& newO2 ) { O2 = newO2; };
	void Set_O2( const int k, const double newO2 ) { O2[k] = newO2; };

	void Set_store( const vector<double>& newstore ) { store = newstore; };
	void Set_store( const int k, const double newstore ) { store[k] = newstore; };

	void Set_A( const vector<double>& newA ) { A = newA; };
	void Set_A( const int k, const double newA ) { A[k] = newA; };

	void Set_AcL( const vector<double>& newAcL ) { AcL = newAcL; };
	void Set_AcL( const int k, const double newAcL ) { AcL[k] = newAcL; };
	
	void Set_pHi( const vector<double>& newpHi ) { pHi = newpHi; };
	void Set_pHi( const int k, const double newpHi ) { pHi[k] = newpHi; };

	// void Set_H( const double newH ) { H = newH; };
	// void Set_CO2( const double newCO2 ) { CO2 = newCO2; };
	
	
	void Set_G_extra( const vector<double>& newG_extra ) { G_extra = newG_extra; };
	void Set_G_extra( const int k, const double newG_extra ) { G_extra[k] = newG_extra; };

	void Set_A_extra( const vector<double>& newA_extra ) { A_extra = newA_extra; };
	void Set_A_extra( const int k, const double newA_extra ) { A_extra[k] = newA_extra; };

	void Set_AcL_extra( const vector<double>& newAcL_extra ) { AcL_extra = newAcL_extra; };
	void Set_AcL_extra( const int k, const double newAcL_extra ) { AcL_extra[k] = newAcL_extra; };


	void Set_pH( const vector<double>& newpH ) { pH = newpH; };
	void Set_pH( const int k, const double newpH ) { pH[k] = newpH; };
	
	
	void Set_protein( const vector<double>& newprot ) { protein = newprot; }; 
	void Set_protein( const int k, const double newprot ) { protein[k] = newprot; }; 

	void Set_prot_rate( const vector<double>& newprot_rate ) { prot_rate = newprot_rate; };
	void Set_prot_rate( const int k, const double newprot_rate ) { prot_rate[k] = newprot_rate; };

	void Set_DNA_rate( const vector<double>& newDNA_rate ) { DNA_rate = newDNA_rate; }; 
	void Set_DNA_rate( const int k, const double newDNA_rate ) { DNA_rate[k] = newDNA_rate; }; 

	void Set_pRb ( const vector<double>& newpRb ) { pRb = newpRb; };
	void Set_pRb ( const int k, const double newpRb ) { pRb[k] = newpRb; };

	void Set_cyclinD( const vector<double>& newcyclinD ) { cyclinD = newcyclinD; };
	void Set_cyclinD( const int k, const double newcyclinD ) { cyclinD[k] = newcyclinD; };

	void Set_cyclinE( const vector<double>& newcyclinE ) { cyclinE = newcyclinE; };
	void Set_cyclinE( const int k, const double newcyclinE ) { cyclinE[k] = newcyclinE; };

	void Set_cyclinX( const vector<double>& newcyclinX ) { cyclinX = newcyclinX; };
	void Set_cyclinX( const int k, const double newcyclinX ) { cyclinX[k] = newcyclinX; };

	void Set_concS( const vector<double>& newconcS ) { ConcS = newconcS; };
	void Set_concS( const int k, const double newconcS ) { ConcS[k] = newconcS; };

	
	void Set_DNA( const vector<double>& newDNA ) { DNA = newDNA; };
	void Set_DNA( const int k, const double newDNA ) { DNA[k] = newDNA; };

	void Set_DNA_spread( const vector<double>& newDNA_spread ) { DNA_spread = newDNA_spread; };
	void Set_DNA_spread( const int k, const double newDNA_spread ) { DNA_spread[k] = newDNA_spread; };

	
	void Set_M_T( const vector<double>& newM_T ) { M_T = newM_T; };
	void Set_M_T( const int k, const double newM_T ) { M_T[k] = newM_T; };
	

	void Set_GAbsRate( const vector<double>& newGAbsRate ) { GAbsRate = newGAbsRate; }; 
	void Set_GAbsRate( const int k, const double newGAbsRate ) { GAbsRate[k] = newGAbsRate; }; 

	void Set_GConsRate( const vector<double>& newGConsRate ) { GConsRate = newGConsRate; };
	void Set_GConsRate( const int k, const double newGConsRate ) { GConsRate[k] = newGConsRate; };

	void Set_AAbsRate( const vector<double>& newAAbsRate ) { AAbsRate = newAAbsRate; }; 
	void Set_AAbsRate( const int k, const double newAAbsRate ) { AAbsRate[k] = newAAbsRate; }; 

	void Set_AConsRate( const vector<double>& newAConsRate ) { AConsRate = newAConsRate; };
	void Set_AConsRate( const int k, const double newAConsRate ) { AConsRate[k] = newAConsRate; };

	void Set_StoreFillRate( const vector<double>& newStoreFillRate ) { StoreFillRate = newStoreFillRate; };
	void Set_StoreFillRate( const int k, const double newStoreFillRate ) { StoreFillRate[k] = newStoreFillRate; };

	void Set_StoreConsRate( const vector<double>& newStoreConsRate ) { StoreConsRate = newStoreConsRate; };
	void Set_StoreConsRate( const int k, const double newStoreConsRate ) { StoreConsRate[k] = newStoreConsRate; };

	void Set_AcLRate( const vector<double>& newAcLRate ) { AcLRate = newAcLRate; };
	void Set_AcLRate( const int k, const double newAcLRate ) { AcLRate[k] = newAcLRate; };

	void Set_AcLOutRate( const vector<double>& newAcLOutRate ) { AcLOutRate = newAcLOutRate; };
	void Set_AcLOutRate( const int k, const double newAcLOutRate ) { AcLOutRate[k] = newAcLOutRate; };


	void Set_ATP_Ox( const vector<double>& newATP_Ox ) { ATP_Ox = newATP_Ox; };
	void Set_ATP_Ox( const int k, const double newATP_Ox ) { ATP_Ox[k] = newATP_Ox; };

	void Set_ATP_NOx( const vector<double>& newATP_NOx ) { ATP_NOx = newATP_NOx; };
	void Set_ATP_NOx( const int k, const double newATP_NOx ) { ATP_NOx[k] = newATP_NOx; };

	void Set_ATP2( const vector<double>& newATP2 ) { ATP2 = newATP2; };
	void Set_ATP2( const int k, const double newATP2 ) { ATP2[k] = newATP2; };

	void Set_ATP3( const vector<double>& newATP3 ) { ATP3 = newATP3; };
	void Set_ATP3( const int k, const double newATP3 ) { ATP3[k] = newATP3; };

	void Set_ConsATP( const vector<double>& newConsATP ) { ConsATP = newConsATP; };
	void Set_ConsATP( const int k, const double newConsATP ) { ConsATP[k] = newConsATP; };

	void Set_ConsATP_1( const vector<double>& newConsATP_1 ) { ConsATP_1 = newConsATP_1; };
	void Set_ConsATP_1( const int k, const double newConsATP_1 ) { ConsATP_1[k] = newConsATP_1; };

	void Set_ConsATP_2( const vector<double>& newConsATP_2 ) { ConsATP_2 = newConsATP_2; };
	void Set_ConsATP_2( const int k, const double newConsATP_2 ) { ConsATP_2[k] = newConsATP_2; };

	void Set_ConsATP_3( const vector<double>& newConsATP_3 ) { ConsATP_3 = newConsATP_3; };
	void Set_ConsATP_3( const int k, const double newConsATP_3 ) { ConsATP_3[k] = newConsATP_3; };

	void Set_ConsATP_4( const vector<double>& newConsATP_4 ) { ConsATP_4 = newConsATP_4; };
	void Set_ConsATP_4( const int k, const double newConsATP_4 ) { ConsATP_4[k] = newConsATP_4; };

	void Set_ConsATP_5( const vector<double>& newConsATP_5 ) { ConsATP_5 = newConsATP_5; };
	void Set_ConsATP_5( const int k, const double newConsATP_5 ) { ConsATP_5[k] = newConsATP_5; };

	void Set_ATPtot( const vector<double>& newATPtot ) { ATPtot = newATPtot; };
	void Set_ATPtot( const int k, const double newATPtot ) { ATPtot[k] = newATPtot; };



	// setters che inizializzano i contatori di produzione dell'ATP (normalmente vanno chiamati solo al momento della mitosi)
	void Set_ATPstart(  ) { for(unsigned long int k=0; k<ncells; k++) ATPstart[k] = (double)ATPp[k]; };	
	void Set_ATPstart( const unsigned long int k ) { ATPstart[k] = (double)ATPp[k]; };	

	void Set_ATPprod( const vector<double>& newATPprod ) { ATPprod = newATPprod; };
	void Set_ATPprod( const unsigned long int k, const double newATPprod ) { ATPprod[k] = newATPprod; };

	void Set_ATPcons( const vector<double>& newATPcons ) { ATPcons = newATPcons; };
	void Set_ATPcons( const unsigned long int k, const double newATPcons ) { ATPcons[k] = newATPcons; };
	



	// getters
	
	vector<unsigned long> Get_name() { return name; };
	unsigned long Get_name( const unsigned long int k ) { return name[k]; };

	vector<int> Get_mark() { return mark; };
	int Get_mark( const unsigned long int k ) { return mark[k]; };
	
	vector<CellType*> Get_type() { return type; };
	CellType* Get_type( const unsigned long int k ) { return type[k]; };
	
	vector<double> Get_Temperature() { return Temperature; };
	double Get_Temperature( const unsigned long int k ) { return Temperature[k]; };
	
	vector<CellPhase> Get_phase() { return phase; };	
	CellPhase Get_phase( const unsigned long int k ) { return phase[k]; };	

	vector<int> Get_death_condition() { return death_condition; };
	int Get_death_condition( const unsigned long int k ) { return death_condition[k]; };

	vector<float> Get_age() { return age; };
	float Get_age( const unsigned long int k ) { return age[k]; };

	vector<float> Get_phase_age() { return phase_age; };
	float Get_phase_age( const unsigned long int k ) { return phase_age[k]; };

	vector<float> Get_age_mother() { return age_mother; };
	float Get_age_mother( const unsigned long int k ) { return age_mother[k]; };

	vector<int> Get_n_mitosis() { return n_mitosis; };
	int Get_n_mitosis( const unsigned long int k ) { return n_mitosis[k]; };

	vector<vbl::BloodVessel> Get_BloodVesselVector() { return BloodVesselVector; };
	int Get_nbv() { return nbv; };

	vector<double> Get_x() { return x; };
	double Get_x( const unsigned long int k ) { return x[k]; };

	vector<double> Get_y() { return y; };
	double Get_y( const unsigned long int k ) { return y[k]; };

	vector<double> Get_z() { return z; };
	double Get_z( const unsigned long int k ) { return z[k]; };

	vector<double> Get_vx() { return vx; };
	double Get_vx( const unsigned long int k ) { return vx[k]; };

	vector<double> Get_vy() { return vy; };
	double Get_vy( const unsigned long int k ) { return vy[k]; };

	vector<double> Get_vz() { return vz; };
	double Get_vz( const unsigned long int k ) { return vz[k]; };

	vector<double> Get_r() { return r; };
	double Get_r( const unsigned long int k ) { return r[k]; };

	vector<double> Get_surface() { return surface; };
	double Get_surface( const unsigned long int k ) { return surface[k]; };

	vector<double> Get_volume() { return volume; };
	double Get_volume( const unsigned long int k ) { return volume[k]; };

	vector<double> Get_mass() { return mass; };
	double Get_mass( const unsigned long int k ) { return mass[k]; };


	vector<double> Get_volume_extra() { return volume_extra; };
	double Get_volume_extra( const unsigned long int k ) { return volume_extra[k]; };



	vector<int> Get_neigh() { return neigh; };
	int Get_neigh( const unsigned long int k ) { return neigh[k]; };

	// questi getters esistono solo nella forma per singole cellule (restituiscono vettori di lunghezza variabile)
	vector<int> Get_vneigh( const unsigned long int k ) { return vneigh[k]; };
	vector<double> Get_vdist( const unsigned long int k ) { return vdist[k]; };
	vector<double> Get_vcsurf( const unsigned long int k ) { return vcsurf[k]; };
	vector<double> Get_gnk( const unsigned long int k ) { return gnk[k]; };
	// fine della parte dei getters non standard

	vector<double> Get_contact_surf() { return contact_surf; };
	double Get_contact_surf( const unsigned long int k ) { return contact_surf[k]; };
	
	vector<bool> Get_isonCH() { return isonCH; };
	bool Get_isonCH( const unsigned long int k ) { return isonCH[k]; };

	vector<bool> Get_isonAS() { return isonAS; };
	bool Get_isonAS( const unsigned long int k ) { return isonAS[k]; };

	vector<int> Get_isonBV() { return isonBV; };
	int Get_isonBV( const unsigned long int k ) { return isonBV[k]; };

	vector<double> Get_env_surf() { return env_surf; }; 
	double Get_env_surf( const unsigned long int k ) { return env_surf[k]; }; 

	vector<double> Get_g_env() { return g_env; }; 
	double Get_g_env( const unsigned long int k ) { return g_env[k]; }; 

	vector<double> Get_bv_surf() { return bv_surf; }; 
	double Get_bv_surf( const unsigned long int k ) { return bv_surf[k]; }; 

	vector<double> Get_g_bv() { return g_bv; }; 
	double Get_g_bv( const unsigned long int k ) { return g_bv[k]; }; 


	vector<double> Get_M() { return M; }; 
	double Get_M( const unsigned long int k ) { return M[k]; }; 


	vector<double> Get_G() { return G; };
	double Get_G( const unsigned long int k ) { return G[k]; };

	vector<double> Get_G6P() { return G6P; };
	double Get_G6P( const unsigned long int k ) { return G6P[k]; };

	vector<double> Get_O2() { return O2; };
	double Get_O2( const unsigned long int k ) { return O2[k]; };

	vector<double> Get_store() { return store; };
	double Get_store( const unsigned long int k ) { return store[k]; };

	vector<double> Get_A() { return A; };
	double Get_A( const unsigned long int k ) { return A[k]; };

	vector<double> Get_AcL() { return AcL; };
	double Get_AcL( const unsigned long int k ) { return AcL[k]; };


	vector<double> Get_h() { return h; };
	double Get_h( const unsigned long int k ) { return h[k]; };
	
	vector<double> Get_pHi() { return pHi; };
	double Get_pHi( const unsigned long int k ) { return pHi[k]; };

	// double Get_H() { return H; };
	// double Get_CO2() { return CO2; };

	vector<double> Get_GAbsRate() { return GAbsRate; };
	double Get_GAbsRate( const unsigned long int k ) { return GAbsRate[k]; };

	vector<double> Get_GConsRate() { return GConsRate; };
	double Get_GConsRate( const unsigned long int k ) { return GConsRate[k]; };

	vector<double> Get_AAbsRate() { return AAbsRate; };
	double Get_AAbsRate( const unsigned long int k ) { return AAbsRate[k]; };

	vector<double> Get_AConsRate() { return AConsRate; };
	double Get_AConsRate( const unsigned long int k ) { return AConsRate[k]; };

	vector<double> Get_StoreFillRate() { return StoreFillRate; };
	double Get_StoreFillRate( const unsigned long int k ) { return StoreFillRate[k]; };

	vector<double> Get_StoreConsRate() { return StoreConsRate; };
	double Get_StoreConsRate( const unsigned long int k ) { return StoreConsRate[k]; };

	vector<double> Get_AcLRate() { return AcLRate; };
	double Get_AcLRate( const unsigned long int k ) { return AcLRate[k]; };

	vector<double> Get_AcLOutRate() { return AcLOutRate; };
	double Get_AcLOutRate( const unsigned long int k ) { return AcLOutRate[k]; };


	vector<double> Get_ATP_St() { return ATP_St; };
	double Get_ATP_St( const unsigned long int k ) { return ATP_St[k]; };

	vector<double> Get_ATP_Ox() { return ATP_Ox; };
	double Get_ATP_Ox( const unsigned long int k ) { return ATP_Ox[k]; };

	vector<double> Get_ATP_NOx() { return ATP_NOx; };
	double Get_ATP_NOx( const unsigned long int k ) { return ATP_NOx[k]; };

	vector<double> Get_ATP2() { return ATP2; };
	double Get_ATP2( const unsigned long int k ) { return ATP2[k]; };

	vector<double> Get_ATP3() { return ATP3; };
	double Get_ATP3( const unsigned long int k ) { return ATP3[k]; };

	vector<double> Get_ConsATP() { return ConsATP; };
	double Get_ConsATP( const unsigned long int k ) { return ConsATP[k]; };

	vector<double> Get_ConsATP_1() { return ConsATP_1; };
	double Get_ConsATP_1( const unsigned long int k ) { return ConsATP_1[k]; };

	vector<double> Get_ConsATP_2() { return ConsATP_2; };
	double Get_ConsATP_2( const unsigned long int k ) { return ConsATP_2[k]; };

	vector<double> Get_ConsATP_3() { return ConsATP_3; };
	double Get_ConsATP_3( const unsigned long int k ) { return ConsATP_3[k]; };

	vector<double> Get_ConsATP_4() { return ConsATP_4; };
	double Get_ConsATP_4( const unsigned long int k ) { return ConsATP_4[k]; };

	vector<double> Get_ConsATP_5() { return ConsATP_5; };
	double Get_ConsATP_5( const unsigned long int k ) { return ConsATP_5[k]; };

	vector<double> Get_ATPtot() { return ATPtot; };
	double Get_ATPtot( const unsigned long int k ) { return ATPtot[k]; };

	vector<double> Get_ATPp() { return ATPp; };	
	double Get_ATPp( const unsigned long int k ) { return ATPp[k]; };	

	vector<double> Get_ATPmin() { return ATPmin; };	
	double Get_ATPmin( const unsigned long int k ) { return ATPmin[k]; };	


	vector<double> Get_ATPstart() { return ATPstart; };
	double Get_ATPstart( const unsigned long int k ) { return ATPstart[k]; };

	vector<double> Get_ATPprod() { return ATPprod; };
	double Get_ATPprod( const unsigned long int k ) { return ATPprod[k]; };

	vector<double> Get_ATPcons() { return ATPcons; };
	double Get_ATPcons( const unsigned long int k ) { return ATPcons[k]; };


	vector<double> Get_G_extra() { return G_extra; };
	double Get_G_extra( const unsigned long int k ) { return G_extra[k]; };

	vector<double> Get_A_extra() { return A_extra; };
	double Get_A_extra( const unsigned long int k ) { return A_extra[k]; };

	vector<double> Get_AcL_extra() { return AcL_extra; };
	double Get_AcL_extra( const unsigned long int k ) { return AcL_extra[k]; };

	
	vector<double> Get_pH() { return pH; };
	double Get_pH( const unsigned long int k ) { return pH[k]; };

	// double Get_H_extra() { return H_extra; };
	// double Get_CO2_extra() { return CO2_extra; };
	
	vector<double> Get_SensO2() { return SensO2; };
	double Get_SensO2( const unsigned long int k ) { return SensO2[k]; };

	vector<double> Get_ConsO() { return ConsO; };
	double Get_ConsO( const unsigned long int k ) { return ConsO[k]; };

	// double Get_ProdCO2() { return ProdCO2; };

	vector<double> Get_DNA_spread() { return DNA_spread; };
	double Get_DNA_spread( const unsigned long int k ) { return DNA_spread[k]; };

	
	vector<double> Get_M_T() { return M_T; };
	double Get_M_T( const int k ) { return M_T[k]; };


	vector<double> Get_DNA() { return DNA; };
	double Get_DNA( const unsigned long int k ) { return DNA[k]; };

	vector<double> Get_ConcS() { return ConcS; };
	double Get_ConcS( const unsigned long int k ) { return ConcS[k]; };


	vector<double> Get_protein() { return protein; };
	double Get_protein( const unsigned long int k ) { return protein[k]; };

	vector<double> Get_prot_rate() { return prot_rate; };
	double Get_prot_rate( const unsigned long int k ) { return prot_rate[k]; };

	vector<double> Get_DNA_rate() { return DNA_rate; };
	double Get_DNA_rate( const unsigned long int k ) { return DNA_rate[k]; };

	vector<double> Get_pRb() { return pRb; };
	double Get_pRb( const unsigned long int k ) { return pRb[k]; };

	vector<double> Get_cyclinD() { return cyclinD; };
	double Get_cyclinD( const unsigned long int k ) { return cyclinD[k]; };

	vector<double> Get_cyclinE() { return cyclinE; };
	double Get_cyclinE( const unsigned long int k ) { return cyclinE[k]; };

	vector<double> Get_cyclinX() { return cyclinX; };
	double Get_cyclinX( const unsigned long int k ) { return cyclinX[k]; };



	vector<double> Get_NpRbk() { return NpRbk; };
	double Get_NpRbk( const unsigned long int k ) { return NpRbk[k]; };

	

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
        if(M[k] < -(numeric_limits<double>::epsilon( )))
            {
            cout << "Inconsistent value of M in the cell " << name[k] << ": M=" << M[k] << endl;
            return_code = -1;
            }
        if(volume[k] < type[k]->Vmin)
            {
            cout << "Inconsistent volume value in the cell " << name[k] << ": Vmin=" << type[k]->Vmin << endl;
            return_code = -2;
            }
        if(ATPp[k] < -(numeric_limits<double>::epsilon( )))
            {
            cout << "Inconsistent value of ATPp in the cell " << name[k] << ": ATPp=" << ATPp[k] << endl;
            return_code = -3;
            }
        if( fabs( ATPmin[k] - ((type[k]->fATPmin)*((type[k]->C2) * M[k])/(type[k]->C1)) ) > (numeric_limits<double>::epsilon( )) )
            {
            cout << "Inconsistent value of ATPmin in the cell " << name[k] << ": ATPmin=" << ATPmin[k] << endl;
            cout << scientific << "ATPmin[k] = " << ATPmin[k] << 
                                  "type[k]->fATPminVmin " << type[k]->fATPmin <<
                                  "type[k]->C2 = " <<     type[k]->C2<<
                                  "type[k]->C1 = " <<     type[k]->C1<<
                                  "M[k] = "        <<     M[k] << endl;
            return_code = -4;
            }
        if( fabs( volume[k] - (type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1) )/volume[k] > 1.e-6 )
            {
            cout << "Inconsistent values ​​of ATPp, M, volume, in the cell " << name[k] << endl;
            cout << scientific << "volume = " << volume[k] << " = type->Vmin * (1.+DNA) + type->C2 * M + ATPp * type->C1 = " << (type[k]->Vmin * (1.+DNA[k]) + type[k]->C2 * M[k] + ATPp[k] * type[k]->C1) << endl;
            cout << "DNA: " << DNA[k] << endl;
            cout << "M: " << M[k] << endl;
            cout << "ATPp: " << ATPp[k] << endl;
            return_code = -5;
            }
        }
		
	if(return_code < 0) 
		{
		cout << endl;
		count++;
		}

	
	// dopo 20 errori di questo tipo il programma si ferma
	if(count > 20) exit(-1);
	
	
	
	return return_code;

	}
	

// *** fine dei metodi per la gestione della parte biofisica *** 

};


// distanza tra la cellula con indice j e la cellula con indice k
//
inline double CellsSystem::Distance(const int j, const int k)
{
	return sqrt( SQR(x[j]-x[k]) + SQR(y[j]-y[k]) + SQR(z[j]-z[k]) );
}

#endif //#ifndef CELLSSYSTEM_H

