//
// header principale del programma di simulazione
//
// EM 19/1/2008
//
// **********************************************************************************


// **********************************************************************************
//
// standard includes
//
// si noti che in tutto il programma si assume sempre che il namespace std sia definito
// in modo da rendere il codice piu' leggibile
//
// **********************************************************************************

#ifndef SIM_H
#define SIM_H  // header guard
namespace vbl{
// #include <iostream>			// inclusione della libreria standard di I/O
// #include <iomanip>			// inclusione dei manipolatori per l'I/O
// #include <fstream>			// dichiarazione headers che servono per l'I/O su file 
// #include <string>			// stringhe della STL
// #include <limits>			// limiti numerici e funzioni collegate
// 
// #include <list>				// inclusione della libreria per la manipolazione delle liste STL
// #include <vector>			// inclusione della libreria per la manipolazione dei vettori STL

// using namespace std;
// 
// #include <cmath>			// inclusione della libreria matematica standard del C
// #include <cstdio>			// inclusione della libreria stdio del C
// #include <cstdlib>			// inclusione della stdlib del C
// #include <time.h>			// inclusione delle utilities di timing del C
// #include <string.h>			// inclusione delle utilities di stringhe del C
// 
// #ifdef _OPENMP
// #include <omp.h>			// inclusione delle utilities di OpenMP
// #endif

// #include "nr.h"			// header di Numerical Recipes


// **********************************************************************************
//
// general inline functions
//
// **********************************************************************************

// quadrato di un long double
inline long double SQR(const long double x) { return x*x; }
// quadrato di un double
inline double SQR(const double x) { return x*x; }
// quadrato di un long int
inline long int SQR(const long int x) { return x*x; }
// quadrato di un int
inline long int SQR(const int x) { return x*x; }


// **********************************************************************************
//
// costanti matematiche e chimico-fisiche non modificabili
//
// **********************************************************************************

const double PI = 3.1415926;

const double NAV = 6.022e23;           // costante di Avogadro
const double Faraday = 96485.34;		// costante di Faraday in C

// ATTENZIONE! TUTTI I PESI MOLECOLARI SONO IN GRAMMI !!
const double PM_G = 180.;              // peso molecolare glucosio in g
const double PM_ATP = 507.;            // peso molecolare ATP in g
const double PM_O2 = 32.;              // peso molecolare dell'O2 in g
const double PM_AcL = 90.;             // peso molecolare dell'acido lattico in g
const double PM_G6P = 270;             // peso molecolare del G6P in g
const double PM_H = 1.;                // peso molecolare dell'idrogeno atomico in g
const double PM_pRb = 110000.;         // peso molecolare della pRb in g
const double PM_A = 146.;              // peso molecolare degli altri nutrienti (glutammina) in g
const double PM_prot = 66476.;         // peso molecolare medio delle proteine (albumina) in g
const double PM_cyclinD = 33729.;      // peso molecolare medio della ciclina D in g
const double PM_cyclinE = 47077.;      // peso molecolare medio della ciclina e in g
const double PM_cyclinX = 52000.;      // peso molecolare medio delle cicline A + B in g

const double Gibbs_ATP = 3.1e4;		// energia libera di Gibbs per l'idrolisi dell'ATP in J/mole

// coefficienti di diffusione in acqua e negli spazi extracellulari: DA RICONTROLLARE CON ATTENZIONE!!!

const double Diff_W_G = 7.e2;          // coefficiente di diffusione del glucosio in acqua in micron^2/s
const double Diff_ES_G = 2.e2;         // coefficiente di diffusione del glucosio negli spazi extracellulari in micron^2/s
const double Diff_Env_G = Diff_ES_G;   // the diffusion coefficient in the environment is 
									   // equal to that in the extracellular spaces in
									   // this simulation (same below)
const double Diff_BV_G = Diff_W_G;	   // the diffusion coefficient within blood vessels 
									   // is the same as in water in this simulation 
									   // (same below)

const double Diff_W_A = 3.e2;          // coefficiente di diffusione della glutammina in acqua in micron^2/s
const double Diff_BV_A = Diff_W_A;
const double Diff_ES_A = 7.e1;         // coefficiente di diffusione della glutammina negli spazi extracellulari in micron^2/s
const double Diff_Env_A = Diff_ES_A;

const double Diff_W_AcL = 3.e2;		// coefficiente di diffusione dell'acido lattico in acqua in micron^2/s
const double Diff_BV_AcL = Diff_W_AcL;
const double Diff_ES_AcL = 1.4e1;      // coefficiente di diffusione dell'acido lattico negli spazi extracellulari in micron^2/s
const double Diff_Env_AcL = Diff_ES_AcL;

const double Diff_W_O2 = 3.2e3;		// coefficiente di diffusione dell'ossigeno in acqua in micron^2/s
const double Diff_BV_O2 = Diff_W_O2;
const double Diff_ES_O2 = 7e1;         // coefficiente di diffusione dell'ossigeno nelle cellule in micron^2/s
const double Diff_Env_O2 = Diff_ES_O2;

const double Diff_W_H = 7e3;           // coefficiente di diffusione dei protoni in acqua in micron^2/s
const double Diff_ES_H = 7e3;          // coefficiente di diffusione dei protoni negli spazi extracellulari in micron^2/s

// const double DiffCO2 = 2.2e3;		// coefficiente di diffusione dell'anidride carbonica sia nel mezzo nutritivo che nelle cellule in micron^2/s


// **********************************************************************************
//
// definizioni relative alle fasi cellulari
//
// **********************************************************************************

// numero di tipi di fase validi
const int Nphase = 7;

// nomi delle fasi cellulari
enum CellPhase {undefined_phase = -1, G0_phase = 0, G1m_phase, G1p_phase, S_phase, G2_phase, M_phase, dead};


// **********************************************************************************
//
// altri enums
//
// **********************************************************************************

// tipo di ambiente
enum EnvironmentTypeSelector {NullEnvironment=-1, Standard, TerminalInput};

// tipi di segnale 
enum EnvironmentalSignalType { NullSignal=-1, ConstantSignal, SineSignal, SquareSignal, Pulse, UserDefined };

// tipo di simulazione
enum SimType { Disperse=0, Full3D, FixedConfig };

// tipo di distribuzione iniziale delle cellule
enum DistType { OneCell=1, Dist2D, Dist3D };

// enum per il timer
enum timer_button { Start_timer=0, Start_intertime, Stop_intertime, Restart_timer, Clear_intertime };

// **********************************************************************************
// 
// parametri che regolano il comportamento del programma
//
// **********************************************************************************

// parametri per l'inizializzazione dello stato cellulare
const double ATPp_STANDARD = 2.6;      // pg ATPp iniziale (grossolanamente uguale a meta' del valore all'equilibrio per cellule da 5 micron di raggio)
const double pRb_STANDARD = 0.1;		// pg pRb iniziale (grossolanamente uguale al valore all'equilibrio)
const double pHi_STANDARD = 7.2;		// pH interno alle cellule in condizioni normali
const double alpha_R_default[Nphase] = { 0, 0.43, 0.43, 0.17, 0.77, 0.77, 0 };	// valori di default del coefficiente alpha della legge lineare-quadratica
const double beta_R_default[Nphase] = { 0, 0.017, 0.017, 0.02, 0, 0, 0 };		// valori di default del coefficiente beta della legge lineare-quadratica

// qui  si definiscono le variabili ambientali iniziali
const double T_ENV = 37;               // temperatura ambientale in °C
//const double O2st = 3.392e-5;		// pg*micron^-3
const double O2st = 7e-6;              // pg*micron^-3
const double XMIN_ENV = -500.;         // definizione del parallelepipedo standard (micron) (corrisponde a 1 microlitro)
const double XMAX_ENV = 500.;
const double YMIN_ENV = -500.;
const double YMAX_ENV = 500.;
const double ZMIN_ENV = -500.;
const double ZMAX_ENV = 500.;

const double G_ENV = 0.9e-3;           // pg*micron^-3 concentrazione glucosio

const double O_ENV = 7e-6;             // pg*micron^-3 concentrazione ossigeno (atmosfera standard) (ATTENZIONE, QUESTO E' ANCHE UGUALE A O2St, controllare in caso di cambiamenti)

const double CO2_ENV = 3.4e-4;         // kg*m^-3 concentrazione anidride carbonica (atmosfera standard)

const double A_ENV = 0.4e-3;           // pg*micron^-3 concentrazione altri nutrienti
const double A_CELL = 1.e-6;           // pg quantita' iniziale degli altri nutrienti nelle cellule
const double AcL_ENV = 0.;             // pg*micron^-3 concentrazione esterna dell'acido lattico
const double STOCK_MAX = 0.018;		// pg contenuto massimo dello stock di glucosio per cellula

// T.F. revised in agreement with Roberto and Edoardo
//const double BufCapEnv = 0.19953e-3;   // buffering capacity dell'ambiente (pg/micron^3)

const double BufCapEnv = 0.584e-3;     // buffering capacity dell'ambiente (pg/micron^3)

const double DOSERATE_ENV = 0;         // radiazione ambientale (Gy/s)

const double VISCOSITY_ENV = 2.e11;     // pg/(micron s) viscosita' dell'ambiente (test di viscosità ambientale molto elevata, uguale a quella cellulare)


const int MAXREPEATS = 10000;                // numero max di ripetizioni del metodo di diffusione
const int NCONV_TEST = 12;                  // numero di tests di convergenza (v. CellSystem-B)

const int TOL = 1e-7;                       // tolleranza minima (in pg) sulla determinazione della quantità di sostanze che determinano la condizione di stop del loop in CellSystem::Diff
//const int TOL = 0.;                       // tolleranza minima (in pg) sulla determinazione della quantità di sostanze che determinano la condizione di stop del loop in CellSystem::Diff


const double ALPHAVALUE = 100.;		// alpha per il calcolo delle alpha shapes uguale a (10micron)^2

const long int RESERVE = 2000000;           // definizione della riserva dinamica per i vettori
const long int RESERVE_BV = 2000000;        // dynamic reserve for the BloodVesselVector

// standard concentration values for blood vessels

const double G_BV = 0.9e-3;           // pg*micron^-3 concentrazione glucosio

const double O2_BV = 7e-6;            // pg*micron^-3 concentrazione ossigeno (atmosfera standard) (ATTENZIONE, QUESTO E' ANCHE UGUALE A O2St, controllare in caso di cambiamenti)

const double CO2_BV = 3.4e-4;         // kg*m^-3 concentrazione anidride carbonica (atmosfera standard)

const double A_BV = 0.4e-3;           // pg*micron^-3 concentrazione altri nutrienti

}//namespace vbl{
#endif// header guard
