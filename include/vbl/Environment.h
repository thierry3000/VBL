// 
// definizione della classe Enviroment che contiene le informazioni sull'ambiente in cui sono immerse le cellule
//
// EM 20/1/2008
//
// **********************************************************************************
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H  // header guard
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>      // std::setprecision  std::scientific
namespace vbl{
class Environment 
{
// friends

	friend class CellsSystem;
	
// update dell'ambiente  (env e' l'ambiente da modificare, env_delta e' la lista dei cambiamenti)
	friend void EnvironmentChange(Environment& env, Environment& env_delta);
	
// funzione di stampa di tutti i dati della classe
	friend std::ostream& operator<<(std::ostream& s, Environment& cEnvironment);
	
private:

	double T;			// temperatura dell'ambiente
	double G;			// quantita' (kg) di glucosio ambientale
	double O2;			// quantita' (kg) di ossigeno ambientale
	double CO2;			// quantita' (kg) di anidride carbonica ambientale
	double A;			// quantita' (kg) di altri nutrienti nell'ambiente
	double AcL;			// quantita' (kg) di acido lattico nell'ambiente
	double pH;			// pH ambientale
	double xmin;		// definizione della forma dell'ambiente (un parallelepipedo con spigoli nelle posizioni definite da xmin, ... )
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	double volume0;		// volume iniziale dell'ambiente
	double volume;		// volume dell'ambiente
	double DoseRate;	// dose irraggiamento per unita' di tempo (idealmente Sv, attualmente Gy, perche' non c'e' una relative effectiveness definita)

public:

// costruttori
// costruttore di default (corrisponde al mezzo di coltura standard senza circolazione del fluido)
	Environment();
// costruttore con un ambiente predefinito
	Environment(EnvironmentTypeSelector environmentType);
// costruttore che prende i dati da un file
	Environment(const std::string filename);
// costruttore copia
	Environment(const Environment& e);
	
// stampa dell'ambiente in formato spreadsheet
	void PrintEnvironmentData(std::ofstream& stream, long int nrec);
	
// nessun distruttore, utilizzo il default del compilatore
	
// overloaded =
	Environment& operator=(const Environment& et);
	
// setters
	void SetEnvironmentT(const double newT) { T = newT; };
	void SetEnvironmentG(const double newG) { G = newG; };
	void SetEnvironmentO2(const double newO2) { O2 = newO2; };
	void SetEnvironmentCO2(const double newCO2) { CO2 = newCO2; };
	void SetEnvironmentA(const double newA) { A = newA; };
	void SetEnvironmentpH() { pH = 7.5443-(AcL/volume)/BufCapEnv;  };	
	void SetEnvironmentAcL(const double newAcL) { AcL = newAcL; };
	void Set_xmin(const double newxmin) { xmin = newxmin; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	void Set_xmax(const double newxmax) { xmax = newxmax; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	void Set_ymin(const double newymin) { ymin = newymin; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	void Set_ymax(const double newymax) { ymax = newymax; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	void Set_zmin(const double newzmin) { zmin = newzmin; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	void Set_zmax(const double newzmax) { zmax = newzmax; volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin); };
	// void SetEnvironmentvolume0(double newvolume0) { volume0 = newvolume0; };
	void SetEnvironmentvolume(double newvolume) { volume = newvolume; };
	void SetEnvironmentDoseRate(double newDoseRate) { DoseRate = newDoseRate; };

// getters
	double GetEnvironmentT() { return T; }; 
	double GetEnvironmentG() { return G; }; 
	double GetEnvironmentO2() { return O2; }; 
	double GetEnvironmentCO2() { return CO2; }; 
	double GetEnvironmentA() { return A; }; 
	double GetEnvironmentpH() { return pH; }; 
	double GetEnvironmentAcL() { return AcL; }; 
	double Get_xmin() { return xmin; };
	double Get_xmax() { return xmax; };
	double Get_ymin() { return ymin; };
	double Get_ymax() { return ymax; };
	double Get_zmin() { return zmin; };
	double Get_zmax() { return zmax; };
	double GetEnvironmentvolume0() { return volume0; }; 
	double GetEnvironmentvolume() { return volume; }; 
	double GetEnvironmentDoseRate() { return DoseRate; }; 

// read and write binario
	void WriteEnvironment( std::ofstream& stream );
	void ReadEnvironment( std::ifstream& stream );


};
}//namespace vbl{
#endif //#ifndef ENVIRONMENT_H
