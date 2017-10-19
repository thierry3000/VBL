// 
// contiene le definizioni delle classi che rappresentano segnali ambientali
//
// EM 19/2/2008
//
// **********************************************************************************

#ifndef ENVIRONMENTSIGNAL_H
#define ENVIRONMENTSIGNAL_H  // header guard

#include <iostream>
#include <fstream>
#include <cmath>
class EnvironmentalSignal
{

// overloaded << 
friend std::ostream& operator<<(std::ostream& stream, EnvironmentalSignal& cSignal);

private:
	EnvironmentalSignalType type;
	double ampMin;
	double ampMax;
	double period;
	double tstart;
	double tstop;
	double tON;
	double tOFF;

public:

// costruttori
	EnvironmentalSignal() { type = NullSignal; ampMin = ampMax = period = tstart = tstop = tON = tOFF = 0.; };
	EnvironmentalSignal( const double newamp ) { type = ConstantSignal; ampMin = newamp; };
	EnvironmentalSignal( const EnvironmentalSignalType ctype, const double cAmpMin, const double cAmpMax, const double ctA, const double ctB );
	EnvironmentalSignal( const EnvironmentalSignalType ctype, const double cAmpMin, const double cAmpMax, const double ctA, const double ctB, const double ctC );

// value 
	double SignalValue(const double t);
	
// integrale in nell'intervallo (t1,t2)
	double SignalIntegral(const double t1, const double t2);
		
// setters
	void Set_ampMin(const double campMin) { ampMin = campMin; };
	void Set_ampMax(const double campMax) { ampMax = campMax; };
	void Set_period(const double cperiod) { period = cperiod; };
	void Set_tstart(const double ctstart) { tstart = ctstart; };
	void Set_tstop(const double ctstop) { tstop = ctstop; };
	void Set_tON(const double ctON) { tON = ctON; period = tON+tOFF; };
	void Set_tOFF(const double ctOFF) { tOFF = ctOFF; period = tON+tOFF; };
	
// getters
	EnvironmentalSignalType Get_type() { return type; };
	double Get_ampMin() { return ampMin; };
	double Get_ampMax() { return ampMax; };
	double Get_period() { return period; };
	double Get_tstart() { return tstart; };
	double Get_tstop() { return tstop; };
	double Get_tON() { return tON; };
	double Get_tOFF() { return tOFF; };

// overloaded =
	EnvironmentalSignal& operator=(const EnvironmentalSignal& es);
	
// read and write binario
	void WriteEnvironmentalSignal( std::ofstream& stream);
	void ReadEnvironmentalSignal( std::ifstream& stream);

	
	
};

#endif //#ifndef ENVIRONMENTSIGNAL_H
