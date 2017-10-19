// 
// contiene i metodi per le classi che rappresentano segnali ambientali
//
// EM 19/2/2008
//
// **********************************************************************************

#include "sim.h"
#include "EnvironmentalSignals.h"

// costruttore (caso ConstantSignal, SineSignal, Pulse, e per il momento UserDefined, che pero' e' vuoto)
EnvironmentalSignal::EnvironmentalSignal( const EnvironmentalSignalType ctype, const double cAmpMin, const double cAmpMax, const double ctA, const double ctB )
{
	type = ctype;
	
	if( cAmpMin > cAmpMax) std::cout << "\nWrong amplitude assignement in EnvironmentalSignal::EnvironmentalSignal\n" << std::endl;
	
	switch(type)
		{
		case ConstantSignal:
		ampMin = cAmpMin;
		break;
		
		case SineSignal:
		ampMin = cAmpMin;
		ampMax = cAmpMax;
		tstart = ctA;
		period = ctB;
		break;
				
		case SquareSignal:
		std::cout << "EnviromentalSignal: wrong constructor\n";
		break;
		
		case Pulse:
		if( tstop < tstart) std::cout << "\nWrong timing assignement in EnvironmentalSignal::EnvironmentalSignal\n" << std::endl;
		ampMin = cAmpMin;
		ampMax = cAmpMax;
		tstart = ctA;
		tstop = ctB;
		break;
		
		case UserDefined:
		break;
		
		default:
		break; 
		 
		}

}

// costruttore (caso ConstantSignal, SquareSignal)
EnvironmentalSignal::EnvironmentalSignal( const EnvironmentalSignalType ctype, const double cAmpMin, const double cAmpMax, const double ctA, const double ctB, const double ctC )
{
	type = ctype;

	switch(type)
		{
		case ConstantSignal:
		ampMin = cAmpMin;
		break;
		
		case SquareSignal:
		ampMin = cAmpMin;
		ampMax = cAmpMax;
		tstart = ctA;
		tON = ctB;
		tOFF = ctC;		
		period = tON + tOFF;
		break;
		
		default:
		std::cout << "EnviromentalSignal: wrong constructor\n";
		break;
		
		}
	
}

// calcola il valore attuale del segnale
double EnvironmentalSignal::SignalValue(const double t)
{
	double val=0;
	double tt;
	
	switch(type)
		{
		case NullSignal:
		val = 0.;
		break;
		
		case ConstantSignal:
		val = ampMin;
		break;
		
		case SineSignal:
		val = ampMin + (ampMax-ampMin) * 0.5*( 1. + sin( 2*PI*(t-tstart)/period ) );
		break;
		
		case SquareSignal:
		tt = (t-tstart) - period * floor((t-tstart)/period);
		val = tt<tON ? ampMax : ampMin;
		break;
		
		case Pulse:
		val = (t>tstart && t<=tstop) ? ampMax : ampMin;
		break;
		
		case UserDefined:
		break;
		
		default: 
		val = 0.;
		break; 
		
		}
	
	return val;

}


// calcola il valore dell'integrale del segnale nell'intervallo (t1,t2)  (ATTENZIONE: si assume t2 >= t1 )
//
double EnvironmentalSignal::SignalIntegral(const double t1, const double t2)
{
	double val=0;
	// double tt;
	double n1;
	double int1;
	double tover1;
	double n2;
	double int2;
	double tover2;
	
	double t1a = t1;
	double t2a = t2;
	
	switch(type)
		{
		case NullSignal:
		val = 0.;
		break;
		
		case ConstantSignal:
		val = ampMin*(t2-t1);
		break;
		
		case SineSignal:
		val = 0.5*(ampMax+ampMin)*(t2-t1) - 0.5*(ampMax-ampMin)*(period/(2.*PI))*( cos( 2.*PI*(t2-tstart)/period ) - cos( 2.*PI*(t1-tstart)/period ) );
		break;
		
		case SquareSignal:
		val = ampMin*(t2-t1);
		
		t1a = t1+(period-tstart);
		
		n1 = floor(t1a/period);
		int1 = (ampMax - ampMin)*n1*tON;
		tover1 = (t1a - n1*period);
		if( tover1 < tON )
			int1 += (ampMax - ampMin)*tover1;
		else 
			int1 += (ampMax - ampMin)*tON;
			
		t2a = t2+(period-tstart);

		n2 = floor(t2a/period);
		int2 = (ampMax - ampMin)*n2*tON;
		tover2 = (t2a - n2*period);
		if( tover2 < tON )
			int2 += (ampMax - ampMin)*tover2;
		else 
			int2 += (ampMax - ampMin)*tON;
			
		val += (int2-int1);
		break;
		
		case Pulse:
		if( tstop < t1 || tstart > t2 )						// impulso fuori dall'intervallo (t1,t2)
			val = ampMin*(t2-t1);
		else if( tstart <= t1 && tstop >= t2 )				// l'impulso copre completamente l'intervallo
			val = ampMax*(t2-t1);
		else if( tstart >= t1 && tstart <= t2 && tstop >= t1 && tstop <= t2 )	// l'impulso e' completamente contenuto nell'intervallo
			val = ampMin*(t2-t1) + (ampMax - ampMin)*(tstop-tstart);
		else if( tstart <= t1 && tstop >= t1 && tstop <= t2 )	// impulso a cavallo del bordo sinistro dell'intervallo
			val = ampMin*(t2-t1) + (ampMax - ampMin)*(tstop-t1);
		else if( tstart >= t1 && tstart <= t2 && tstop >= t2 )	// impulso a cavallo del bordo destro dell'intervallo
			val = ampMin*(t2-t1) + (ampMax - ampMin)*(t2-tstart);
		break;
		
		case UserDefined:
		break;
		
		default: 
		val = 0.;
		break; 
		
		}
	
	return val;

}

// overloaded =
EnvironmentalSignal& EnvironmentalSignal::operator=(const EnvironmentalSignal& e)
{

	type = e.type;
	ampMin = e.ampMin;
	ampMax = e.ampMax;
	period = e.period;
	tstart = e.tstart;
	tstop = e.tstop;
	tON =  e.tON;
	tOFF =  e.tOFF;
	
	return *this;
}

// overloading dell'operatore di output su file
std::ostream& operator<<(std::ostream& stream, EnvironmentalSignal& cSignal)
{

	EnvironmentalSignalType iType; 
	
	const std::string SType[] = {"segnale nullo", "segnale costante", "sinusoide", "onda quadra", "singolo impulso rettangolare"};
	
	iType = cSignal.Get_type();
	
	stream << "Tipo del segnale (codice): " << (int)iType << std::endl;
	stream << "Tipo del segnale: " <<  SType[ (int)iType + 1 ]  << std::endl;

	switch(iType)
		{
		case ConstantSignal:
		stream << "Ampiezza: " << cSignal.Get_ampMin() << std::endl;
		break;
		
		case SineSignal:
		stream << "Ampiezza minima: " << cSignal.Get_ampMin() << std::endl;
		stream << "Ampiezza massima: " << cSignal.Get_ampMax() << std::endl;
		stream << "Periodo (s): " << cSignal.Get_period() << std::endl;
		stream << "tstart (s): " << cSignal.Get_tstart() << std::endl;
		break;
		
		case SquareSignal:
		stream << "Ampiezza minima: " << cSignal.Get_ampMin() << std::endl;
		stream << "Ampiezza massima: " << cSignal.Get_ampMax() << std::endl;
		stream << "tstart (s): " << cSignal.Get_tstart() << std::endl;
		stream << "tON (s): " << cSignal.Get_tON() << std::endl;
		stream << "tOFF (s): " << cSignal.Get_tOFF() << std::endl;
		break;
		
		case Pulse:
		stream << "Ampiezza minima: " << cSignal.Get_ampMin() << std::endl;
		stream << "Ampiezza massima: " << cSignal.Get_ampMax() << std::endl;
		stream << "tstart (s): " << cSignal.Get_tstart() << std::endl;
		stream << "tstop (s): " << cSignal.Get_tstop() << std::endl;
		break;
		
		case UserDefined:
		break;
		
		default: 
		break; 
		
		}
					
	return stream;

}

void EnvironmentalSignal::WriteEnvironmentalSignal( std::ofstream& stream)
{
	stream.write( (char*)(&ampMin), sizeof(double) );
	stream.write( (char*)(&ampMax), sizeof(double) );
	stream.write( (char*)(&period), sizeof(double) );
	stream.write( (char*)(&tstart), sizeof(double) );
	stream.write( (char*)(&tstop), sizeof(double) );
	stream.write( (char*)(&tON), sizeof(double) );
	stream.write( (char*)(&tOFF), sizeof(double) );

}

void EnvironmentalSignal::ReadEnvironmentalSignal( std::ifstream& stream)
{
	stream.read( (char*)(&ampMin), sizeof(double) );
	stream.read( (char*)(&ampMax), sizeof(double) );
	stream.read( (char*)(&period), sizeof(double) );
	stream.read( (char*)(&tstart), sizeof(double) );
	stream.read( (char*)(&tstop), sizeof(double) );
	stream.read( (char*)(&tON), sizeof(double) );
	stream.read( (char*)(&tOFF), sizeof(double) );

}

