// 
// definizione dei metodi per la classe Enviroment che contiene le informazioni sull'ambiente in cui sono immerse le cellule
//
// EM 20/1/2008
//
// **********************************************************************************

#include "sim.h"

#include "InputFromFile.h"
#include "Environment.h"
#include <cstdlib>
using namespace vbl;

// costruttore di default (corrisponde al mezzo di coltura standard senza circolazione del fluido)
Environment::Environment()
{
	xmin = XMIN_ENV;
	xmax = XMAX_ENV;
	ymin = YMIN_ENV;
	ymax = YMAX_ENV;
	zmin = ZMIN_ENV;
	zmax = ZMAX_ENV;
	
	volume0 = (XMAX_ENV-XMIN_ENV)*(YMAX_ENV-YMIN_ENV)*(ZMAX_ENV-ZMIN_ENV);
	
	T = T_ENV;
	G = G_ENV * volume0;
	O2 = O_ENV * volume0;
	A = A_ENV * volume0;
	AcL = AcL_ENV * volume0;
	pH =  7.5443-(AcL/volume0)/BufCapEnv;			// pH ambientale (si assume la linearita' della relazione
													// che lega la conc. di acido lattico al pH)
	volume = volume0;
	DoseRate = 0;									// normalmente non c'e' irraggiamento
}


// costruttore con un ambiente predefinito
// 
// per ora sono definiti: 
// -1. Ambiente vuoto (i parametri sono tutti uguali a 0)
//  0. Standard (equivalente al costruttore default)
//  1. TerminalInput (input dei parametri da terminale)
// 
Environment::Environment(EnvironmentTypeSelector environmentType)
{
	double Gc, O2c, CO2c, Ac, AcLc;
	double v0;
	double dose;
	
	switch(environmentType)
		{
		case NullEnvironment:
			T = 0;
			G = 0;
			O2 = 0;
			CO2 = 0;
			A = 0;
			AcL = 0;
			pH =  0;
			xmin = 0;
			xmax = 0;
			ymin = 0;
			ymax = 0;
			zmin = 0;
			zmax = 0;
			volume0 = 0;
			volume = 0;
			DoseRate = 0;
			break;
		
		case Standard:
			xmin = XMIN_ENV;
			xmax = XMAX_ENV;
			ymin = YMIN_ENV;
			ymax = YMAX_ENV;
			zmin = ZMIN_ENV;
			zmax = ZMAX_ENV;
			
			volume0 = (XMAX_ENV-XMIN_ENV)*(YMAX_ENV-YMIN_ENV)*(ZMAX_ENV-ZMIN_ENV);
			
			T = T_ENV;
			G = G_ENV * volume0;
			O2 = O_ENV * volume0;
			CO2 = CO2_ENV * volume0;
			A = A_ENV * volume0;
			AcL = AcL_ENV * volume0;
			pH =  7.5443-(AcL/volume0)/BufCapEnv;
			volume = volume0;
			DoseRate = DOSERATE_ENV;
			break;
			
		case TerminalInput:
			std::cout << "Input dei parametri ambientali (dando in input un numero negativo si sceglie il valore di default)\n" << std::endl;
			
			std::cout << "Temperatura dell'ambiente (deg C): ";
			std::cin >> T;
			if(T < 0)
				T = T_ENV;
			std::cout << "Temperatura dell'ambiente: " << T << " C"<< std::endl;
				
			std::cout << "Concentrazione del glucosio (pg/micron^3): ";
			std::cin >> Gc;
			if(Gc < 0) 
				Gc = G_ENV;
			std::cout << "La concentrazione del glucosio e' " << Gc << " pg/micron^3" << std::endl;
			
			std::cout << "Concentrazione dell' ossigeno (O2) (pg/micron^3): ";
			std::cin >> O2c;
			if(O2c < 0) 
				O2c = O_ENV;
			std::cout << "La concentrazione dell' ossigeno (O2) e' " << O2c << " pg/micron^3" <<  std::endl;
			
			std::cout << "Concentrazione dell' anidride carbonica (CO2) (pg/micron^3): ";
			std::cin >> CO2c;
			if(CO2c < 0) 
				CO2c = CO2_ENV;
			std::cout << "La concentrazione dell'anidride carbonica (CO2) e' " << CO2c << " pg/micron^3" <<  std::endl;

			std::cout << "Concentrazione degli altri nutrienti (pg/micron^3): ";
			std::cin >> Ac;
			if(Ac < 0) 
				Ac = A_ENV;
			std::cout << "La concentrazione degli altri nutrienti e' " << Ac << " pg/micron^3" <<  std::endl;
			
			std::cout << "Concentrazione iniziale dell'acido lattico (pg/micron^3): ";
			std::cin >> AcLc;
			if(AcLc < 0) 
				AcLc = AcL_ENV;
			std::cout << "La concentrazione iniziale dell'acido lattico e' " << AcLc << " pg/micron^3" <<  std::endl;
			
			
			std::cout << "Posizioni degli spigoli del parallelepipedo che definisce il volume (in mm): " << std::endl;
			
			std::cout << "xmin: ";	// X
			std::cin >> xmin;
			xmin *= 1e-3;	// conversione da mm a m
			
			std::cout << "xmax: ";
			std::cin >> xmax;
			xmax *= 1e-3;	// conversione da mm a m
			
			std::cout << "ymin: ";	// Y
			std::cin >> ymin;
			ymin *= 1e-3;	// conversione da mm a m
			
			std::cout << "ymax: ";
			std::cin >> ymax;
			ymax *= 1e-3;	// conversione da mm a m
			
			std::cout << "zmin: ";	// Z
			std::cin >> zmin;
			zmin *= 1e-3;	// conversione da mm a m
			
			std::cout << "zmax: ";
			std::cin >> zmax;
			zmax *= 1e-3;	// conversione da mm a m
			
			v0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
			
			std::cout << "Il volume dell'ambiente e' " << v0*(1e9) << " microlitri" << std::endl;
			
			
			std::cout << "Dose ambientale di radiazione per unita' di tempo (Gy/s): ";
			std::cin >> dose;
			if(dose < 0) 
				dose = DOSERATE_ENV;
			std::cout << "La dose di radiazione nell'ambiente corrisponde a " << dose << "Gy/s" << std::endl;
			
			G = Gc * v0;
			O2 = O2c * v0;
			CO2 = CO2c * v0;
			A = Ac * v0;
			AcL = AcLc * v0;
			pH =  7.5443-AcLc/BufCapEnv;
			std::cout << "Il pH ambientale e' " << pH << std::endl;
			volume0 = v0;
			volume = v0;
			DoseRate = dose;
			break;
			
		default:
			xmin = XMIN_ENV;
			xmax = XMAX_ENV;
			ymin = YMIN_ENV;
			ymax = YMAX_ENV;
			zmin = ZMIN_ENV;
			zmax = ZMAX_ENV;
			
			volume0 = (XMAX_ENV-XMIN_ENV)*(YMAX_ENV-YMIN_ENV)*(ZMAX_ENV-ZMIN_ENV);
			
			T = T_ENV;
			G = G_ENV * volume0;
			O2 = O_ENV * volume0;
			CO2 = CO2_ENV * volume0;
			A = A_ENV * volume0;
			AcL = AcL_ENV * volume0;
			pH =  7.5443-(AcL/volume0)/BufCapEnv;
			volume = volume0;
			DoseRate = DOSERATE_ENV;
		}
	
}

// costruttore che prende i dati da un file (alla fine non da' output su terminale)
Environment::Environment(const std::string filename)
{

	std::ifstream EnvironmentFile( filename.c_str() );
	if( !EnvironmentFile ) 
		{
		std::cout << "Error opening environment file " << filename << std::endl;
		std::exit(-1);
		}
	
	xmin = InputRealPar(EnvironmentFile);
	xmax = InputRealPar(EnvironmentFile);
	ymin = InputRealPar(EnvironmentFile);
	ymax = InputRealPar(EnvironmentFile);
	zmin = InputRealPar(EnvironmentFile);
	zmax = InputRealPar(EnvironmentFile);
	volume0 = (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
	
	T = InputRealPar(EnvironmentFile);
	G = InputRealPar(EnvironmentFile);				// il file contiene la concentrazione di G
	O2 = InputRealPar(EnvironmentFile);				// il file contiene la concentrazione di O2
	CO2 = InputRealPar(EnvironmentFile);			// il file contiene la concentrazione di CO2
	A = InputRealPar(EnvironmentFile);				// il file contiene la concentrazione di A
	AcL = InputRealPar(EnvironmentFile);			// il file contiene la concentrazione di AcL
	pH =  7.5443-AcL/BufCapEnv;						// pH ambientale (si assume la linearita' della relazione
													// che lega la conc. di acido lattico al pH)

	volume = volume0;
	DoseRate = InputRealPar(EnvironmentFile);
	
	G *= volume0;									// dalle concentrazioni alle quantità ambientali ... 
	O2 *= volume0;
	CO2 *= volume0;
	A *= volume0;
	AcL *= volume0;

}

// costruttore copia
Environment::Environment(const Environment& e)
{
	xmin = e.xmin;
	xmax = e.xmax;
	ymin = e.ymin;
	ymax = e.ymax;
	zmin = e.zmin;
	zmax = e.zmax;
	T = e.T;
	G = e.G;
	O2 = e.O2;
	CO2 = e.CO2;
	A = e.A;
	AcL = e.AcL;
	pH =  e.pH;
	volume0 = e.volume0;
	volume = e.volume;
	DoseRate = e.DoseRate;

}

// overloaded =
Environment& Environment::operator=(const Environment& e)
{
	xmin = e.xmin;
	xmax = e.xmax;
	ymin = e.ymin;
	ymax = e.ymax;
	zmin = e.zmin;
	zmax = e.zmax;
	T = e.T;
	G = e.G;
	O2 = e.O2;
	CO2 = e.CO2;
	A = e.A;
	AcL = e.AcL;
	pH =  e.pH;
	volume0 = e.volume0;
	volume = e.volume;
	DoseRate = e.DoseRate;
	
	return *this;
}

// update
void EnvironmentChange(Environment& env, Environment& env_delta)
{
  env.Set_xmin(env.Get_xmin()+env_delta.Get_xmin());
  env.Set_xmax(env.Get_xmax()+env_delta.Get_xmax());
  env.Set_ymin(env.Get_ymin()+env_delta.Get_ymin());
  env.Set_ymax(env.Get_ymax()+env_delta.Get_ymax());
  env.Set_zmin(env.Get_zmin()+env_delta.Get_zmin());
  env.Set_zmax(env.Get_zmax()+env_delta.Get_zmax());
// 	env.xmin += env_delta.xmin;
// 	env.xmax += env_delta.xmax;
// 	env.ymin += env_delta.ymin;
// 	env.ymax += env_delta.ymax;
// 	env.zmin += env_delta.zmin;
// 	env.zmax += env_delta.zmax;
	
  env.SetEnvironmentT(env.GetEnvironmentT()+env_delta.GetEnvironmentT());
  if(env.GetEnvironmentT()<0)
  {
    env.SetEnvironmentT(0);
    std::cout << "Errore in Environment::EnvironmentChange: variazione di temperatura troppo grande e negativa\n" << std::endl;
  }
// 	env.T += env_delta.T;	
// 	if(env.T < 0) 
// 		{
// 		env.T = 0;
// 		std::cout << "Errore in Environment::EnvironmentChange: variazione di temperatura troppo grande e negativa\n" << std::endl;
// 		}
  env.SetEnvironmentG(env.GetEnvironmentG()+env_delta.GetEnvironmentG());
  if(env.GetEnvironmentG()<0)
  {
    env.SetEnvironmentG(0);
    std::cout << "Errore in Environment::EnvironmentChange: variazione di glucosio troppo grande e negativa\n" << std::endl;
  }
// 	env.G += env_delta.G;
// 	if(env.G < 0) 
// 		{
// 		env.G = 0;
// 		std::cout << "Errore in Environment::EnvironmentChange: variazione di glucosio troppo grande e negativa\n" << std::endl;
// 		}
#if 0
	env.O2 += env_delta.O2;
	if(env.O2 < 0) 
		{
		env.O2 = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di ossigeno troppo grande e negativa\n" << std::endl;
		}
	env.CO2 += env_delta.CO2;
	if(env.CO2 < 0) 
		{
		env.CO2 = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di anidride carbonica troppo grande e negativa\n" << std::endl;
		}
	env.A += env_delta.A;
	if(env.A < 0) 
		{
		env.A = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di glutamina troppo grande e negativa\n" << std::endl;
		}
	env.AcL += env_delta.AcL;
	if(env.AcL < 0) 
		{
		env.AcL = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di AcL troppo grande e negativa\n" << std::endl;
		}

	env.volume0 += env_delta.volume0;
	if(env.volume0 < 0) 
		{
		env.volume0 = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di volume0 troppo grande e negativa\n" << std::endl;
		}
	env.volume += env_delta.volume;
	if(env.volume < 0) 
		{
		env.volume = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di volume troppo grande e negativa\n" << std::endl;
		}
	env.DoseRate += env_delta.DoseRate;
	if(env.DoseRate < 0) 
		{
		env.DoseRate = 0;
		std::cout << "Errore in Environment::EnvironmentChange: variazione di dose-rate troppo grande e negativa\n" << std::endl;
		}

	env.pH =  7.5443-(env.AcL/env.volume)/BufCapEnv;
#endif
}


// overloading dell'operatore di output su file
namespace vbl{
std::ostream& operator<<(std::ostream& stream, Environment& cEnvironment)
{

	stream << "T: " <<  cEnvironment.GetEnvironmentT() << " deg C" << std::endl;																// temperatura
	stream << "G: " <<  cEnvironment.GetEnvironmentG()/cEnvironment.GetEnvironmentvolume() << " g/cm^3 -> "  <<  cEnvironment.GetEnvironmentG()*1e-9 << " mg" << std::endl;			// quantita' (mg) di glucosio ambientale
#if 0
	stream << "O2: " <<  cEnvironment.O2/cEnvironment.volume << " g/cm^3 -> "   <<  cEnvironment.O2*1e-9 << " mg" <<  std::endl;		// quantita' (mg) di ossigeno ambientale
	stream << "CO2: " <<  cEnvironment.CO2/cEnvironment.volume << " g/cm^3 -> "   <<  cEnvironment.CO2*1e-9 << " mg" <<  std::endl;	// quantita' (mg) di anidride carbonica ambientale
	stream << "A: " <<  cEnvironment.A/cEnvironment.volume << " g/cm^3 -> "   <<  cEnvironment.A*1e-9 << " mg" <<  std::endl;		// quantita' (mg) di altri nutrienti nell'ambiente
	stream << "AcL: " <<  cEnvironment.AcL/cEnvironment.volume << " g/cm^3 -> "   <<  cEnvironment.AcL*1e-9 << " mg" <<  std::endl;	// quantita' (mg) di acido lattico nell'ambiente
	stream << "pH: " <<  cEnvironment.pH << std::endl;																			// pH ambientale
	stream << "xmin: " << 1e-3*cEnvironment.xmin << " mm" << std::endl;															// limiti spaziali dell'ambiente			
	stream << "xmax: " << 1e-3*cEnvironment.xmax << " mm" << std::endl;															
	stream << "ymin: " << 1e-3*cEnvironment.ymin << " mm" << std::endl;															
	stream << "ymax: " << 1e-3*cEnvironment.ymax << " mm" << std::endl;															
	stream << "zmin: " << 1e-3*cEnvironment.zmin << " mm" << std::endl;															
	stream << "zmax: " << 1e-3*cEnvironment.zmax << " mm" << std::endl;															
	stream << "volume0: " <<  1e-9*cEnvironment.volume0 << " µl" << std::endl;													// volume iniziale dell'ambiente
	stream << "volume: " <<  1e-9*cEnvironment.volume << " µl" << std::endl;														// volume dell'ambiente
	stream << "Radiazione di fondo: " <<  cEnvironment.DoseRate << " Gy/s" << std::endl;										// rate di radiazione nell'ambiente (Gy/s)
#endif
	return stream;
}
}//namespace vbl{

void Environment::PrintEnvironmentData(std::ofstream& stream, long int nrec)
{

	static bool first_print=true;

// stampa degli header nel caso che questa sia la prima volta che si stampa
	if(first_print)
		{
		first_print = false;
		
		stream << "n \t" \
			<< "T \t" \
			<< "G \t" \
			<< "O2 \t" \
			<< "CO2 \t" \
			<< "A \t" \
			<< "AcL \t" \
			<< "pH \t" \
			<< "xmin \t" \
			<< "xmax \t" \
			<< "ymin \t" \
			<< "ymax \t" \
			<< "zmin \t" \
			<< "zmax \t" \
			<< "volume0 \t" \
			<< "volume \t" \
			<< "DoseRate "<< std::endl;
		}

	stream << std::setprecision(8) << std::scientific \
			<< nrec << "\t" \
			<< T << "\t" \
			<< G << "\t" \
			<< O2 << "\t" \
			<< CO2 << "\t" \
			<< A << "\t" \
			<< AcL << "\t" \
			<< pH << "\t" \
			<< xmin << "\t" \
			<< xmax << "\t" \
			<< ymin << "\t" \
			<< ymax << "\t" \
			<< zmin << "\t" \
			<< zmax << "\t" \
			<< volume0 << "\t" \
			<< volume << "\t" \
			<< DoseRate << std::endl;


}

void vbl::Environment::assign(const boost::property_tree::ptree& pt)
{
#define DOPT(name_buffer)name_buffer = pt.get<double>(#name_buffer)
  DOPT(T);
  DOPT(G);
  DOPT(O2);
  DOPT(CO2);
  DOPT(A);
  DOPT(AcL);
  DOPT(pH);
  DOPT(xmin);
  DOPT(xmax);
  DOPT(ymin);
  DOPT(ymax);
  DOPT(zmin);
  DOPT(zmax);
  DOPT(volume0);
  DOPT(volume);
  DOPT(DoseRate);
#undef DOPT
}

boost::property_tree::ptree vbl::Environment::as_ptree() const
{
  boost::property_tree::ptree pt;
#define DOPT(name_buffer) pt.put(#name_buffer, name_buffer)
  
  //for development 
//   for(auto entry:vector_of_env_parameters)
//   {
//     std::cout << "DOPT(" << entry << ");"<< std::endl;
//     //DOPT(entry);
//   }
  
  DOPT(T);
  DOPT(G);
  DOPT(O2);
  DOPT(CO2);
  DOPT(A);
  DOPT(AcL);
  DOPT(pH);
  DOPT(xmin);
  DOPT(xmax);
  DOPT(ymin);
  DOPT(ymax);
  DOPT(zmin);
  DOPT(zmax);
  DOPT(volume0);
  DOPT(volume);
  DOPT(DoseRate);
#undef DOPT
  return pt;
}

void Environment::WriteEnvironment( std::ofstream& stream )
{
	stream.write( (char*)(&T), sizeof(double) );
	stream.write( (char*)(&G), sizeof(double) );
	stream.write( (char*)(&O2), sizeof(double) );
	stream.write( (char*)(&CO2), sizeof(double) );
	stream.write( (char*)(&A), sizeof(double) );
	stream.write( (char*)(&AcL), sizeof(double) );
	stream.write( (char*)(&pH), sizeof(double) );
	stream.write( (char*)(&xmin), sizeof(double) );
	stream.write( (char*)(&xmax), sizeof(double) );
	stream.write( (char*)(&ymin), sizeof(double) );
	stream.write( (char*)(&ymax), sizeof(double) );
	stream.write( (char*)(&zmin), sizeof(double) );
	stream.write( (char*)(&zmax), sizeof(double) );
	stream.write( (char*)(&volume0), sizeof(double) );
	stream.write( (char*)(&volume), sizeof(double) );
	stream.write( (char*)(&DoseRate), sizeof(double) );

}

void Environment::ReadEnvironment( std::ifstream& stream )
{
	stream.read( (char*)(&T), sizeof(double) );
	stream.read( (char*)(&G), sizeof(double) );
	stream.read( (char*)(&O2), sizeof(double) );
	stream.read( (char*)(&CO2), sizeof(double) );
	stream.read( (char*)(&A), sizeof(double) );
	stream.read( (char*)(&AcL), sizeof(double) );
	stream.read( (char*)(&pH), sizeof(double) );
	stream.read( (char*)(&xmin), sizeof(double) );
	stream.read( (char*)(&xmax), sizeof(double) );
	stream.read( (char*)(&ymin), sizeof(double) );
	stream.read( (char*)(&ymax), sizeof(double) );
	stream.read( (char*)(&zmin), sizeof(double) );
	stream.read( (char*)(&zmax), sizeof(double) );
	stream.read( (char*)(&volume0), sizeof(double) );
	stream.read( (char*)(&volume), sizeof(double) );
	stream.read( (char*)(&DoseRate), sizeof(double) );

}
