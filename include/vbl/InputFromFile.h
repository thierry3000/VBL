//
// funzioni di input da files che hanno il formato standard
// le linee nel file hanno tutte la forma 
// name(string) \t value(number) \t description(string) 
//
// EM 19/1/2008
//
// **********************************************************************************

#ifndef INPUT_FROM_FILE_H
#define INPUT_FROM_FILE_H  // header guard

#include <fstream>

inline int InputIntPar(std::ifstream& ParameterFile)
{
	int ivalue;
	std::string name, description;

	getline(ParameterFile,name,'\t');
	ParameterFile >> ivalue;
	getline(ParameterFile,description,'\n');	
	
	return ivalue;
}

inline double InputRealPar(std::ifstream& ParameterFile)
{
	double value;
	std::string name, description;

	getline(ParameterFile,name,'\t');
	ParameterFile >> value;
	getline(ParameterFile,description,'\n');	
	
	return value;
}

#endif //#ifndef INPUT_FROM_FILE_H
