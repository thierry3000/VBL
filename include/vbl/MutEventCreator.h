/*
 * MutEventCreator.h
 *
 *  Created on: 11/giu/2018
 *      Author: sabry
 */

#ifndef MUTEVENTCREATOR_H_
#define MUTEVENTCREATOR_H_

#include <iostream>			// inclusione della libreria standard di I/O
#include <iomanip>			// inclusione dei manipolatori per l'I/O
#include <fstream>			// dichiarazione headers che servono per l'I/O su file

#include "Utilities.h"

//class CellsSystem;

using namespace std;

//list of the possible mutation modalily
enum MutModality{NoMutation = 0, SingleMutEvent = 1 , MutPHThreshold = 2} ;

const int SINGLE_MUTATION = 1; //numero massimo di mutazioni eseguibili

class MutEventCreator {


	//friend class CellsSystem;
private:
	double m_eventTime; //tempo associato ad un singolo evento speciale
	double m_pAlt;       // probabilità dell'evento speciale "transizione al tipo alternativo"
	double m_pHThreshold;

	int m_NEvent;  //number of induced mutation

	//Flags to set the result of I/O routine
    MutModality m_MutModality;


public:
	MutEventCreator();


	virtual ~MutEventCreator();

	//Getter
	double Get_eventTime() {return m_eventTime;}
	double Get_pAlt() {return m_pAlt;}
	double Get_pHThreshold() {return m_pHThreshold;}
	MutModality Get_MutModality(){return m_MutModality;}

	//Setter
	void Set_eventTime(double eventTime){ m_eventTime = eventTime; return;}
	void Set_pAlt(double pAlt){m_pAlt = pAlt; return;}
	void Set_pHThreshold(double pHThreshold){m_pHThreshold = pHThreshold; return;}


	bool Generate_SingleMutEvent(double treal);
  bool Generate_SingleMutPHThreshold(double pHExtra);


    //I/O routine (to be modified if a new function Generator is included)
    void InitMutEvent(bool terminal, std::ifstream& commands);
    void PrintMutEvent(std::ofstream& run_data);  // write mutation option on run_data.txt
    void WriteMutEvent(ofstream& stream);
    void ReadMutEvent (ifstream& stream);
};

#endif /* MUTEVENTCREATOR_H_ */
