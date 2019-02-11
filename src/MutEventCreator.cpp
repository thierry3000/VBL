/*
 * MutEventCreator.cpp
 *
 *  Created on: 11/giu/2018
 *      Author: sabry
 */


#include "MutEventCreator.h"

MutEventCreator::MutEventCreator() {
	m_eventTime = 0;
	m_pAlt = 0;
	m_pHThreshold = 0;
	m_NEvent = 0;
	m_MutModality = NoMutation;
}

MutEventCreator::~MutEventCreator() {
	// TODO Auto-generated destructor stub
}



bool MutEventCreator::Generate_SingleMutEvent(double treal)
{
	bool temp;
  auto myRandomNumber = ran2();
  cout << "debug random: " << myRandomNumber << endl;
	if(treal < m_eventTime || (ran2() > m_pAlt || m_pAlt < 2 ) ){
		temp = false;
	}

	else
  {
		temp = true;
		cout<< endl << "****************mutation occurred!" << endl;
		cout <<  "m_eventTime: " << m_eventTime << "  treal: " << treal << endl;
		if(m_pAlt == 2) m_pAlt = 0;
	}
	 return temp;

}


bool MutEventCreator::Generate_SingleMutPHThreshold(double pHExtra){

	bool temp;

	if(pHExtra > m_pHThreshold || m_NEvent == SINGLE_MUTATION) {  //no mutation condition
		temp = false;
	}
	else {
		temp = true;
		cout<< endl << "****************mutation occurred!" << endl;
		cout<< "pH extracellular:" << pHExtra << endl;
		m_NEvent ++;
	}



	return temp;
}



void MutEventCreator::InitMutEvent(bool terminal, std::ifstream& commands){

	bool eventON;
	int temp;

	if( terminal )
	        {
	        cout << "\nIs there a mutation event for this run? (NO = 0, YES = 1) ";
	        cin >> eventON;
	        }
	    else
	        {
	        commands >> eventON;
	        }

	    if( eventON )
	        {
	        cout << "There is a mutation event for this run:" << endl;
	        if( terminal )
	            {
	            cout<<"Select:"<< endl <<"1: Single mutation event at fixed time" << endl <<"2: Single mutation driven by PH "<<endl;
	            cin >> temp;

	            switch(temp){

	            	case(1):
	            		m_MutModality = SingleMutEvent;
	                 	cout << "Mutation event time (s) (real time): ";
			        	cin >> m_eventTime;
			            cout << "Mutation event probability ( select 2 = for a single mutation event): ";
			            cin >> m_pAlt;

			        break;

	            	case(2):
						m_MutModality = MutPHThreshold;
						cout << "Set the PH extracellular value:";
		                cin >> m_pHThreshold;
					break;

	            }

	            }
	        else
	            {
	            commands >> temp;

	            switch(temp){

	            case(1):
					  m_MutModality = SingleMutEvent;
	                  cout<< "Event selected:" << m_MutModality;
		              commands >> m_eventTime;
	            	  cout << "\nMutation event time: " << m_eventTime << endl;

	            	  commands >> m_pAlt;
	            	  cout << "\nMutation event probability (2 = single mutation event): " << m_pAlt << endl;

	            break;

	            case(2):
						m_MutModality = MutPHThreshold;
			            cout<< "Event selected:" << m_MutModality;
	                    commands >> m_pHThreshold;
	                    cout << "\n The PH extracellular value for the mutation event is" << m_pHThreshold;
	            break;


	            }
	            }
	        } //close if(eventON)

	    else
	        {
	        m_MutModality = NoMutation; //no mutation event
	        }


  return;
}

//scrive su output file: run_data.txt

void MutEventCreator::PrintMutEvent(std::ofstream& run_data){


	switch(m_MutModality)
	{

		case(NoMutation):
	        run_data<< "There is no special event for this run." << endl;

		break;

		case(SingleMutEvent):
		    run_data<< "There is a mutation event for this run" << endl;
			run_data<< "Event selected: SingleMutEvent" << endl;
			run_data << "Mutation event time: " << m_eventTime << endl;
			run_data << "Mutation event probability (2 = single mutation event): " << m_pAlt << endl;
		break;

		case(MutPHThreshold):
			run_data<< "There is a mutation event for this run" << endl;
			run_data<< "Event selected: MutPHThreshold" << endl;
			run_data<< "The PH extracellular value for the mutation event is" << m_pHThreshold <<endl;

			break;
	}

	return;
}


// Write binary for run continuation
void MutEventCreator::WriteMutEvent(ofstream& stream)
{

	stream.write( (char*)(&m_eventTime), sizeof(double) );
	stream.write( (char*)(&m_pAlt), sizeof(double) );
	stream.write( (char*)(&m_pHThreshold), sizeof(double) );
	stream.write( (char*)(&m_NEvent), sizeof(int) );
	stream.write( (char*)(&m_MutModality), sizeof(MutModality) );

      return;
}


// Read binary for  run continuation
void MutEventCreator::ReadMutEvent(ifstream& stream)
{
	stream.read( (char*)(&m_eventTime), sizeof(double) );
	stream.read( (char*)(&m_pAlt), sizeof(double) );
	stream.read( (char*)(&m_pHThreshold), sizeof(double) );
	stream.read( (char*)(&m_NEvent), sizeof(int) );
	stream.read( (char*)(&m_MutModality), sizeof(MutModality) );

	return;

}




