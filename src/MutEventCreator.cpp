/*
 * MutEventCreator.cpp
 *
 *  Created on: 11/giu/2018
 *      Author: sabry
 */


#include "MutEventCreator.h"

namespace vbl 
{

vbl::MutEventCreator::MutEventCreator() 
{
	m_eventTime = 0;
	m_pAlt = 0;
	m_pHThreshold = 0;
	m_NEvent = 0;
	m_MutModality = NoMutation;
}

// vbl::MutEventCreator::MutEventCreator(const vbl::MutEventCreator& mutev)
// {
//   m_eventTime = mutev.m_eventTime;
//   m_pAlt = mutev.m_pAlt;
//   m_pHThreshold = mutev.m_pHThreshold;
//   m_NEvent = mutev.m_NEvent;
//   m_MutModality = mutev.m_MutModality;
// }
vbl::MutEventCreator::~MutEventCreator()
{
#ifdef ENABLE_MUTATION_DEBUG_OUTPUT
  cout << "MutEventCreator destructor called" << endl;
#endif
}
// // overloaded =
// MutEventCreator& MutEventCreator::operator=(const MutEventCreator& mutev)
// {
// 	m_eventTime = mutev.m_eventTime;
//   m_pAlt = mutev.m_pAlt;
//   m_pHThreshold = mutev.m_pHThreshold;
//   m_NEvent = mutev.m_NEvent;
//   m_MutModality = mutev.m_MutModality;
// 	
// 	return *this;
// }
// from sabry
// MutEventCreator::~MutEventCreator() 
// {
// 	// TODO Auto-generated destructor stub
//   cout << "MutEventCreator destructor called" << endl;
// }



bool MutEventCreator::Generate_SingleMutEvent(double treal)
{
	bool temp;
  double myRandomNumber = ran2();
#ifdef ENABLE_MUTATION_DEBUG_OUTPUT
  cout << "myRandomNumber: " << myRandomNumber << endl;
  cout << "treal: " << treal << endl;
  cout << "m_eventTime: " << m_eventTime << endl;
  cout << "m_pAlt: " << m_pAlt << endl;
#endif
	if(treal < m_eventTime || (myRandomNumber > m_pAlt || m_pAlt < 2 ) )
  {
		temp = false;
	}

	else
  {
		temp = true;
		cout<< endl << "****************mutation occurred!" << endl;
#ifdef ENABLE_MUTATION_DEBUG_OUTPUT
		cout <<  "m_eventTime: " << m_eventTime << "  treal: " << treal << endl;
#endif
		if(m_pAlt == 2)
      m_pAlt = 0;
	}
	
  return temp;
}

bool MutEventCreator::Generate_RateMutPHThreshold(double pHExtra)
{
  bool temp;
  double myRandomNumber = ran2();
  if( pHExtra < m_pHThreshold && myRandomNumber < m_pAlt)
  {
    temp = true;
    m_NEvent ++;
#ifdef ENABLE_MUTATION_DEBUG_OUTPUT
    cout << "pH_ex: " << pHExtra << "<" << m_pHThreshold << " and rand2(): " << myRandomNumber << "<" << m_pAlt << endl;
    cout << "m_NEvent: " << m_NEvent << std::endl;
#endif
  }
  else 
  {
    temp = false;
  }
  return temp;
}

bool MutEventCreator::Generate_SingleMutPHThreshold(double pHExtra)
{
	bool temp;

	if(pHExtra > m_pHThreshold || m_NEvent == SINGLE_MUTATION) 
  {  //no mutation condition
		temp = false;
	}
	else 
  {
		temp = true;
		cout<< endl << "****************mutation occurred!" << endl;
		cout<< "pH extracellular:" << pHExtra << endl;
		m_NEvent ++;
	}

	return temp;
}



void MutEventCreator::InitMutEvent(bool terminal, std::ifstream& commands){

	bool eventON;
	int eventType;

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
      cout<<"Select:"<< endl <<"1: Single mutation event at fixed time" << endl 
      <<"2: Single mutation driven by PH "<<endl
      <<"3: Multiple mutation driven by PH with rate m_pAlt "<<endl;
      cin >> eventType;

      switch(eventType)
      {
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
        
        case(3):
          m_MutModality = MutPHRate;
          cout << "Rate of mutation from ph Value:";
          cin >> m_pHThreshold;
          cout << "Rate probability (between 0.0 and 1.0)";
          cin >> m_pAlt;
        break;
      }
    }
    else
    {
      //command file case
      commands >> eventType;
      switch(eventType)
      {
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
        
        case(3):
          m_MutModality = MutPHRate;
          cout << "Rate of mutation from ph Value: ";
          commands >> m_pHThreshold;
          cout << m_pHThreshold << endl;
          cout << "Rate probability (between 0.0 and 1.0): ";
          commands >> m_pAlt;
          cout << m_pAlt << endl;
        break;
      }
    }//close else (terminal)
  } //close if (eventON)
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

void MutEventCreator::assign(const boost::property_tree::ptree& pt)
{
  /** T.F. why the hell temp parameters?
   */
  //double m_eventTime; //tempo associato ad un singolo evento speciale
	//double m_pAlt;       // probabilità dell'evento speciale "transizione al tipo alternativo"
	//double m_pHThreshold;

	//int m_NEvent;  //number of induced mutation

	//Flags to set the result of I/O routine
  //MutModality m_MutModality;
#define DOPT(name_buffer)name_buffer = pt.get<double>(#name_buffer)
  DOPT(m_eventTime);
  DOPT(m_pAlt);
  DOPT(m_pHThreshold);
#undef DOPT
  
#define DOPT(name_buffer)name_buffer = pt.get<int>(#name_buffer)
  DOPT(m_NEvent);
  int buffer = pt.get<int>("m_MutModality");
  switch(buffer)
  {
    case NoMutation:
      m_MutModality = NoMutation;
      break;
    case SingleMutEvent:
      m_MutModality =SingleMutEvent;
      break;
    case MutPHThreshold:
      m_MutModality = MutPHThreshold;
      break;
    case MutPHRate:
      m_MutModality = MutPHRate;
      break;
    default:
      cout << "unknown mutoation option " << buffer << endl;
  }
#undef DOPT
}

boost::property_tree::ptree MutEventCreator::as_ptree() const
{
  boost::property_tree::ptree pt;
  
#define DOPT(name_buffer) pt.put(#name_buffer, name_buffer)
  DOPT(m_eventTime);
  DOPT(m_pAlt);
  DOPT(m_pHThreshold);
  DOPT(m_NEvent);
  DOPT(m_MutModality);
#undef DOPT
  
  return pt;
}


}///end namespace vbl
