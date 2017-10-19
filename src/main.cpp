//
// includes: attenzione l'ordine degli includes e' importante ... 
// 


#include <../include/vbl.h>
using namespace std;
using namespace vbl;


// **********************************************************************************
// main program 
// 
// Edoardo Milotti - aprile 2009 - febbraio 2011
// 
// **********************************************************************************
//

int main( int argc, char* argv[] )
{
  int run_type = 0;	// tipo di run 
			// 0 = run iniziale, input da terminale
			// 1 = run iniziale, input da command file
			// 2 = continuazione di un run precedente
  bool terminal = false; // nel caso sia un run iniziale questa var. booleana dice se l'input iniziale e' da terminale
  
  int idum=-1; 	// seme iniziale per ran2 (Numerical Recipes)
  
  string run_name;
  
  CellsSystem CellsSystem;	// Standard allocation of the CellsSystem (in this case, the initial dynamic reserve is 2000000)
  //CellsSystem.Set_BV_reserve(RESERVE_BV);
      
  if( argc < 3 )
  {
    cout << "\nparametri insufficienti: \nUsage: \ninput da terminale: /.Sim3D <random seed (intero neg.)> \n";
    cout << "in general: /.Sim3D <t/c> <seed> or ./Sim3D r <run number>\n" << endl;
    exit(-1);		// non ci sono abbastanza parametri
  }
  
  if( argv[1][0] == 't' || argv[1][0] == 'T')	// terminale? 
  {
    terminal = true;
    idum = atoi(argv[2]);		// The seed is the third parameter
    if(idum >= 0) idum = -1;		// Currently using the ran2 generator of Numerical Recipes, so the seed must be negative
    run_type = 0;
  }
  else if( argv[1][0] == 'c' || argv[1][0] == 'C' )	// command file?
  {
    terminal = false;
    idum = atoi(argv[2]);				// The third parameter is the random number seed
    if(idum > 0) idum = -1;
    run_type = 1;
    if( strlen(argv[1]) > 1 )
	  CellsSystem.Set_Commands( argv[1] );		// In this case, the command file is specified
    else
	  CellsSystem.Set_Commands( "commands.txt" );	// In this case, it gets its default name
  }
  else if( argv[1][0] == 'r' || argv[1][0] == 'R' )	// Continuation of a previous run?
  {
    terminal = false;
    run_type = 2;
    run_name = argv[2];				// The third parameter is the running name that continues
  }
  else 
  {
    cout << "\ncomando non valido: \nUsage: \ninput da terminale: /.Sim3D <random seed (intero neg.)> \n";
    cout << "in generale: /.Sim3D <t/c> <seed> oppure ./Sim3D r <run number>\n" << endl;
    exit(-1);							// comando non valido
  }
      
      
  cout << endl;
  
  if(run_type == 0)
	  cout << "*** terminal input ***" << endl;
  else if (run_type == 1)
	  {
	  cout << "*** command file input ***" << endl;
	  cout << "command file: " << CellsSystem.Get_Commands() << endl;
	  }
  else if (run_type == 2)
	  cout << "*** continue run " << run_name << " ***" << endl;
	      
      
      
      
  if(run_type == 0 || run_type == 1)// In case the run is not the continuation of a previous run, some standard operations are performed
  {
    CellsSystem.Set_idum( idum );// Seed of random number generator

    if(argc < 4)
      CellsSystem.Set_CellTypeFile( "CellType.txt" );
    else// If there are at least three arguments then the third is the name of the CellType file (which is in the same dir)
      CellsSystem.Set_CellTypeFile( argv[3] );

    if(argc < 5)
      CellsSystem.Set_EnvironmentFile( "Environment.txt" );
    else// If there are at least 4 arguments then the fourth is the name of the Environment file (which is in the same executable)
      CellsSystem.Set_EnvironmentFile( argv[4] );

    if(argc < 6)
      // If there are less than 5 arguments then the alternate type is the same as the initial type
      if(argc <4 ) 
	CellsSystem.Set_CellTypeFileAlt( "CellType.txt" );
      else
	CellsSystem.Set_CellTypeFileAlt( argv[3] );
    else// If there are at least 5 arguments then the fifth is the name of the file CellTypeAlt (which is in the same dir of the executable)
	CellsSystem.Set_CellTypeFileAlt( argv[5] );
	  
	  
    cout << "\n*** VBL - Virtual Biophysics Lab simulation program ***";
    cout << "\n\nInitialized CellsSystem." << endl;
    cout << "Cell system is of size " << sizeof(CellsSystem) << " bytes" << endl;
    cout << "CellType is taken from file " << CellsSystem.Get_CellTypeFile( ) << endl;
    cout << "CellTypeAlt is taken from file " << CellsSystem.Get_CellTypeFileAlt( ) << endl;
    cout << "The environment is taken from file " << CellsSystem.Get_EnvironmentFile( ) << endl;

    CellsSystem.InitializeCellsSystem( terminal );
    cout << "Initialization completed" << endl;	
    CellsSystem.RunDefinition( );// Run number and output directory output directory & output file opening for metabolism
    CellsSystem.Set_nconfiguration( 0 ); // The configuration number is initialized to 0
  }
  else if( run_type == 2 )   
  {
    cout << "\n*** VBL - Virtual Biophysics Lab simulation program ***";
    cout << "\n\nInitialized CellsSystem" << endl;
    cout << "Cell system is of size " << sizeof(CellsSystem) << " bytes" << endl;
    cout << "Continuing ...  " << run_name << endl;
    
    CellsSystem.RunDefinition( run_name );//Setup in case of continuation of a run
    cout << "filenames definiti ...\n" << endl;;
    CellsSystem.ReadCellsSystem( );//Reading the saved configuration
    cout << "configurazione letta ... \n" << endl;
    CellsSystem.Set_ready2start( true );// There is no initialization in this case ...
    cout << "ready2start = true ... \n" << endl;
  }
      
  cout << "\nComplete file definitions" << endl;
      
      
  if( CellsSystem.Get_sim_type() == 1 )	// If this is a 3D simulation then an initial calculation of geometry is made
  {
    cout << "\nInitial calculation of geometry" << endl;
    CellsSystem.Geometry( );// Initial calculation of cluster geometry
    cout << "Initial calculation of complete geometry" << endl;
    CellsSystem.Set_time_from_CGAL(0.);	// Timer reset from last call to CGAL
  }
      
  cout << "\nTimer setup and statistics" << endl;
  
  // Control prints of the first, second, and last cells on the logfile
  if(run_type == 0 || run_type == 1)
    CellsSystem.Print2logfile("Cell status at the end of initialization");
  else if (run_type == 2)
    CellsSystem.Print2logfile("Cell status at restart of simulation");
  

  CellsSystem.CPU_timer(timer_button::Start_timer);		// start del CPU timer (e reset del timer degli intertempi)
  CellsSystem.Timing( true );				// reset del timer
  CellsSystem.StepStat( true );			// reset delle statistiche (azzera anche il vettore convergence_fail)

  cout << "\nStartup completed" << endl;
      
  // ********************** main simulation loop **************************

  unsigned int returnValue = CellsSystem.runMainLoop(boost::none);

  // **************************** fine loop principale ******************************


  cout << "\nFinishing run" << endl;
  
  CellsSystem.Print2logfile("Cells at the end of the run");
  CellsSystem.WriteCellsSystem( );					// dump of the final configuration
  
  cout << "Final configuration written on file" << endl;
  
  CellsSystem.CloseOutputFiles();						// Closing output files
      

  return 0;
      
}
