// MfiaHysteresis.cpp -- Hysteresis simulations
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//

#include"MC_Hysteresis.hpp"
#include"MC_Walker.hpp"
#include"Meas_Trace.hpp"
#include"Mfia_Hamiltonian.hpp"
#include"Random.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>


void DoIt(int& argc, char* argv[])
{
   int iproc = 0;
   int nproc = 1;
#  if USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   ProgramOptions options("MfiaHysteresis","Hysteresis simulations for the MFIA model.");

   // Construct a hamiltonian object 
   typedef Mfia_Hamiltonian HAMILTONIAN;
   HAMILTONIAN hamilton;
   hamilton.D     = 2;
   hamilton.L     = 32;
   hamilton.JSIGN = +1;
   hamilton.A     = 0;
   options.add_option("dim",    "dimension of spin lattice",     'D', &(hamilton.D) );
   options.add_option("length", "extent of lattice in each dim", 'L', &(hamilton.L) );
   options.add_option("jsign",  "sign in front of exchange",     'J', &(hamilton.JSIGN) );
   options.add_option("A",      "magnitude of mean field",       'A', &(hamilton.A) );

   // Construct the simulation object
   MC_Hysteresis simulator;
   simulator.NStep0 = 1000;          // Number of measurements to take
   simulator.NStep1 = 1000;          // Number of measurement blocks before begin averaging
   simulator.NStep  = 5;             // Number of steps to take between measurements to minimize correlations
   simulator.NCycle = 100;           // Number of trials or hysteresis loops
   simulator.kT     = 1.0;           // Temperature of simulation
   simulator.H0     = 3.0;
   simulator.H1     =-3.0;
   options.add_option("nstep0",  "time of equil or first part",         ' ', &(simulator.NStep0) ); 
   options.add_option("nstep1",  "time of second part",                 ' ', &(simulator.NStep1) ); 
   options.add_option("nstep",   "time between measurements",           ' ', &(simulator.NStep) ); 
   options.add_option("ncycle",  "number of cycles or trials",          ' ', &(simulator.NCycle) ); 
   options.add_option("numwalk", "number of independent random walks",  ' ', &(simulator.NWalkPerProcess) ); 
   options.add_option("kT",      "temperature of the simulation",       ' ', &(simulator.kT) ); 
   options.add_option("H0",      "first field value of the simulation", ' ', &(simulator.H0) ); 
   options.add_option("H1",      "second field value of the simulation",' ', &(simulator.H1) ); 

   int waveform = 0;
   options.add_option("waveform","0=quench,1=sweep,2=sine,3=square,4=sawtooth", ' ', &waveform);

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
   bool verbose = options.get_value<bool>("verbose") && (iproc==0);

   switch(waveform)
   {
      case 0: simulator.waveform = MC_Hysteresis::quench;   break;
      case 1: simulator.waveform = MC_Hysteresis::sweep;    break;
      case 2: simulator.waveform = MC_Hysteresis::sine;     break;
      case 3: simulator.waveform = MC_Hysteresis::square;   break;
      case 4: simulator.waveform = MC_Hysteresis::sawtooth; break;
   }
   hamilton.H = simulator.H0;
   hamilton.init(verbose);

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   typedef MC_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   hamilton.initial_mixed(walkerpool[0].sigma);
   // Flip +/- FM depending on applied field
   if( hamilton.H<0 )
      for(int ispin=0; ispin<walkerpool[0].sigma.size(); ispin++) walkerpool[0].sigma[ispin] *= -1;
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   simulator.mp = MPI_Struct::world();
   simulator.InitPool(walkerpool);
   // Flip the configuration of half the walks so evenly sample +/-FM configurations (for low fields)
   if( walkerpool.size()>1 && std::fabs(hamilton.H)<0.1 )
   {
      for(int ispin=0; ispin<walkerpool[1].sigma.size(); ispin++) walkerpool[1].sigma[ispin] *= -1;
      hamilton.calc_observable(walkerpool[1].sigma,walkerpool[1].now);
      for(int iwalk=3; iwalk<walkerpool.size(); iwalk+=2)
         walkerpool[iwalk].copy(walkerpool[1]);
   }

   // seed the random number geneator
   simulator.urng.seed( ParallelSeed(SeedFromClock()) );

   // Object for accumulating equilibrium statistics
   Meas_Trace measure_obj;
   measure_obj.mp = simulator.mp;
   if( simulator.waveform==MC_Hysteresis::quench || simulator.waveform==MC_Hysteresis::sweep ) measure_obj.cyclic = false;
   long NStepCycle = (simulator.NStep0+simulator.NStep1);
   long NMeasCycle = NStepCycle/simulator.NStep;
   if( (NStepCycle-NMeasCycle*simulator.NStep)>0 ) NMeasCycle++;
   measure_obj.init(NStepCycle,NMeasCycle);

   // Do the simulation
   simulator.DoSample(hamilton,walkerpool,measure_obj); 
   options.write();
   measure_obj.write();

}


int main(int argc, char* argv[])
{

   MPI_Init(&argc,&argv);
   int iproc=0;
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
    
   const bool very_verbose = false;
   time_t start_time_main;
   if( iproc==0 )
   {
      if( very_verbose )
      {
         std::cout << "cmdline:";
         for(int i=0; i<argc; i++) std::cout << " " << argv[i];
         std::cout << std::endl;
      }
      time(&start_time_main);
      struct tm * timeinfo = localtime(&start_time_main);
      std::cout << "start time = " << asctime(timeinfo); // << std::endl;
      int pid = 0;
      if(true) pid = getpid();
      std::cout << "process id = " << pid << std::endl;
   }

   DoIt(argc,argv);

   MPI_Finalize();

}
