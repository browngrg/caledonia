// MfiaMetrop.cpp -- Drive Replica Exchange Monte Carlo
//
//
// History:
//
// June 29, 2016 Convert to IsingCaledonia for simple Ising simulation
//
// July 17, 2015 Convert MfiaMetrop1.cpp to Replica Exchange driver
//
// Sept 29, 2014 Convert Latgas1.cpp to MfiaMetrop.cpp for Metropolic sampling
//               of Histograms
//
// June 16, 2014 Make Latgas_Hamiltonian a model of the Hamiltonian concept
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// Wham3.cpp -- Generate the trajectory of one Metropolis Monte Carlo walker
//
// June 4, 2014    Refactored for larger MC scheme
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//

#include"MC_Metropolis.hpp"
#include"RE_Walker.hpp"
#include"Ising_Measure.hpp"
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

   ProgramOptions options("MfiaMetrop","Metropolis sampling for the MFIA model.");

   // Construct a hamiltonian object 
   typedef Mfia_Hamiltonian HAMILTONIAN;
   HAMILTONIAN hamilton;
   hamilton.D     = 2;
   hamilton.L     = 16;
   hamilton.JSIGN = +1;
   hamilton.H     = 0;
   hamilton.A     = 0;
   options.add_option("dim",    "dimension of spin lattice",     'D', &(hamilton.D) );
   options.add_option("length", "extent of lattice in each dim", 'L', &(hamilton.L) );
   options.add_option("jsign",  "sign in front of exchange",     'J', &(hamilton.JSIGN) );
   options.add_option("H",      "magnitude of magnetic field",   'H', &(hamilton.H) );
   options.add_option("A",      "magnitude of mean field",       'A', &(hamilton.A) );

   // Construct the simulation object
   MC_Metropolis simulator;
   simulator.write_mcssec = true;
   simulator.NMeas  = 10000;         // Number of measurements to take
   simulator.NTherm = 100;           // Number of measurement blocks before begin averaging
   simulator.NStep  = 10;            // Number of steps to take between measurements to minimize correlations
   simulator.kTmin  = 1.3;           // Temperature of simulation
   simulator.kTmax  = 3.3;           // Temperature of simulation
   options.add_option("nmeas",  "number of measurements",         ' ', &(simulator.NMeas) ); 
   options.add_option("ntherm", "number of thermalize moves",     ' ', &(simulator.NTherm) ); 
   options.add_option("nstep",  "number of MC moves between meas",' ', &(simulator.NStep) ); 
   options.add_option("kTmin",  "temperature of simulation",      ' ', &(simulator.kTmin) );
   options.add_option("kTmax",  "temperature of simulation",      ' ', &(simulator.kTmax) );
   options.add_option("nwalk",  "number of walkers per process",  ' ', &(simulator.NWalkPerProcess) );
   options.add_option("nexchg", "number of meas between replica exchange", ' ', &(simulator.NExchg) );

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
   bool verbose = options.get_value<bool>("verbose") && (iproc==0);

   hamilton.init(verbose);

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   typedef RE_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   hamilton.initial_mixed(walkerpool[0].sigma);
   // Flip +/- FM depending on applied field
   const bool prepare_metastable = false;
   if( hamilton.H<0 || prepare_metastable )
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
      {
         walkerpool[iwalk].sigma = walkerpool[1].sigma;
         walkerpool[iwalk].now = walkerpool[1].now;
      }
   }
   // Read configurations from disk, if they exist
   simulator.read_config(hamilton,walkerpool);
   // seed the random number geneator
   simulator.urng.seed( ParallelSeed(SeedFromClock()) );
#  if 0
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(simulator.urng,rng_output,std::cout);
   if( rng_failed )
   {
      std::cout << "Problem detected with random number generator" << std::endl;
      return;
   }
#  endif

   // Object for accumulating equilibrium statistics
   Ising_Measure measure_obj;
   measure_obj.add_header( hamilton.header() ); 
   measure_obj.add_header( simulator.header() ); 
   measure_obj.mp = simulator.mp;
   measure_obj.init( simulator.kTlist );
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
