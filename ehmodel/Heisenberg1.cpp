// Heisenberg1.cpp -- WangLandau sampling of Heisenberg model
// Aug 15, 2014
//
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// History:
//
// Mfia2.cpp -- WangLandau sampling of Mfia model (Borrows heavily from Latgas3.cpp)
// August 7, 2014
//
// Mfia1.cpp -- Generate the Metropolis trajectory of one Mfia Monte Carlo walker
// June 28, 2014  Copied from Wham3.cpp
//
// Wham3.cpp -- Generate the trajectory of one Metropolis Monte Carlo walker
// June 4, 2014    Refactored for larger MC scheme
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//

#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"Heisenberg_Hamiltonian.hpp"
#include"Random.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>




int main(int argc, char* argv[])
{
   if( true )
   {
      std::cout << "cmdline:";
      for(int i=0; i<argc; i++) std::cout << " " << argv[i];
      std::cout << std::endl;
   }

   int iproc = 0;
   int nproc = 1;
#  ifdef USE_MPI
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   using namespace std;
   ProgramOptions options("Heisenberg","Wang-Landau Sampling.");
 
   // Construct a hamiltonian object 
   typedef Heisenberg_Hamiltonian HAMILTONIAN;
   HAMILTONIAN hamilton;
   hamilton.L     = 64;
   hamilton.H     = 0;
   options.add_option("length", "extent of lattice in each dim", 'L', &(hamilton.L) );
   options.add_option("H",      "magnitude of magnetic field",   'H', &(hamilton.H) );

   // Construct the simulation object
   MC_WangLandau wanglandau;
   options.add_option( "Elo",     "lower side of energy window",    ' ', &(wanglandau.Elo));
   options.add_option( "Ehi",     "lower side of energy window",    ' ', &(wanglandau.Ehi));
   options.add_option( "numwin",  "number of windows",              ' ', &(wanglandau.NWindow));
   options.add_option( "numwalk", "number of walkers per window",   ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "overlap", "fractional overlap of windows",  ' ', &(wanglandau.fwinover));
   options.add_option( "numbin",  "number of bins per window",      ' ', &(wanglandau.nbin));
   options.add_option( "configout","write some configs to disk",    ' ', &(wanglandau.output_configs));


   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return 0;
   bool verbose = options.get_value<bool>("verbose") && (iproc==0);

   // Re-initialize objects
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Reinitialize begun" << endl;
   hamilton.init(verbose);
   const int NGrid = hamilton.nspin;
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " NGrid=" << NGrid << endl;
   //simulator.NStep = static_cast<int>( static_cast<float>(NGrid)*NStep_MCSS );

   // seed the random number geneator
   wanglandau.urng.seed( ParallelSeed(SeedFromClock()) );
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(wanglandau.urng,rng_output,std::cout);
   if( rng_failed )
   {
      cout << "Problem detected with random number generator" << endl;
      return 1;
   }
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Random number generator created" << endl;

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Create one walker" << endl;
   typedef WL_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   walkerpool[0].sigma.resize(NGrid);
   hamilton.initial_ferro(walkerpool[0].sigma);
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   wanglandau.verbose = options.get_value<bool>("verbose");
   wanglandau.init_pool(walkerpool);

   if(verbose) cout << __FILE__ << ":" << __LINE__ << " walker energy=" << walkerpool[0].now.E << endl;

   wanglandau.DoWLConverge(hamilton,walkerpool);

#  ifdef USE_MPI
   MPI_Finalize();
#  endif

}

