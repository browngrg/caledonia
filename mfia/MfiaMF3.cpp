// MfiaMF3.cpp -- WangLandau Transtion-Matrix sampling of totally mean-field Mfia model
//                Changed exchange to a mean-field value
// Sept7, 2014
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// History:
// 
// Mfia3.cpp -- WangLandau Transtion-Matrix sampling of Mfia model
//              Brought Mfia Hamiltonian back in
// Sept 5, 2014
//
// Heisenberg2.cpp -- Changed from WL update to TM update
// Aug 29, 2014
//
// Heisenberg1.cpp -- WangLandau sampling of Heisenberg model
// Aug 15, 2014
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
#include"MfiaMF_Hamiltonian.hpp"
#include"Random.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>


void DoHeisenberg(int& argc, char* argv[])
{

   //  Get basic MPI information
   int iproc = 0;
   int nproc = 1;
#  ifdef USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   using namespace std;
   ProgramOptions options("MfiaMF3","Transition-Matrix Wang-Landau Sampling for complete MF model.");
 
   // Construct a hamiltonian object 
   typedef MfiaMF_Hamiltonian HAMILTONIAN;
   HAMILTONIAN hamilton;
   hamilton.D     = 2;
   hamilton.L     = 64;
   hamilton.JSIGN = -1;
   hamilton.H     = 0;
   hamilton.A     = 0;
   options.add_option("dim",    "dimension of spin lattice",     'D', &(hamilton.D) );
   options.add_option("length", "extent of lattice in each dim", 'L', &(hamilton.L) );
   options.add_option("jsign",  "sign in front of exchange",     'J', &(hamilton.JSIGN) );
   options.add_option("H",      "magnitude of magnetic field",   'H', &(hamilton.H) );
   options.add_option("A",      "magnitude of mean field",       'A', &(hamilton.A) );

   // Construct the simulation object
   MC_WangLandau wanglandau;
   options.add_option( "Elo",     "lower side of energy window",    ' ', &(wanglandau.Elo));
   options.add_option( "Ehi",     "lower side of energy window",    ' ', &(wanglandau.Ehi));
   options.add_option( "numwin",  "number of windows",              ' ', &(wanglandau.NWindow));
   options.add_option( "numwalk", "number of walkers per window",   ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "nchunk",  "number of steps per test",       ' ', &(wanglandau.NChunk));
   options.add_option( "overlap", "fractional overlap of windows",  ' ', &(wanglandau.fwinover));
   options.add_option( "numbin",  "number of bins per window",      ' ', &(wanglandau.nbin));

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
   bool verbose = options.get_value<bool>("verbose") && (iproc==0);

   // Re-initialize objects
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Reinitialize begun" << endl;
   hamilton.init(verbose);
   const int NGrid = hamilton.V;
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " NGrid=" << NGrid << endl;

   // seed the random number geneator
   wanglandau.urng.seed( ParallelSeed(SeedFromClock()) );
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(wanglandau.urng,rng_output,std::cout);
   if( rng_failed )
   {
      cout << "Problem detected with random number generator" << endl;
      return;
   }
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Random number generator created" << endl;

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Create one walker" << endl;
   typedef WL_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   walkerpool[0].sigma.resize(NGrid);
   hamilton.initial_mixed(walkerpool[0].sigma);
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   wanglandau.verbose = options.get_value<bool>("verbose");
   wanglandau.init_pool(walkerpool);

   if(verbose) cout << __FILE__ << ":" << __LINE__ << " walker energy=" << walkerpool[0].now.E << endl;
   wanglandau.DoTMConverge(hamilton,walkerpool);

}



int main(int argc, char* argv[])
{
   if( true )
   {
      std::cout << "cmdline:";
      for(int i=0; i<argc; i++) std::cout << " " << argv[i];
      std::cout << std::endl;
   }

#  ifdef USE_MPI
   MPI_Init(&argc,&argv);
#  endif
   DoHeisenberg(argc,argv);
#  ifdef USE_MPI
   MPI_Finalize();
#  endif

}

