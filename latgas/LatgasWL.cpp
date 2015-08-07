// Latgas3.cpp -- Wang-Landau convergence of Lattice Gas model
//
// invoke compiler with -DUSE_MPI 
// Look at ProgamOptions options object to see command-line options for setting windows, etc
//
// History:
//
// Aug 6, 2014 Latgas3.cpp -- Add ParallelSeed for random numbers
// Greg Brown (gbrown@fsu.edu,browngrg@comcst.net)
//
// June 17, 2014 Latgas2.cpp -- Add Wang-Landau functionality
// Greg Brown (gbrown@fsu.edu,browngrg@comcst.net)
//
// June 16, 201 4Latgas1.cpp -- Make Latgas_Hamiltonian a model of the Hamiltonian concept
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// Wham3.cpp -- Generate the trajectory of one Metropolis Monte Carlo walker
//
// June 4, 2014    Refactored for larger MC scheme
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//


#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"Latgas_Hamiltonian.hpp"
#include"Latgas_ReadConfig.hpp"
#include"Random.hpp"
#include"ProgramOptions.hpp"


#include<stdio.h>


void WL_Latgas(int argc, char* argv[])
{
   // Need a local scope so MPI cleanup in destructors get called before MPI Finalize

   using namespace std;
   ProgramOptions options("Latgas","Wang-Landau Monte Carlo of lattice-gas models");

   int iproc = 0;
   int nproc = 1;
#  ifdef USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   // Construct the simulation object
   MC_WangLandau wanglandau;
   options.add_option( "Elo",     "lower side of global energy window",    ' ', &(wanglandau.Elo));
   options.add_option( "Ehi",     "lower side of global energy window",    ' ', &(wanglandau.Ehi));
   options.add_option( "Ebin",    "Width of energy bins",                  ' ', &(wanglandau.Ebin));
   options.add_option( "numwin",  "number of windows",                     ' ', &(wanglandau.NWindow));
   options.add_option( "numwalk", "number of walkers per window",          ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "overlap", "fractional overlap of windows",         ' ', &(wanglandau.fwinover));
   options.add_option( "nchunk",  "number of steps per test",              ' ', &(wanglandau.NChunk));
   options.add_option( "Q",        "target convergence factor",            ' ', &(wanglandau.Qquit));
   options.add_option( "maxupdate","maximum number of iterations",         ' ', &(wanglandau.MaxUpdate));
   options.add_option( "dosinterp","linear interpolation of dos",          ' ', &(wanglandau.LinearInterp));
   options.add_option( "wleta",    "weighting between WL and ITTM",        ' ', &(wanglandau.wleta));
   options.add_option( "wlgamma",  "starting value of WL parameter",       ' ', &(wanglandau.wlgamma_start));
   options.add_option( "configout","output configurations to disk",        ' ', &(wanglandau.output_configs));

   // Get user input
   options.parse_command_line(argc,argv);

   // seed the random number geneator
   wanglandau.urng.seed( ParallelSeed(SeedFromClock()) );
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(wanglandau.urng,rng_output,std::cout);
   if( rng_failed )
   {
      cout << "Problem detected with random number generator" << endl;
      return;
   }

   // Construct a hamiltonian object 
   typedef Latgas_Hamiltonian HAMILTONIAN;
   Latgas_Hamiltonian hamilton("Latgas_ModelIn.txt");

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   typedef WL_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   ReadStrOut("str.out",hamilton,walkerpool[0].sigma);
   int NGrid = hamilton.grid.size();
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   wanglandau.verbose = options.get_value<bool>("verbose");
   wanglandau.init_pool(walkerpool);

   wanglandau.DoWLConverge(hamilton,walkerpool);
}


int main(int argc, char* argv[])
{
   // Main only takes care of MPI and gaurantees that destructors (which make clean
   // up MPI by freeing groups and communicators) get called before MPI_Finalize
   MPI_Init(&argc,&argv);
   WL_Latgas(argc,argv);         
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
}
