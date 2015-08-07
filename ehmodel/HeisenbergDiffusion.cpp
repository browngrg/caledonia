// HeisenbergDiffuse.cpp -- WangLandau sampling of the Heisenberg model
//
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// History
//
// HiesenbergDiffuse.cpp -- Added Meas_Diffusion object to measure diffusion of walkers.
// Aug 5, 2015
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
#include"Heisenberg_Hamiltonian.hpp"
#include"Meas_Diffusion.hpp"
#include"Random.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>

void ReadLngEst(std::string lng_est_fn, std::vector<double>& energy, std::vector<double>& lng_est)
{
      if( lng_est_fn[0]==0 ) return;
      if( true )
      {
         energy.resize(0);
         lng_est.resize(0);
         double E,lng;
         std::string buff;
         std::ifstream fin(lng_est_fn);
         while( fin && fin.is_open() && !fin.eof() )
         {
            std::getline(fin,buff);
            if( buff.size()>0 && buff[0]!='#' )
            {
               sscanf( buff.c_str(),"%lf %lf",&E,&lng);
               energy.push_back(E);
               lng_est.push_back(lng);
            }
         }
      }
}

void TrimLng(double Elo, double Ehi, std::vector<double>& energy, std::vector<double>& lng_est)
{
      if( Elo<std::numeric_limits<double>::max() )
      {
         int i=0; 
         while( i<energy.size() && energy[i]<Elo ) i++;
         if( i<energy.size() )
         { 
            std::cout << "Elo = " << Elo << std::endl;
            std::cout << "Triming energy from " << energy.front() << " with " << energy.size() << " elements" << std::endl;
            energy.erase(energy.begin(),energy.begin()+i);
            lng_est.erase(lng_est.begin(),lng_est.begin()+i);
            energy[0] = Elo;
            std::cout << "To energy from " << energy.front() << " with " << energy.size() << " elements" << std::endl;
         }
      }
      if( Ehi>-std::numeric_limits<double>::max() )
      {
         int i=0; 
         while( i<energy.size() && energy[i]<Ehi ) i++;
         if( i<energy.size() )
         { 
            std::cout << "Ehi = " << Ehi << std::endl;
            std::cout << "Triming energy ending at " << energy.back() << " with " << energy.size() << " elements" << std::endl;
            energy.erase(energy.begin()+i,energy.end());
            lng_est.erase(lng_est.begin()+i,lng_est.end());
            energy[i-1] = Ehi;
            std::cout << "To energy ending at " << energy.back() << " with " << energy.size() << " elements" << std::endl;
         }
      }
}


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
   options.add_option( "Ebin",    "width of energy bins",           ' ', &(wanglandau.Ebin));
   options.add_option( "numwin",  "number of windows",              ' ', &(wanglandau.NWindow));
   options.add_option( "numwalk", "number of walkers per window",   ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "overlap", "fractional overlap of windows",  ' ', &(wanglandau.fwinover));
   options.add_option( "nchunk",  "number of steps per iteration",  ' ', &(wanglandau.NChunk));
   options.add_option( "maxupdate","maximum number of iterations",  ' ', &(wanglandau.MaxUpdate));
   options.add_option( "dosinterp","linear interpolation of dos",   ' ', &(wanglandau.LinearInterp));
   options.add_option( "wleta",    "weighting between WL and ITTM", ' ', &(wanglandau.wleta));
   options.add_option( "wlgamma",  "starting value of WL parameter",' ', &(wanglandau.wlgamma_start));
   options.add_option( "Q",        "target convergence factor",     ' ', &(wanglandau.Qquit));
   options.add_option( "configout","output configurations to disk", ' ', &(wanglandau.output_configs));

   Meas_Diffusion measure;
   options.add_option( "delay",    "target convergence factor",     ' ', &(measure.delay));

   // Local options
   char lng_est_fn[512]; lng_est_fn[0]=0;
   options.add_option( "estdos",   "estimated dos file",            ' ', lng_est_fn);

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
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
   hamilton.initial_ferro(walkerpool[0].sigma);
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   wanglandau.verbose = options.get_value<bool>("verbose");
   wanglandau.partition_windows();
   wanglandau.init_pool(walkerpool);
   std::vector<double> energy;
   std::vector<double> est_lng;
   ReadLngEst(std::string(lng_est_fn),energy,est_lng);
   if( energy.size()>0 )
   {
      TrimLng(wanglandau.Elo,wanglandau.Ehi,energy,est_lng); 
      wanglandau.partition_windows(energy,est_lng);
      for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
         walkerpool[iwalk].set_fixed(energy,est_lng);
   }

   measure.init(wanglandau.Elo,wanglandau.Ehi,wanglandau.Ebin);

   wanglandau.DoConverge(hamilton,walkerpool,measure);
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

