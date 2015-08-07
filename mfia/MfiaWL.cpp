// Mfia7.cpp -- WangLandau methods for Ising AFM with LR Term
//
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//
// History:
// 
// Mfia6.cpp -- Wang-Landau sampling of Ising models with LR Term
// June 2015
//
// Mfia5.cpp -- Fixed-magnetization sampling (never implemented)
// Used for convergence testing (mfiaconverge)
// May 2015 
//
// Mfia4.cpp -- WangLandau Transtion-Matrix sampling of Mfia model
// Using for more testing of Transition-Matrix methods
// October 29, 2014
//
// Mfia3.cpp -- WangLandau Transtion-Matrix sampling of Mfia model
// Brought Mfia Hamiltonian back in
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
//

#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"Mfia_Hamiltonian.hpp"
#include"Random.hpp"
#include"EMX_Measure.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>
#include<ctime>


void DoIt(int& argc, char* argv[])
{

   //  Get basic MPI information
   int iproc = 0;
   int nproc = 1;
#  ifdef USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   using namespace std;
   ProgramOptions options("Mfia6","Wang-Landau Sampling of Ising model with LR term.");
 
   // Construct a hamiltonian object 
   typedef Mfia_Hamiltonian HAMILTONIAN;
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

   // Construct the sampling object
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

   // Local options
   char lng_est_fn[512]; lng_est_fn[0]=0;
   options.add_option( "dos",   "dos file to use in sampling",      ' ', lng_est_fn);

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
   bool verbose = options.get_value<bool>("verbose") && (iproc==0);

   // Re-initialize objects
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Reinitialize begun" << endl;
   hamilton.init(verbose);
   const int NGrid = hamilton.V;
   std::cout << hamilton.header();
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " NGrid=" << NGrid << endl;

   // seed the random number geneator
   wanglandau.urng.seed( ParallelSeed(SeedFromClock()) );
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(wanglandau.urng,rng_output,std::cout);
   if( rng_failed )
   {
      cout << __FILE__ << ":" << __LINE__ << "Problem detected with random number generator (returns 1?)" << endl;
   }
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Random number generator created" << endl;

   // The measurement object
   EMX_Measure measure_obj;

   // Read in the density of state
   std::vector<double> energy,lng_est;
   if( lng_est_fn[0]!=0 )
   {
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
   if( energy.size()==0 )
   {
      const int npt = 100;
      energy.resize(npt);
      lng_est.resize(npt);
      double de = (wanglandau.Ehi-wanglandau.Elo)/static_cast<double>(npt-1);
      for(int i=0; i<npt; i++)
      {
         energy[i] = wanglandau.Elo + i*de;
         lng_est[i] = 0;
      }
   }
   if( wanglandau.Elo<std::numeric_limits<double>::max() )
   {
      int i=0; 
      while( i<energy.size() && energy[i]<wanglandau.Elo ) i++;
      if( i<energy.size() )
      { 
         std::cout << "Elo = " << wanglandau.Elo << std::endl;
         std::cout << "Triming energy from " << energy.front() << " with " << energy.size() << " elements" << std::endl;
         energy.erase(energy.begin(),energy.begin()+i);
         lng_est.erase(lng_est.begin(),lng_est.begin()+i);
         energy[0] = wanglandau.Elo;
         std::cout << "To energy from " << energy.front() << " with " << energy.size() << " elements" << std::endl;
      }
   }
   if( wanglandau.Ehi>-std::numeric_limits<double>::max() )
   {
      int i=0; 
      while( i<energy.size() && energy[i]<wanglandau.Ehi ) i++;
      if( i<energy.size() )
      { 
         std::cout << "Ehi = " << wanglandau.Ehi << std::endl;
         std::cout << "Triming energy ending at " << energy.back() << " with " << energy.size() << " elements" << std::endl;
         energy.erase(energy.begin()+i,energy.end());
         lng_est.erase(lng_est.begin()+i,lng_est.end());
         energy[i-1] = wanglandau.Ehi;
         std::cout << "To energy ending at " << energy.back() << " with " << energy.size() << " elements" << std::endl;
      }
   }

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
   wanglandau.partition_windows(energy,lng_est);
   wanglandau.init_pool(walkerpool);   
   for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
      walkerpool[iwalk].set_fixed(energy,lng_est);
   //measure_obj.init(walkerpool[0].window.Elo,walkerpool[0].window.Ehi,walkerpool[0].window.Ebin);
   measure_obj.init(wanglandau.Elo,wanglandau.Ehi,wanglandau.Ebin);
   measure_obj.mp_window = wanglandau.mp_window; 
   wanglandau.DoConverge(hamilton,walkerpool,measure_obj);
   //wanglandau.DoTMConverge(hamilton,walkerpool);
}



int main(int argc, char* argv[])
{
   const bool very_verbose = false;
   time_t start_time_main;
   if( very_verbose )
   {
      std::cout << "cmdline:";
      for(int i=0; i<argc; i++) std::cout << " " << argv[i];
      std::cout << std::endl;
      time(&start_time_main);
      struct tm * timeinfo = localtime(&start_time_main);
      std::cout << "start time = " << asctime(timeinfo); // << std::endl;
      int pid = 0;
      if(true) pid = getpid();
      std::cout << "process id = " << pid << std::endl;
   }



#  ifdef USE_MPI
   MPI_Init(&argc,&argv);
#  endif
   DoIt(argc,argv);
#  ifdef USE_MPI
   MPI_Finalize();
#  endif


   time_t stop_time_main;
   if( very_verbose )
   {
      time(&stop_time_main);
      double secs = difftime(stop_time_main,start_time_main);
      clock_t tics = clock();
      double psecs = static_cast<double>(tics)/static_cast<double>(CLOCKS_PER_SEC);
      struct tm* timeinfo = localtime(&stop_time_main);
      std::cout << "end time = " << asctime(timeinfo); // << std::endl;
      std::cout << "wall time = " << secs << std::endl;
      std::cout << "process time = " << psecs << std::endl; 
#     if 1
      int tSize = 0, resident = 0, share = 0;   // In kilobytes
      if(true) 
      {
         std::ifstream buffer("/proc/self/statm");
         buffer >> tSize >> resident >> share;
         buffer.close();
         long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
         resident = resident * page_size_kb;
         share = share * page_size_kb;
         std::cout << "resident memory = " << resident << " KBytes" << std::endl;
      }
#     endif
   }


}
