// MfiaWLfixed.cpp -- WangLandau methods for Ising AFM with LR Term
//
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
//

#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"Mfia_Hamiltonian.hpp"
#include"Random.hpp"
#include"EMX_Measure.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>
#include<iomanip>
#include<fstream>
#include<ctime>


void DoIt(int& argc, char* argv[])
{

   //  Get basic MPI information
   MPI_Struct world = MPI_Struct::world();
#  ifdef USE_MPI
   MPI_Comm_rank(world.comm,&world.iproc);
   MPI_Comm_size(world.comm,&world.nproc);
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

   // Construct the sampling object
   MC_WangLandau wanglandau;
   wanglandau.wleta = 0;
   wanglandau.wlgamma_start = 1;
   wanglandau.Qquit = 0.05;
   wanglandau.MaxUpdate = 15;
   wanglandau.NWindow = 1;
   wanglandau.wall_limit = 2*24*60*60;
   options.add_option( "numwalk",   "number of walkers per window",   ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "nstep",     "number of steps per iteration",  ' ', &(wanglandau.NStep));
   options.add_option( "maxupdate", "maximum number of iterations",   ' ', &(wanglandau.MaxUpdate));
   options.add_option( "wleta",     "weighting between WL and ITTM",  ' ', &(wanglandau.wleta));
   options.add_option( "wlgamma",   "starting value of WL parameter", ' ', &(wanglandau.wlgamma_start));
   options.add_option( "Q",         "target convergence factor",      ' ', &(wanglandau.Qquit));
   options.add_option( "wall_limit","wall_limit limit for run",       ' ', &(wanglandau.wall_limit));
   wanglandau.LinearInterp = 0;

   // Local options
   int NRes = 32;
   options.add_option("res", "Number sampling pts in NA,NB", ' ', &(NRes));

   // Get parameters
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;
   bool verbose = options.get_value<bool>("verbose") && (world.iproc==0);

   // Re-initialize objects
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Reinitialize begun" << endl;
   hamilton.init(verbose);
   const int NGrid = hamilton.V;
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " NGrid=" << NGrid << endl;
   
   // Set the energy ranges
   hamilton.op_constrained = true;
   wanglandau.Elo = -2*NGrid - 1;  // Center bin on possible energy
   wanglandau.Ehi = +2*NGrid + 1;
   wanglandau.Ebin = 2;

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
   Null_Measure measure_obj;

   // Construct a zero density of states
   std::vector<double> energy,lng_est;
   const int npt = 100;
   energy.resize(npt);
   lng_est.resize(npt);
   double de = (wanglandau.Ehi-wanglandau.Elo)/static_cast<double>(npt-1);
   for(int i=0; i<npt; i++)
   {
      energy[i] = wanglandau.Elo + i*de;
      lng_est[i] = 0;
   }

   // Potential threading
   int ithread = 0;
   int nthread = 1;
#  ifdef USE_MPI
   {
      MPI_Struct world = MPI_Struct::world();
      ithread = world.iproc;
      nthread = world.nproc;
   }
#  endif
#  ifdef USE_CONDOR
   sscanf(argv[1],"%d",&ithread);
   sscanf(argv[2],"%d",&nthread);
#  endif

   std::vector< std::pair<int,int> > work;
   char work_fn[100] = "MfiaWLfixed-todo.csv";
   if( access(work_fn,F_OK)!=-1 )
   {
      std::ifstream fin(work_fn);
      while( !fin.eof() )
      {
         std::string buffer;
         int ThetaA,ThetaB;
         fin >> ThetaA >> ThetaB;
         if( !fin.eof() )
         {
            work.push_back( std::pair<int,int>(ThetaB,ThetaA ) );
            std::getline(fin,buffer);
         }
      }
   }
   else 
   {
      int stride = NGrid/4/NRes;
      for(int ThetaB=stride; ThetaB<=(NGrid/4); ThetaB+=stride)
      {
         for(int ThetaA=ThetaB; ThetaA<=(NGrid/4); ThetaA+=stride)
         {
            work.push_back( std::pair<int,int>(ThetaB,ThetaA ) );
         }
      }
   }

   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Create one walker" << endl;
   typedef WL_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;

   wanglandau.mp_window.set_pool(MPI_Struct::local());
   std::vector<Walker> walkerpool(1);
   walkerpool[0].sigma.resize(NGrid);
   for(int iwork=ithread; iwork<work.size(); iwork+=nthread)
   {
      int ThetaB = work[iwork].first;
      int ThetaA = work[iwork].second;
      char fname[100];
      sprintf(fname,"MfiaWLfixedHist-%03d-%03d.csv",ThetaB,ThetaA);
      std::ofstream fout(fname);
      fout << "# WL Histograms in the M-Phi plane" << std::endl;
      fout << "# L=" << hamilton.L << " N=" << NGrid << " H=" << hamilton.H << " A=" << hamilton.A << std::endl;
      fout << "# Column 1: N1A=ThetaA, upspins on A Sublattice" << std::endl;
      fout << "# Column 2: N1B=ThetaB, upspins on B Sublattice" << std::endl;
      fout << "# Column 3: M = -N + 2(ThetaA+ThetaB)" << std::endl;
      fout << "# Column 4: Phi = 2(ThetaA-ThetaB)" << std::endl;
      fout << "# Column 5: Neighbor Energy 2N-2U (U unlike bonds)" << std::endl;
      fout << "# Column 6: Inferred Ungerade pairs (U=N-E/2)" << std::endl;
      fout << "# Column 7: lnh, normalized to max=1" << std::endl;
      fout << "# Column 8: lnh, normalized to sum=1" << std::endl;
      Walker& walker(walkerpool[0]);
      std::fill(walker.sigma.begin(),walker.sigma.end(),-1);
      int nflip=0;
      for(int isite=0; isite<NGrid && nflip<ThetaA; isite++)
      {
         if( hamilton.stagmask[isite]>0 )
         {
            walker.sigma[isite]=+1;
            nflip++;
         }
      }
      nflip=0;
      for(int isite=0; isite<NGrid && nflip<ThetaB; isite++)
      {
         if( hamilton.stagmask[isite]<0 )
         {
            walker.sigma[isite]=+1;
            nflip++;
         }
      }
      hamilton.calc_observable(walker.sigma,walker.now);
      std::cout << "M=" << walker.now.M << "=" << -NGrid+2*(ThetaA+ThetaB) << " P="  << walker.now.P << "=" <<  2*(ThetaA-ThetaB) << std::endl;
      wanglandau.partition_windows(energy,lng_est);
      wanglandau.init_pool(walkerpool);   
      for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
         walkerpool[iwalk].set_fixed(energy,lng_est);
      wanglandau.DoConverge(hamilton,walkerpool,measure_obj);
      // Save results to a master file
      std::vector<double> Sadd( walker.S.size() );
      for(int i=0; i<walker.window.NBin; i++) Sadd[i] = walker.S[i] + walker.Sfixed[i];
      double Smax = -std::numeric_limits<double>::max();
      for(int i=0; i<walker.window.NBin; i++) 
         if( walker.h[i]>0 ) Smax = std::max(Smax,Sadd[i]);
      for(int i=0; i<walker.window.NBin; i++) Sadd[i] -= Smax;
      double hnorm = 0;
      for(int i=0; i<walker.window.NBin; i++) 
         if( walker.h[i]>0 ) hnorm += std::exp(Sadd[i]);
      hnorm = std::log(hnorm);
      for(int i=0; i<walker.window.NBin; i++)
      {
         if( walker.h[i]>0 )
         {
            int E = static_cast<int>(walker.window.unbin(i));
            int U = NGrid - E/2;
            fout << ThetaA << " " << ThetaB << " " << walker.now.M << " " << walker.now.P << " " << E << " " << U << " " 
                 << std::setprecision(12) <<  Sadd[i] << " " << std::setprecision(12) << Sadd[i]-hnorm << std::endl;
         }
      }
   }
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
