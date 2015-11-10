// HeisenbergDiffuse.cpp -- Measure diffusion properties of WangLandau walkers
//
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)

#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"MPI_Struct.hpp"
#include"Heisenberg_Hamiltonian.hpp"
#include"Meas_Diffusion.hpp"
#include"Null_Measure.hpp"
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
   hamilton.L     = 16;
   hamilton.H     = 0;
   options.add_option("length", "extent of lattice in each dim", 'L', &(hamilton.L) );
   options.add_option("H",      "magnitude of magnetic field",   'H', &(hamilton.H) );

   // Construct the simulation object
   MC_WangLandau wanglandau;
   wanglandau.Elo   = -11000;
   wanglandau.Ehi   = +11000;
   wanglandau.Ebin  = 1;
   wanglandau.wall_limit = 5*24*60*60;
   options.add_option( "Elo",     "lower side of energy window",    ' ', &(wanglandau.Elo));
   options.add_option( "Ehi",     "lower side of energy window",    ' ', &(wanglandau.Ehi));
   options.add_option( "Ebin",    "width of energy bins",           ' ', &(wanglandau.Ebin));
   options.add_option( "numwin",  "number of windows",              ' ', &(wanglandau.NWindow));
   options.add_option( "numwalk", "number of walkers per window",   ' ', &(wanglandau.NWalkPerProcess));
   options.add_option( "overlap", "fractional overlap of windows",  ' ', &(wanglandau.fwinover));
   options.add_option( "nstep",    "number of steps per iteration", ' ', &(wanglandau.NStep));
   options.add_option( "maxupdate","maximum number of iterations",  ' ', &(wanglandau.MaxUpdate));
   options.add_option( "dosinterp","linear interpolation of dos",   ' ', &(wanglandau.LinearInterp));
   options.add_option( "wleta",    "weighting between WL and ITTM", ' ', &(wanglandau.wleta));
   options.add_option( "convh",    "converge using visit histogrma",' ', &(wanglandau.convh));
   options.add_option( "wlgamma",  "starting value of WL parameter",' ', &(wanglandau.wlgamma_start));
   options.add_option( "Q",        "target convergence factor",     ' ', &(wanglandau.Qquit));
   options.add_option( "configout","output configurations to disk", ' ', &(wanglandau.output_configs));
   options.add_option( "walllimit","maximum run time",              ' ', &(wanglandau.wall_limit));

#  undef MEASURE_DIFFUSION
#  ifdef MEASURE_DIFFUSION
   Meas_Diffusion measure;
   options.add_option( "diff_delay","delay in measuring diffusion", ' ', &(measure.delay));
#  else
   Null_Measure measure;
#  endif

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
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Random number generator created" << endl;

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Create one walker" << endl;
   typedef WL_Walker<HAMILTONIAN::Config,HAMILTONIAN::Observables> Walker;
   std::vector<Walker> walkerpool(1);
   walkerpool[0].sigma.resize(NGrid);
   hamilton.initial_random(walkerpool[0].sigma,wanglandau.urng);
   hamilton.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
   wanglandau.verbose = verbose;
   wanglandau.mp_window.pool = MPI_Struct::world();
   wanglandau.partition_windows();
   wanglandau.init_pool(walkerpool);
   std::vector<double> energy;
   std::vector<double> est_lng;
   ReadLngEst(std::string(lng_est_fn),energy,est_lng);
   if( false && energy.size()<3 )
   {
      energy.resize(1000);
      est_lng.resize(1000);
      // Set the fixed DOS to spin-wave part (N-1)*( ln(1+x)+ln(1-x) )
      double N = static_cast<double>( walkerpool[0].now.V );    // Volume of the system
      double Eb = -3*N;
      double Et = +3*N;
      double de = (Et-Eb)/1000; 
      for(int i=0; i<1000; i++)
      {
          double x = -1 + 2*(static_cast<double>(i+1)/1001.);
          if( x>-1 && x<1 )
          {
             energy[i] = Eb + (x+1.)*(Et-Eb)/2.;
             est_lng[i] = (N-1.)*(std::log(1+x)+std::log(1-x));
          }
      }
   }
   if( energy.size()>0 )
   {
      TrimLng(wanglandau.Elo,wanglandau.Ehi,energy,est_lng); 
      wanglandau.partition_windows(energy,est_lng);
      for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
         walkerpool[iwalk].set_fixed(energy,est_lng);
   }
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Done create one walker" << endl;

#  ifdef MEASURE_DIFFUSION
   measure.init(wanglandau.Elo,wanglandau.Ehi,wanglandau.Ebin);
   if(verbose) cout << __FILE__ << ":" << __LINE__ << " Measurement object initialized" << endl;
#  endif

   wanglandau.DoConverge(hamilton,walkerpool,measure);

   options.write();
}



int main(int argc, char* argv[])
{

   int iproc = 0;
#  ifdef USE_MPI
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
#  endif
   if( iproc==0 )
   {
      std::cout << "cmdline:";
      for(int i=0; i<argc; i++) std::cout << " " << argv[i];
      std::cout << std::endl;
   }
   DoHeisenberg(argc,argv);
#  ifdef USE_MPI
   MPI_Finalize();
#  endif

}

