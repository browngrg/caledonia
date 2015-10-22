// IsingWL.cpp -- WangLandau convergence of the Ising model
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)

#include"MC_WangLandau.hpp"
#include"WL_Walker.hpp"
#include"../mfia/Mfia_Hamiltonian.hpp"
#include"../ehmodel/EHModel_Hamiltonian.hpp"
#include"MPI_Struct.hpp"
#include"Random.hpp"
#include"EMX_Measure.hpp"
#include"Null_Measure.hpp"
#include"ProgramOptions.hpp"

#include<stdio.h>
#include<ctime>
#include<algorithm>


template<typename HAMILTON, typename MEASURE>
void WangLandauDriver(ProgramOptions& options, int& argc, char* argv[])
{

   //  Get basic MPI information
   MPI_Struct world;
   world = MPI_Struct::world();
#  ifdef USE_MPI
   MPI_Comm_rank(world.comm,&world.iproc);
   MPI_Comm_size(world.comm,&world.nproc);
#  endif
 
   bool verbose = options.get_value<bool>("verbose") && (world.iproc==0);
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Construction objects" << std::endl;

   // Construct a hamiltonian object 
   HAMILTON hamilton;
   hamilton.add_options(options);

   // Construct the sampling object
   MC_WangLandau simulation;
   simulation.add_options(options);

   typedef WL_Walker<typename HAMILTON::Config,typename HAMILTON::Observables> Walker;
   std::vector<Walker> walkerpool(1);

   // The measurement object
   // EMX_Measure measure_obj;
   MEASURE measure_obj;

   // Get parameters
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Parsing parameters" << std::endl;
   options.parse_command_line(argc,argv);
   if( options.get_value<bool>("help") ) return;

   // seed the random number geneator
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Creating random number generator" << std::endl;
   simulation.urng.seed( ParallelSeed(SeedFromClock()) );
   bool rng_output = false;
   bool rng_failed = RNGTestMoments(simulation.urng,rng_output,std::cout);
   if( false && rng_failed )
      std::cout << __FILE__ << ":" << __LINE__ << "Problem detected with random number generator (returns 1?)" << std::endl;

   // Initialize objects based on ProgramOptions
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Call Init routines" << std::endl;
   hamilton.init(verbose);
   hamilton.initial(walkerpool[0].sigma);
   simulation.init(hamilton,walkerpool,verbose);
   measure_obj.init(simulation);

   // Construct a representation of the model
   // This includes the "microscopic" configuration sigma_i
   // and the macroscopic quantities like magnetization and energy
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Beginning simulation" << std::endl;
   simulation.DoConverge(hamilton,walkerpool,measure_obj);

   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Saving options" << std::endl;
   options.write();
}



int main(int argc, char* argv[])
{
   const bool very_verbose = true;
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

   //  Get basic MPI information
   MPI_Struct world;
#  ifdef USE_MPI
   MPI_Init(&argc,&argv);
   world = MPI_Struct::world();
#  endif

   ProgramOptions options("caledonia","Monte Carlo simulations");
   char model_name[128] = "ising";
   char sim_name[128] = "wanglandau";
   options.add_option("model","Hamiltonian type", ' ', model_name);
   options.add_option("sim",  "Monte Carlo simulation type", ' ', sim_name);
   options.parse_command_line2(argc,argv);
   if( strcmp(sim_name,"wanglandau")==0 )
   {
      if( strcmp(model_name,"ising")==0 )   {  WangLandauDriver<Mfia_Hamiltonian,EMX_Measure>(options,argc,argv); }
      else if( strcmp(model_name,"ehmodel")==0 ) {  WangLandauDriver<EHModel_Hamiltonian,EMX_Measure>(options,argc,argv); }
      else { if(world.iproc==0) std::cout << "model \"" << model_name << "\" not recognized for \"" << sim_name << "\"" << std::endl; }
   }


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
