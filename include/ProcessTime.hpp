#ifndef PROCESS_TIME_HPP_
#define PROCESS_TIME_HPP_

#define PROCESSTIME_USETIME_
#ifdef PROCESSTIME_USETIME_
#include<ctime>
#endif

#ifdef USE_MPI
#include"mpi.h"
#endif

#include "MPI_Struct.hpp"


class ProcessTime
{
public:

   ProcessTime() { mp=MPI_Struct::world(); start(); }

   void start();

   double elapsed();

   MPI_Struct mp;

private:

#  ifdef PROCESSTIME_USETIME_
   time_t start_time;
   time_t now_time;
#  endif
   double start_time_mpi;
   double elapsed_sec;

};


void ProcessTime::start()
{
#  ifdef USE_MPI
      start_time_mpi = MPI_Wtime();
#  else
#     ifdef PROCESSTIME_USETIME_
         time(&start_time);
#     endif
#  endif
}


double ProcessTime::elapsed()
{
   elapsed_sec = 0;
   // Get time
#  ifdef USE_MPI
      double now_time_mpi = MPI_Wtime();
      elapsed_sec = now_time_mpi - start_time_mpi;
      MPI_Bcast(&(elapsed_sec),1,MPI_DOUBLE,0,mp.comm);
#  else
#     ifdef PROCESSTIME_USETIME_
#     ifdef 1
         // Wall clock time
         time(&now_time);
         elapsed_sec = difftime(now_time,start_time);
#     else
         // Process time
         clock_t tics = clock();
         elapsed_sec = static_cast<double>(tics)/static_cast<double>(CLOCKS_PER_SEC);
#     endif
#     else
         std::cout << __FILE__ << ":" << __LINE__ << " No method for accessing time provided"
#     endif
#  endif
   return elapsed_sec;
}



#endif
