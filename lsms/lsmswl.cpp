// lsmswl.cpp -- WangLandau driver for LSMS DFT energy calculator
//
// October 13, 2015  Started importing LSMS_Class.hpp to this wrapper
//
// Greg Brown (browngrg@comcast.net,gbrown@fsu.edu)
//


#include"LSMS_Hamiltonian.hpp"
#include"ProgramOptions.hpp"

#include<iostream>
#include<string.h>

void DoLSMS(int& argc, char* argv[])
{
   // Get basic MPI information
   int iproc = 0;
   int nproc = 1;
#  ifdef USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#  endif

   ProgramOptions options("lsmswl","Wang-Landau Sampling");

   // Construct the hamiltonian object
   typedef LSMS_Hamiltonian HAMILTONIAN;
   HAMILTONIAN hamilton;
   hamilton.add_options(options);

   options.parse_command_line(argc,argv);

   hamilton.old_dosim();
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
   DoLSMS(argc,argv);
#  ifdef USE_MPI
   MPI_Finalize();
#  endif
}
