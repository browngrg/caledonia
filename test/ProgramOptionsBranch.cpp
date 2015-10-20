
#include"ProgramOptions.hpp"

#include<cmath>
#include<string>
#include<iostream>
#include"mpi.h"


void do_ehmodel(ProgramOptions& options, int argc, char* argv[])
{
   int iproc_world;
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc_world);
   int nproc_world;
   MPI_Comm_size(MPI_COMM_WORLD,&nproc_world);
   double k = 3;
   options.add_option("ehk","ehmodel K",' ',&k);
   options.parse_command_line(argc,argv);
   std::cout << "ehmodel " << iproc_world << "/" << nproc_world << " argc=" << argc << std::endl;
}


void do_ising(ProgramOptions& options, int argc, char* argv[])
{
   int iproc_world;
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc_world);
   int nproc_world;
   MPI_Comm_size(MPI_COMM_WORLD,&nproc_world);
   int  ij = 30;
   options.add_option("ijk","ising model K",' ',&ij);
   options.parse_command_line(argc,argv);
   std::cout << "Ising " << iproc_world << "/" << nproc_world << " argc=" << argc << std::endl;
}

int main(int argc, char* argv[])
{

   MPI_Init(&argc,&argv);

   ProgramOptions options("Test_ProgramOptions2","Test with branching");

   char s1[20] = "ehmodel";
   char s2[20] = "wanglandau";

   options.add_option("model", "hamiltonian", ' ', s1);
   options.add_option("sim",   "simulation",  ' ', s2);

   options.parse_command_line2(argc,argv);

   if( strcmp(s1,"ehmodel")==0 )
      do_ehmodel(options,argc,argv);
   else if( strcmp(s1,"ising")==0 )
      do_ising(options,argc,argv);
   else
      std::cout << "Did not recognize model=\"" << s1 << "\"" << std::endl;

   MPI_Finalize();
}
