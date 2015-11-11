
#include"ProgramOptions.hpp"

#include<cmath>
#include<string>
#include<iostream>
#include"mpi.h"


int main(int argc, char* argv[])
{

   MPI_Init(&argc,&argv);

   ProgramOptions options("Test_ProgramOptions","Test MPI Bcast of options");

   bool b;
   int i;
   long l;
   float f;
   float d;

   char s1[10] = "smurf",s2[20],s3[10];

   options.add_option("bool",   "boolean variable",      ' ', &b);
   options.add_option("int",    "integer variable",      ' ', &i);
   options.add_option("long",   "long int variable",     ' ', &l);
   options.add_option("float",  "float variable",        ' ', &f);
   options.add_option("double", "double prec variable",  ' ', &d);
   options.add_option("str1",   "first string variable" ,' ', s1);
   options.add_option("str2",   "second string variable",' ', s2);
   options.add_option("str3",   "third string variable" ,' ', s3);

   int iproc_world;
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc_world);
   int nproc_world;
   MPI_Comm_size(MPI_COMM_WORLD,&nproc_world);

   b = ((iproc_world%2)==0);
   i = 10+iproc_world;
   l = 100*(10+iproc_world);
   f = 3.141529*(1+iproc_world);
   d = std::exp(f);

   std::strcpy(s1,"dummy");
   std::strcpy(s2,"dummy");
   std::strcpy(s3,"dummy");
  
   if( iproc_world==0 )
   {
      std::strcpy(s1,"first");
      std::strcpy(s2,"second");
      std::strcpy(s3,"thrid"); 
   }

   options.parse_command_line(argc,argv);

   for(int ip=0; ip<nproc_world; ip++)
   {
      if( ip==iproc_world )
         std::cout << iproc_world << ": " << b << " " << i << " " << l << " " << f << " " << d << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
   }
   for(int ip=0; ip<nproc_world; ip++)
   {
      if( ip==iproc_world )
         std::cout << iproc_world << ": \"" << s1 << "\" \"" << s2 << "\" \"" << s3 << "\"" << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
   }

   MPI_Finalize();
}
