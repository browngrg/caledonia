
#include"ITTM.hpp"

#include<iostream>
#include<vector>


int main(int argc, char* argv[])
{
   // Read the Transition matrix, 
   std::vector<long long> C,Cob;
   int CWIDTH;
   read_C(std::cin,C,Cob,CWIDTH);
   // Calculate the DOS from the TM
   std::vector<double> SITTM;
   calc_SITTM(C,Cob,CWIDTH,SITTM);
   // Output the DOS
   float Ebin = 1;
   if(argc>1) sscanf(argv[1],"%f",&Ebin);
   float Emin = 0;
   if(argc>2) sscanf(argv[2],"%f",&Emin);
   std::cout << "# DOS from TM stored in write_C" << std::endl;
   std::cout << "# ITTM_DOS.exe [Ebin] [Emin] < CFILE.txt > DOSFILE.csv" << std::endl;
   std::cout << "# Column 1: E" << std::endl;
   std::cout << "# Column 2: DOS" << std::endl;
   std::cout << "# Column 3: ibin" << std::endl;
   for(int i=0; i<Cob.size(); i++)
      std::cout << Emin+i*Ebin << " " << SITTM[i] << " " << i << std::endl;
}
