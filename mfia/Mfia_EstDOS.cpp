
#include"Mfia_Hamiltonian.hpp"

#include<iostream>


int main(int argc, char* argv[])
{
   Mfia_Hamiltonian model;

   model.D = 2;
   model.L = 128;
   model.H = 0;
   model.A = 7;
   model.init();

   int nbin = 500;
   std::vector<double> energy,lng;
   model.estimate_lng(energy,lng,nbin);

   for(int i=0; i<nbin; i++)
      std::cout << " " << energy[i] << " " << lng[i] << " " << i << std::endl;
}
