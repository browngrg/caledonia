#ifndef LATGAS_READCONFIG_HPP_
#define LATGAS_READCONFIG_HPP_


#include "Latgas_Hamiltonian.hpp"


void ReadStrOut(std::string fname, Latgas_Hamiltonian& hamilton, Latgas_Hamiltonian::Config& sigma)
{
   const bool verbose = false;
   typedef Vec<float,3> Vec3;
   Vec3 a[3];   // unit cell vector
   Vec3 rL[3];  // extent
   std::ifstream fin;
   try
   {
      fin.open(fname.c_str());
   }
   catch(...)
   {
      std::cout << "ReadStrOut could not open \"" << fname << "\"" << std::endl;
      return;
   }
   if( !fin || fin.eof() )
   {
      std::cout << "ReadStrOut could not open \"" << fname << "\"" << std::endl;
      return;
   }
   // read system box information
   fin >>  a[0][0] >>  a[0][1] >>  a[0][2];
   fin >>  a[1][0] >>  a[1][1] >>  a[1][2];
   fin >>  a[2][0] >>  a[2][1] >>  a[2][2];
   fin >> rL[0][0] >> rL[0][1] >> rL[0][2];
   fin >> rL[1][0] >> rL[1][1] >> rL[1][2];
   fin >> rL[2][0] >> rL[2][1] >> rL[2][2];
   // Test matches grid in all but extent?
   // assume cubic
   hamilton.grid.L[0] = static_cast<int>(rL[0][0]);
   hamilton.grid.L[1] = static_cast<int>(rL[1][1]);
   hamilton.grid.L[2] = static_cast<int>(rL[2][2]);
   int nsite = hamilton.grid.size();
   if(verbose) std::cout << "ReadStrOut nsite=" << nsite << std::endl;
   // Create a list of species
   std::map<std::string,int> ispecies;
   for(int i=0; i<hamilton.species.size(); i++) 
      ispecies.insert( std::make_pair(hamilton.species[i],i) );
   if(verbose) std::cout << "ReadStrOut MSpecies=" << ispecies.size() << std::endl;
   // Read the configuration proper
   sigma.resize(nsite);
   int isite = 0;
   while( isite<nsite && !fin.eof() )
   {
      float r0, r1, r2;
      std::string label;
      fin >> r0 >> r1 >> r2 >> label;
      int hn0 = static_cast<int>(2*r0);
      int hn1 = static_cast<int>(2*r1);
      int hn2 = static_cast<int>(2*r2);
      Grid::NrmVec hn(hn0,hn1,hn2);
      int index = hamilton.grid.index(hn);
      sigma[index] = ispecies[label]; 
      isite++;
  }
  // if(isite!=nsite) std::cout << "ReadStrOut only " << isite << " atoms listed out of " << nsite << std::endl; 
  if(verbose) std::cout << "ReadStrOut config read" << std::endl;
}


#endif
