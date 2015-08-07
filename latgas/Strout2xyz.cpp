// Strout2xyz.cpp -- convert a str.out from ATAT into XYZ format of chemical vizualization
//
// May 7, 2014
// Greg Brown (browngrg@comcast.net,gbrown@fsu.edu)
//


#include "Vec.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>


void ReadStrOut(std::string fname)
{
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
   // 
   std::vector<std::string> atom;
   std::vector<Vec3> npos; 
   while( !fin.eof() )
   {
      float r0, r1, r2;
      std::string label,buff;
      fin >> r0 >> r1 >> r2 >> label;
      std::getline(fin,buff);
      if( !fin.eof() )
      {
         Vec3 pos=r0*a[0]+r1*a[1]+r2*a[2];
         atom.push_back(label);
         npos.push_back(pos);
      }
   }
   int nsite = atom.size();
   std::cout << nsite << std::endl;
   std::cout << "Converted from str.out" << std::endl;
   for(int isite=0; isite<nsite; isite++) 
   {
      std::cout << atom[isite] << " " << npos[isite][0] << " " << npos[isite][1] << " " << npos[isite][2] << std::endl;
   }
}


int main(int argc, char* argv[])
{
   ReadStrOut("str.out");
}
