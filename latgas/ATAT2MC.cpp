// Convert ATAT output to caledonia input
//
// April 24, 2014
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)

// "Multicomponent multisublattice alloys, nonconfigurational entropy
// and other additions to the Alloy Theoretic Automated Toolkit,"
// A. van de Walle   http://arxiv.org/abs/0906.1608v1
// June 8, 2009

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>
#include<map>
#include<string>
#include<unistd.h>

#include"MCModel.hpp"
#include"MCVariables.hpp"
#include"MCFramework.hpp"

int main(int argc, char* argv[])
{
   MCModel model;
   model.ReadATAT();
   if( false ) 
   {
      MCFramework::grid.pbc[0] = false;
      MCFramework::grid.pbc[1] = false;
      MCFramework::grid.pbc[2] = false;
   }
   std::string fname = "str.out";
   if( access( fname.c_str(), F_OK ) != -1 )
   {
      // Get extent, etc from "str.out"
      StoredConfigs<MCModel> configs(model);
      configs.ReadStrOut("str.out");
      if(true)
      {
         int iconfig = configs.index["strout"];
         float E = configs.table[iconfig].macro.E;
         float Eper = E/static_cast<float>(MCFramework::grid.size());
         std::cout << "Energy of str.out is " << E << std::endl;
         std::cout << "Energy per site is " << Eper << std::endl;
      }
   }
   model.Write();
   if(true)
   {
      // Test of read routine
      MCModel model2;
      model2.Read();
      model2.Write("MCModelIn2.txt");
   }
}
