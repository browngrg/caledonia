// GenAlNiCo.cpp -- Generate configurations of this alloy in str.out format
//
// Greg Brown (browngrg@comcast.net,gbrown@fsu.edu)
// July 14, 2014

#include "Grid.hpp"
#include "Vec.hpp"
#include "R1279.hpp"
#include "SeedFromClock.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>


struct AlNiCo
{
public:
   float fAl;
   float fNi;
   float fCo;
   float fCu;
   float fTi;
public:
   AlNiCo(float Al, float Ni, float Co, float Cu=0, float Ti=0)
   {
      fAl = Al;
      fNi = Ni;
      fCo = Co;
      fCu = Cu;
      fTi = Ti;
   }
};


void GenerateAlNiCoTable(std::map<std::string,AlNiCo>& table)
{
   table.clear();
   table["AlNiCo1"]   = AlNiCo(0.12,0.21,0.05,0.03,0.00);
   table["AlNiCo2"]   = AlNiCo(0.10,0.19,0.12,0.03,0.00);
   table["AlNiCo3"]   = AlNiCo(0.12,0.25,0.00,0.03,0.00);
   table["AlNiCo5"]   = AlNiCo(0.08,0.14,0.24,0.03,0.00);
   table["AlNiCo6"]   = AlNiCo(0.08,0.16,0.24,0.03,0.01);
   table["AlNiCo8"]   = AlNiCo(0.07,0.15,0.35,0.04,0.05);
   table["AlNiCo8HC"] = AlNiCo(0.08,0.14,0.38,0.03,0.08);
   table["AlNiCo9"]   = AlNiCo(0.07,0.15,0.35,0.04,0.05);
}



// Populate site vector s with random species of Al, Ni and Fe
void random_alloy(float fAl,  float fNi, float fCo, std::vector<std::string>& s)
{
   float fFe = 1 - fAl - fNi - fCo;
   if( fAl<0 || fAl>1 ) { std::cout << "Fraction fAl=" << fAl << " is not in [0,1]" << std::endl; return; }
   if( fNi<0 || fNi>1 ) { std::cout << "Fraction fNi=" << fNi << " is not in [0,1]" << std::endl; return; }
   if( fCo<0 || fCo>1 ) { std::cout << "Fraction fCo=" << fCo << " is not in [0,1]" << std::endl; return; }
   if( fFe<0 || fFe>1 ) { std::cout << "Fraction fFe=" << fFe << " is not in [0,1]" << std::endl; return; }
   int NSite = s.size();
   int NAl = static_cast<int>(fAl*NSite+0.5);
   int NNi = static_cast<int>(fNi*NSite+0.5);
   int NCo = static_cast<int>(fCo*NSite+0.5);
   if( fNi<0 || fNi>1 ) { std::cout << "Fraction fNi=" << fNi << " is not in [0,1]" << std::endl; return; }
   for(int i=0; i<NSite; i++) s[i] = "Fe";
   R1279 urng;
   urng.seed(SeedFromClock());
   for(int i=0; i<NAl; i++)
   {
      int isite = -1;
      while( isite<0 || isite>=NSite || s[isite].compare("Fe")!=0 )
         isite = static_cast<int>(NSite*urng());
      s[isite] = "Al";
   }
   for(int i=0; i<NNi; i++)
   {
      int isite = -1;
      while( isite<0 || isite>=NSite || s[isite].compare("Fe")!=0 )
         isite = static_cast<int>(NSite*urng());
      s[isite] = "Ni";
   }
   for(int i=0; i<NCo; i++)
   {
      int isite = -1;
      while( isite<0 || isite>=NSite || s[isite].compare("Fe")!=0 )
         isite = static_cast<int>(NSite*urng());
      s[isite] = "Co";
   }
}



int main(int argc, char* argv[])
{

   using namespace std;

   float a[3] = { 2.82, 2.82, 2.82 };

   cout << "How many cubic unit cells (two atoms each)" << endl 
        << "Lx Ly Lz: ";

   int Lx,Ly,Lz;
   cin >> Lx >> Ly >> Lz;

   Grid grid(Lx,Ly,Lz,Sym::CUBIC,Sym::BODYCENTERED);
   int NSite = grid.size();
   vector<string> s(NSite,"Fe");
   cout << "That's " << NSite << " atoms" << endl;

   int istruct;
   cout << "Structure choices:" << endl;
   cout << "   (1) Random Alloy" << endl;
   cout << "What structure? ";
   cin >> istruct;

   float fAl,fNi,fFe;
   switch(istruct)
   {
   case 1: 
      cout << "Fraction Al? "; cin >> fAl;
      cout << "Fraction Ni? "; cin >> fNi;
      fFe = 1. - fAl - fNi;
      cout << "That makes fraction Fe " << fFe;
      random_alloy(fAl,fNi,s);
      break;
   default:
      cout << "Unrecognized choice " << istruct << endl;
      break;
   }

   ofstream fout("GenAlNiFe.out"); 
   fout << setw(10) << a[0] << " " << setw(10) <<   0. << " " << setw(10) <<   0. << endl;
   fout << setw(10) <<   0. << " " << setw(10) << a[1] << " " << setw(10) <<   0. << endl;
   fout << setw(10) <<   0. << " " << setw(10) <<   0. << " " << setw(10) << a[2] << endl;
   fout << setw(10) <<   Lx << " " << setw(10) <<   0  << " " << setw(10) <<   0  << endl;
   fout << setw(10) <<   0  << " " << setw(10) <<   Ly << " " << setw(10) <<   0  << endl;
   fout << setw(10) <<   0  << " " << setw(10) <<   0  << " " << setw(10) <<   Lz << endl;
   for(int isite=0; isite<NSite; isite++)
   {
       Grid::NrmVec rposHN;
       grid.hnorml(isite,rposHN);
       fout << setw(10) << rposHN[0]/2. << " " << setw(10) << rposHN[1]/2. << " " << setw(10) << rposHN[2]/2. << " " << s[isite] << endl;
   }




}
