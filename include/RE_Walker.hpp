#ifndef REWalker_HPP_
#define REWalker_HPP_

#include "MC_Walker.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>


// Expands an MCWalker to include Replica Exchange information
template<typename Config, typename Observables>
class RE_Walker : public MC_Walker<Config,Observables>
{

public:

   typedef MC_Walker<Config,Observables> BASE;

   double   kT;            // Temperature of the walker
   int      kTi;           // Index of temperature
   enum     { MAXK=5 };    // Maximum moments of energy
   double   emom[MAXK];    // Accumulated moments of energy
   long     re_swap;       // Number of Replica Exchange swaps

public:

   // Deep copy of the walker
   void copy(const RE_Walker<Config,Observables>& orig)
   {
      this->BASE::copy(orig);
      kT = orig.kT;
      for(int i=0; i<MAXK; i++) emom[i] = orig.emom[i];
      re_swap = orig.re_swap;
   }

   RE_Walker()
   {
       kT = 1;
       for(int i=0; i<MAXK; i++) emom[i]=0;
       re_swap = 0;
   }

   void sample()
   {
      double E = BASE::now.E;
      double Ek = 1;
      for(int k=0; k<MAXK; k++)
      {
         emom[k] += Ek;
         Ek *= E;
      }
   }

   RE_Walker(const RE_Walker& orig)
   {
       this->copy(orig);
   }

   void save_initial()
   {
      BASE::save_initial();
   }

   void restore_initial()
   {
      BASE::restore_initial();
   }

};


#endif  // REWalker_HPP_
