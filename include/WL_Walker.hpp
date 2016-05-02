#ifndef WLWalker_HPP_
#define WLWalker_HPP_

#include "WL_Window.hpp"
#include "MC_Walker.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>


struct WL_State
{
   int    ibin;
   double Sval;
   void copy(const WL_State& orig) { ibin=orig.ibin; Sval=orig.Sval; }
   void operator=(const WL_State& orig) { copy(orig); }
};


// Expands an MCWalker to include Wang-Landau information
template<typename Config, typename Observables>
class WL_Walker : public MC_Walker<Config,Observables>
{

public:

   typedef MC_Walker<Config,Observables> BASE;

   double               wlgamma;            // Wang-Landau gamma = ln(f)
   WL_State             wl_now,wl_old;      // Wang-Landau information
   WL_Window            window;             // Wang-Landau Window provides limits and binning
   std::vector<double>  S;                  // Entropy = log-density-of-states
   std::vector<double>  Slast;              // Last set of values, used to check convergence
   std::vector<double>  Sfixed;             // Changes only very slowly
   std::vector<double>  h;                  // Visit histogram

   // The sampling window might be smaller than the energy window
   double WEmin,WEmax;

   // Combined WL-Transition Matrix: RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
   // wleta comes from Eq. (24) S_WL[I] <- (1-wleta)*S_WL[I] + wleta*S_TM[I]
   std::vector<double>  Sittm;              // A local version of the Infinite-temperature transition matrix result
   std::vector<double>  Vittm;              // A local version of the Vittm error matrix
   double    EStepSum;                      // Statistics of how big a MC step is
   long long EStepCnt;



public:

   // Clear the visit histogram
   void clear_histogram() 
   { 
      if( h.size()!=window.NBin ) h.resize(window.NBin); 
      fill(h.begin(),h.end(),0); 
   }

   // Clear the estimated entropy
   void clear_entropy() 
   { 
      S.resize(window.NBin); 
      fill(S.begin(),S.end(),0); 
      Sfixed.resize(window.NBin);
      fill(Sfixed.begin(),Sfixed.end(),0);
   } 

   // Set values of member array. Window must already be set
   void set_values(const std::vector<double>& E, const std::vector<double>& y, std::vector<double>& seq)
   {
      const int nread = E.size();
      if( E.size()<2 ) 
      {
         std::cout << __FILE__ << ":" << __LINE__ << " Energy array is length " << nread << std::endl;
         int nbin = window.NBin;
         seq.resize(nbin);
         for(int i=0; i<nbin; i++) seq[i] = 0;
         return;
      }
      if( y.size()<nread ) std::cout << __FILE__ << ":" << __LINE__ << " array y is too small" << std::endl;
      int iread = 0;
      int nbin = window.NBin;
      seq.resize(nbin);
      for(int ibin=0; ibin<nbin; ibin++)
      {
         double Eset = window.unbin(ibin);
         while( iread<nread && E[iread]<Eset ) iread++;
         if(iread>=nread) iread = nread-1;
         if(iread<1) iread = 1;
         double E0 = E[iread];
         double y0 = y[iread];
         double E1 = E[iread-1];
         double y1 = y[iread-1];
         double slope = (y1-y0)/(E1-E0);
         seq[ibin] = y0 + slope*(Eset-E0);
      }
   }

   void set_fixed(const std::vector<double>& E, const std::vector<double>& y)
   {
      this->set_values(E,y,this->Sfixed);
   }

   // Deep copy of the walker
   void copy(const WL_Walker<Config,Observables>& orig)
   {
      this->BASE::copy(orig);
      wlgamma = orig.wlgamma;
      wl_now = orig.wl_now;
      wl_old = orig.wl_old;
      window = orig.window;
      S = orig.S;
      Sfixed = orig.Sfixed;
      Sittm = orig.Sittm;
      Vittm = orig.Vittm;
      h = orig.h;
      EStepSum = orig.EStepSum;
      EStepCnt = orig.EStepCnt;
   }

   WL_Walker()
   {
   }

   WL_Walker(const WL_Walker& orig)
   {
       this->copy(orig);
   }

   // Use this to find the current bin of the walker
   int bin() { return wl_now.ibin = window.bin(BASE::now.E); }

   // Estimate Entropy/lndos using linear interpolation
   double get_lndos(double E, bool LinearInterp=true)
   {
      int ibin = window.bin(E);
      double E0 = window.unbin(ibin);
      double S0 = S[ibin] + Sfixed[ibin];
      if( !LinearInterp ) return S0;
      double S1,E1;
      if( E<E0 ) 
      {
         // look left
         if( ibin==0 )
         {
            // Interpolate to edge half a bin away. Force WL portion to zero
            E1 = window.Elo;
            S1 = S[ibin]+Sfixed[ibin] - 0.5*( Sfixed[ibin+1]-Sfixed[ibin] ); 
         }
         else
         {
            E1 = window.unbin(ibin-1);
            S1 = S[ibin-1] + Sfixed[ibin-1]; 
         }
      }
      else
      {
         // look right
         if( ibin==(window.NBin-1) )
         {
            E1 = window.Ehi; 
            S1 = S[ibin]+Sfixed[ibin] + 0.5*( Sfixed[ibin]-Sfixed[ibin-1] );
         }
         else
         {
            E1 = window.unbin(ibin+1);
            S1 = S[ibin+1] + Sfixed[ibin+1];
         }
      }
      double slope = (S1-S0)/(E1-E0);
      double SE = S0 + slope*(E-E0);
      return SE;
   }

   float flatness() const
   {
      // Find the prefactor that normalizes histogram such that
      // h_ibin = 1 for perfectly flat histogram (divide by mean h_ibin)
      long long htot = 0;
      int hbin = 0;
      for(int i=0; i<this->window.NBin; i++)
      {
	 if(h[i]>0) 
	 { 
	    htot += h[i]; 
	    hbin++; 
	 } 
      }
      double hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
      // Calculate the flatness for the frozen entropy
      // Q = (1/nbin) \sum_ibin (h_ibin-1)^2
      double Q = 0;
      for(int i=0; i<window.NBin; i++) 
      {
	 if( h[i]>0 )
	 {
	    double Qi=(hnorm*h[i]-1);
 	    Q+=Qi*Qi;
	 }
      }
      if( hbin<10 ) Q=1.e+6;   // Not flat if trapped
      return Q;
   }

/*
   bool move_walls(double qflat=0.1)
   {
      //if( move_walls_pileup() ) return true;
      // Calculate the flatness from high energy down
      int ilft = window.bin(WEmin);
      int irgt = window.bin(WEmax);
      const int mincts = 500;   // Minimum vists before move wall past a bin
      int jmax = irgt;
      long long htot = 0;
      int hbin = 0;
      int jbin=jmax;
      for(int kbin=jmax; kbin>0 && kbin>(jmax-5); kbin--)
         if( h[kbin]>1000000 ) jbin=kbin;
      bool flat = true;
      while( jbin>=0 && (flat || h[jbin]>1000000) )
      {
         if( h[jbin]>mincts )
         {
            htot += h[jbin];
            hbin++;
            double hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
            double Q = 0;
            for(int i=jbin; i<=jmax; i++)
            {
               double Q1 = hnorm*h[i]-1;
               Q += Q1*Q1;
            }
            if(hbin>1) Q /= static_cast<double>(hbin);
            flat = (Q<=qflat);
            jbin--;
         }
         else
         {
            flat = false;
         }
      }
      jmax = jbin;
      // Cacluate the flatness from low energy up
      int jmin = ilft;
      htot = 0;
      hbin = 0; 
      jbin = jmin;
      if( jbin<(window.NBin-1) && h[jbin+1]>1000000 ) jbin++;
      flat = true;
      while( jbin<window.NBin && (flat || h[jbin]>1000000) )
      {
         if( h[jbin]>mincts )
         {
            htot += h[jbin];
            hbin++;
            double hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
            double Q = 0;
            for(int i=jmin; i<=jbin; i++)
            {
               double Q1 = hnorm*h[i]-1;
               Q += Q1*Q1;
            }
            if(hbin>1) Q /= static_cast<double>(hbin);
            flat = (Q<=qflat);
            jbin++;
         }
         else
         {
            flat = false;
         }
      }
      const int overlap = 5;
      jmin = jbin;
      flat = (jmin>jmax);
      jmax += overlap;
      if(jmax<irgt) irgt=jmax;
      jmin -= overlap;
      if(jmin>ilft) ilft=jmin;
      WEmin = window.unbin(ilft);
      WEmax = window.unbin(irgt);
      return flat;
   }
*/

public:

   void save_initial()
   {
      BASE::save_initial();
      wl_old = wl_now;
   }

   void restore_initial()
   {
      wl_now = wl_old;
      BASE::restore_initial();
   }

};


// Return integer index of wlgamma assuming wlgamma_i+1 = wlgamma_i/2
int igamma(double wlgamma)
{
   int p2 = static_cast<int>(-std::log(wlgamma)/std::log(2.)+0.01);
   return p2;
}


#endif  // WLWalker_HPP_
