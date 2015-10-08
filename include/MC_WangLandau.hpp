#ifndef MC_WangLandau_HPP
#define MC_WangLandau_HPP


#include "ITTM.hpp"
#include "WL_Walker.hpp"
#include "WL_Window.hpp"
#include "MPI_Struct.hpp"
#include "MPI_Gang.hpp"
#include "Random.hpp"
#include "Null_Measure.hpp"

#include <cstdlib>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdio.h>
#include <unistd.h>


#ifdef USE_MPI
#include<mpi.h>
#include"MPITypeTraits.hpp"
#endif



class MC_WangLandau
{

public:

   MC_WangLandau();

   MC_WangLandau(const MC_WangLandau& orig) { this->copy(orig); }

   // Create the pool 
   template<typename Walker>
   void init_pool(std::vector<Walker>& pool);

   bool verbose;                              // output level

   enum { latgas, ising, heisenberg } stype;  // Describes the order parameter

public:

   MultiplyWithCarry<4> urng;          // The random number generator 

   template<typename Model, typename Walker>
   void DoWLConverge(Model& model, std::vector<Walker>& wlpool);

   template<typename Model, typename Walker>
   void DoTMConverge(Model& model, std::vector<Walker>& wlpool);

   template<typename Model, typename Walker>
   void DoConverge(Model& model, std::vector<Walker>& wlpool);

   template<typename Model, typename Walker, typename Measure>
   void DoConverge(Model& model, std::vector<Walker>& wlpool, Measure& measure);

   template<typename Model, typename Walker, typename Measure>
   void DoSample(Model& model, std::vector<Walker>& wlpool, Measure& measure);

   template<typename Model, typename Walker>
   void DoWangLandau(Model& model, Walker& wlwalk, long long nstep);

   template<typename Model, typename Walker, typename MeasureObj>
   void DoWangLandau(Model& model, Walker& wlwalk, long long nstep, MeasureObj& measure);

   template<typename Model, typename Walker>
   void DoIntoWindow(Model& model, Walker& wlwalk);

   template<typename Model, typename Walker>
   void DoConstrictorSample(Model& model, std::vector<Walker>& wlwalk);

   double WLAnalyze(const std::vector<double>& h, int ibegin, int iend);

   template<typename Walker>
   double WLAnalyze(const Walker& walker);

   template<typename Walker>
   double TMAnalyze(const Walker& walker);

   template<typename WLWalker>
   double DoAnalyze(std::vector<WLWalker>& walkerpool, WLWalker& average, int iwin);

   template<typename WLWalker>
   double DoAnalyzeMPI(std::vector<WLWalker>& walkerpool, WLWalker& global_walker);

   template<typename Walker>
   void WriteWLDOS(const Walker& walker, int iupdate);

   template<typename Walker>
   void WriteWLBoltzmann(const Walker& walker, int iupdate);

   template<typename Walker>
   void WriteWalkers(const std::vector<Walker>& walker, int iupdate);

   template<typename Walker>
   void CombinedDOS(const std::vector<Walker>& walkerpool, Walker& combined);

   ~MC_WangLandau();

   void copy(const MC_WangLandau& orig);

public:

   // Convergence Information
   long NChunk;                        // Number of steps per test
   long MaxUpdate;                     // Maximum number of convergence iterations
   double Qquit;                       // Convergence factor to stop at
   bool restart_ITTM;                  // Used stored results for ITTM
   double wlgamma_start;
   bool LinearInterp;                  // Use linear interpolation for DOS (good for continuous energy domain)
   double wleta;                       // weighting between WL and ITTM. See Eq. (24) of
                                       // RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
   bool convh;                         // Converge via ln(g_i+1) = ln(g_i) + ln(h_i) when wlgamma = 0;
   bool output_configs;

   double kTlo, kThi;                  // Range of output temperatures for WriteWLBoltzmann

   // Windowing information
   int NWindow;                           // Total number of windows
   int NWalkPerProcess;                   // On each node
   double  Elo,Ehi;
   float   fwinover;                      // fractional overlap of windows
   float   Ebin;                          // A better parameter for binning than nbin
   std::vector<double> Ewin;              // (Elo,Ehi)[0], (Elo,Ehi][1], ... for the windows

   std::vector<int> beginwin;             // Grouping of walkers in serial run
   MPI_Gang         mp_window;            // Details of parallel implementation

   ITTM ittm;                             // measures infinite temperature transition matrix

   long re_iter;                          // replica exchange iteration

   static void positive_temps(std::vector<double>& energy, std::vector<double>& est_lng);
   void partition_windows(std::vector<double>& energy, std::vector<double>& est_lng);
   void partition_windows();

   template<typename Model, typename Walker>
   void DoReplicaExchange(Model& model, std::vector<Walker>& walkerpool);

   template<typename Model, typename Walker>
   void DoReplicaExchange_serial(Model& model, std::vector<Walker>& walkerpool);

   template<typename Model, typename Walker>
   void DoReplicaExchange_serial(Model& model, std::vector<Walker>& walkerpool, int ileft, int imid, int iright);

private:

   enum  { verbose_debug = false };          // very verbose tracing

};


MC_WangLandau::MC_WangLandau()
{
   // Default values
   stype   = latgas;
   verbose = false;
   NWindow = 1;
   NWalkPerProcess = 1;
   NChunk = 1000000;
   MaxUpdate = 100;
   restart_ITTM = false;
   output_configs = false;
   Qquit    = 1.e-4;
   wlgamma_start = 1.;
   wleta    = 0;
   convh    = false;
   re_iter  = 0;
   LinearInterp = true;
   Elo      = +std::numeric_limits<double>::max();         // Elo,Ehi,Ebin,nbin, and fwinover depend on each other
   Ehi      = -std::numeric_limits<double>::max();
   Ebin     = 1.;
   fwinover = 0.5;                                         // fractional overlap, maximal for M/(M+1)
   kTlo     = 0.5;
   kThi     = 5.0;
}


MC_WangLandau::~MC_WangLandau()
{
;
}


void MC_WangLandau::copy(const MC_WangLandau& orig)
{
   stype   = orig.stype;
   verbose = orig.verbose;
   Elo = orig.Elo;
   Ehi = orig.Ehi;
   Ebin = orig.Ebin;
   fwinover = orig.fwinover;
   NWindow = orig.NWindow;
   NWalkPerProcess = orig.NWalkPerProcess;
   NChunk = orig.NChunk;
   MaxUpdate = orig.MaxUpdate;
   restart_ITTM = orig.restart_ITTM;
   Qquit = orig.Qquit;
   wlgamma_start = orig.wlgamma_start;
   wleta = orig.wleta;
   convh = orig.convh;
   beginwin = orig.beginwin;
   Ewin = orig.Ewin;
   mp_window = orig.mp_window;
   re_iter = orig.re_iter;
   LinearInterp = orig.LinearInterp;
   ittm = orig.ittm;
   output_configs = orig.output_configs;
   kTlo = orig.kTlo;
   kThi = orig.kThi;
}


////////////////////////////////////////////////////////////////////////////////
// Take unweighted steps until get into window
template<typename Model, typename Walker>
void MC_WangLandau::DoIntoWindow(Model& model, Walker& walker)
{
   const bool verbose = (this->verbose_debug) && (this->verbose) && (mp_window.pool.iproc==0);
   if(walker.sigma.size()<1)
   {
      int NSite = walker.sigma.size();
      std::cout << __FILE__ << ":" << __LINE__ << " " << "NSite=" << NSite << std::endl;
      exit(1);
   }
   model.calc_observable(walker.sigma,walker.now);
   int ibin = walker.bin();                                   // sets walker.wl_now.ibin also
   long istep = 0;
   while( walker.now.E<walker.WEmin && istep<10000000 )
   {
      istep++;
      if( 0==(istep%100000) )
      {
         std::cout << __FILE__ << ":" << __LINE__
                   << " Walker " << walker.iwalk_global 
                   << " can't find initial window [" << walker.window.Elo << "," << walker.window.Ehi << "] from lo side" << std::endl;
      }
      walker.save_initial();
      model.mc_step(walker,urng);
      if( walker.now.E<walker.old.E )
        walker.restore_initial();
      ibin = walker.bin();                                   // sets walker.wl_now.ibin also
   }
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   while( walker.now.E>walker.WEmax && istep<1000000000 )
   {
      istep++;
      if( 0==(istep%100000) )
      {
         std::cout << __FILE__ << ":" << __LINE__ 
                   << " Walker " << walker.iwalk_global 
                   << " can't find initial window [" << walker.window.Elo << "," << walker.window.Ehi << "] from hi side" << std::endl;
      }
      walker.save_initial();
      model.mc_step(walker,urng);
      if( walker.now.E>walker.old.E )
        walker.restore_initial();
      ibin = walker.bin();                                   // sets walker.wl_now.ibin also
   }
   typename Model::Observables direct;
   model.calc_observable(walker.sigma,direct);
   if( direct!=walker.now )
   {
      std::cout << __FILE__ << ":" << __LINE__ << " Energy does not match direct calculation" << std::endl
                << "Direct: E=" << direct.E << " X=" << direct.X << std::endl
                << "Now:    E=" << walker.now.E << " X=" << walker.now.X << std::endl
                << "diff:   E=" << direct.E-walker.now.E << " X=" << direct.X-walker.now.X << std::endl;
      walker.now = direct;
   }
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// Do NSTEP of the Wang-Landau Method
template<typename Model, typename Walker>
void MC_WangLandau::DoWangLandau(Model& model, Walker& walker, long long nstep)
{
   Null_Measure measure;
   DoWangLandau(model,walker,nstep,measure);
}

template<typename Model, typename Walker, typename MeasureObj>
void MC_WangLandau::DoWangLandau(Model& model, Walker& walker, long long nstep, MeasureObj& measure)
{
   const bool verbose = (this->verbose_debug) && (this->verbose) && (mp_window.pool.iproc==0);   // debugging
   if(walker.sigma.size()<1)
   {
      int NSite= walker.sigma.size();
      std::cout << __FILE__ << ":" << __LINE__ << " iwin=" << walker.window.iwindow << " NSite=" << NSite << std::endl;
      exit(1);
   }
   if( walker.S.size()!=walker.Sfixed.size() )
   {
      std::cout << __FILE__ << ":" << __LINE__ << " Sfixed unset" << std::endl;
      walker.Sfixed.resize( walker.S.size() );
      for(int ibin=0; ibin<walker.Sfixed.size(); ibin++)
         walker.Sfixed[ibin] = 0;
   }
   double wlgamma = walker.wlgamma;
   int ibin = walker.bin();                                            // _C and _PD depend on bin
   walker.wl_now.Sval = walker.get_lndos(walker.now.E,LinearInterp);  
   for(long long istep=0; istep<nstep; istep++)
   {
      // Propose a step in the Markov chain
      // Boundary Effects: p. 1461 of GPB Notebook #9; Email from 12/30/2013 to sullivan@kenyon.edu
      // See also BJ Schultz, K Binder, M Muller, DP Landau, Phys Rev E 67, 067102 (2003)
      // The correct thing to do is to REJECT=!ACCEPT moves outside walls (ignoring is the wrong thing!!!)
      long iskip = 0;
      bool inwall = walker.now.E>=walker.WEmin && walker.now.E<=walker.WEmax;
      bool start_inwall = false;
      while( !start_inwall )           // Need at least one call
      {
         start_inwall = inwall;
         walker.save_initial();
         model.mc_step(walker,urng);
         walker.wl_now.Sval = walker.get_lndos(walker.now.E,LinearInterp); 
         inwall = ( walker.now.E>=walker.WEmin && walker.now.E<=walker.WEmax );
         // Drift towards wall
         if( !start_inwall &&
             ( ( walker.old.E< walker.WEmin && walker.old.E> walker.now.E )
            || ( walker.old.E> walker.WEmax && walker.old.E< walker.now.E ) ) )
         {
            walker.restore_initial();
         }
         // Want old to be inside windows (have to allow for steps outside, though)
         iskip++;
         if(iskip>1000000) 
         {
            iskip = 0;
            std::cout << __FILE__ << ":" << __LINE__ << ": 1000000 tries to step outside of stage window E=" 
                      << walker.now.E << " [" << walker.WEmin << "," << walker.WEmax << "]" << std::endl;
         }
      }
      walker.wl_now.Sval = walker.get_lndos(walker.now.E,LinearInterp);
      ittm.add_sample(walker); 
      // A tolerant way to find acceptance vs rejection
      double deltaS = walker.wl_now.Sval - walker.wl_old.Sval;
      bool accept = inwall && ( (deltaS<=0) || (urng()<exp(-deltaS)) ); 
      if( !accept ) walker.restore_initial();
      walker.imcs++;
      // The Wang-Landau histogram and update
      if( inwall )
      {
         ibin = walker.bin();             // sets walker.wl_now.ibin also
         walker.h[ibin]++;
         walker.wl_now.Sval += wlgamma;   // Need to update both
         walker.S[ibin] += wlgamma;
         measure.add_sample(walker,wlgamma==0);   
      }
   }
   walker.MCSS = static_cast<double>(walker.imcs)/walker.now.V;
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// Analyze the flatness of the last walk
double MC_WangLandau::WLAnalyze(const std::vector<double>& h, int ibegin, int iend)
{
   double Q = 1.e+6;
   // Find the prefactor that normalizes histogram such that
   // h_ibin = 1 for perfectly flat histogram (divide by mean h_ibin)
   long long htot = 0;
   int hbin = 0;
   for(int i=ibegin; i<iend; i++)
   {
      if(h[i]>0) 
      { 
         htot += h[i]; 
         hbin++; 
      } 
   }
   if( hbin<1 ) return Q;
   double hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
   if( false )
   {
      // ignore pileup
      const double thresh = 0.90;
      while( ibegin<iend && std::fabs(hnorm*h[ibegin]-1.)>thresh ) ibegin++;
      while( iend>ibegin && std::fabs(hnorm*h[iend  ]-1.)>thresh ) iend--;
      htot = 0;
      hbin = 0;
      for(int i=ibegin; i<iend; i++)
      {
         if(h[i]>0) 
         { 
            htot += h[i]; 
            hbin++; 
         } 
      }
      hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
      //if( ibegin>1 || iend<(h.size()-2) ) std::cout << "WLQ pileup ibegin=" << ibegin << " iend=" << iend << std::endl;
   }
   // Calculate the flatness for the frozen entropy
   // Q = (1/nbin) \sum_ibin (h_ibin-1)^2
   Q = 0;
   for(int i=ibegin; i<iend; i++) 
   {
      if( h[i]>0 )
      {
         double Qi=(hnorm*h[i]-1);
         Q+=Qi*Qi;
      }
   }
   Q /= static_cast<double>(hbin);
   return Q;
}

template<typename Walker>
double MC_WangLandau::WLAnalyze(const Walker& walker)
{
   return WLAnalyze(walker.h,0,walker.window.NBin);
}

template<typename Walker>
double MC_WangLandau::TMAnalyze(const Walker& walker)
{
   // NOTE: mpic++ on cygwin crashes on vector creation (not necessarily, always here)
   const std::vector<double>& vttt(walker.Vittm);
   double Q = 0;
   if( vttt.size()<3 ) return Q;
   for(int i=0; i<(vttt.size()-3); i++)
   {
      double vabs = std::fabs(vttt[i]);
      if( vabs>0 && vabs>Q ) Q = vabs;
   }
   return Q;
}


////////////////////////////////////////////////////////////////////////////////
// Output Routines
template<typename Walker>
void MC_WangLandau::WriteWLDOS(const Walker& walker, int iupdate)
{
   ittm.write();
   const bool global = ( walker.window.iwindow<0 );
   if( global && mp_window.pool.iproc!=0) return;
   if(mp_window.gang.iproc!=0) return;
   // Find normalization due to number of sites
   double V = walker.now.V;   // Volume of system = NSite
   if(V<1) V=1;
   // Find the normalization for the histogram
   long long htot = 0;
   int hbin = 0;
   for(int i=0; i<walker.window.NBin; i++) { if(walker.h[i]>0) 
      { 
         htot += walker.h[i]; 
         hbin++; 
      } 
   }
   double hnorm = static_cast<double>(hbin)/static_cast<double>(htot);
   // Find the effective entropy Swl + Sfixed
   std::vector<double> Sadd( walker.S.size() );
   for(int i=0; i<walker.window.NBin; i++) Sadd[i] = walker.S[i] + walker.Sfixed[i];
   // Find the transition matrix estimate of S
   std::vector<double> sittm = walker.Sittm;
   std::vector<double> vttt = walker.Vittm;
   // Find normalization for the entropy
   double Smax = -std::numeric_limits<double>::max();
   double S2max = -std::numeric_limits<double>::max();
   double S3max = -std::numeric_limits<double>::max();
   double S4max = -std::numeric_limits<double>::max();
   for(int i=0; i<walker.window.NBin; i++) 
   {
      if(walker.h[i]>0)
      {
         Smax  = std::max(Smax,walker.S[i]);
         S2max = std::max(S2max,walker.Sittm[i]);
         S3max = std::max(S3max,walker.Sfixed[i]);
         S4max = std::max(S4max,Sadd[i]);
      }
   }
   double EStepMean = walker.EStepSum/static_cast<double>(walker.EStepCnt);
   double EStepFrac = EStepMean/walker.window.Ebin;
   bool EbinOK = (EStepFrac>1) && (EStepFrac<10);
   // Write the data
   char fname[100];
   if( walker.wlgamma==0 && iupdate<1 ) iupdate = 1;                            // save at least one WL stage
   if( global )
      sprintf(fname,"WangLandau-%02d.csv",iupdate);                             // Global output
   else
      sprintf(fname,"WangLandau-%02d-%02d.csv",iupdate,walker.window.iwindow);  // By window output
   std::ofstream out(fname);
   int icol = 1;
   out << "# Wang-Landau Density of States" << std::endl; 
   out << "# V = NSite = " << V << std::endl;
   out << "# Wang-Landau gamma = " << walker.wlgamma << std::endl;
   out << "# iupdate = " << iupdate << std::endl;
   out << "# Average Energy step = " << EStepMean << " = " << EStepFrac << " of Ebin" << std::endl;
   out << "# Column " << icol++ << ": E central energy of bin" << std::endl;
   out << "# Column " << icol++ << ": S = Entropy = lng" << std::endl;
   out << "# Column " << icol++ << ": s = S/V = entropy per spin" << std::endl;
   out << "# Column " << icol++ << ": h = visit histogram raw" << std::endl;
   out << "# Column " << icol++ << ": h normalized to 1 for flat histogram" << std::endl;
   out << "# Column " << icol++ << ": ibin" << std::endl;
   out << "# Column " << icol++ << ": SITTM = transition matrix estimate of S (ok=" << EbinOK << ")"  << std::endl;
   out << "# Column " << icol++ << ": v = violation of total transition matrix (measure of sampling)" << std::endl;
   out << "# Column " << icol++ << ": Sfixed" << std::endl;
   out << "# Column " << icol++ << ": Swanglandau" << std::endl;
   out.precision(16);
   for(int ibin=0; ibin<walker.window.NBin; ibin++)
   {
      if( walker.h[ibin]>0 || sittm[ibin]!=0 )
      {
         double Swl = walker.S[ibin] - S4max;
         double Sfx = walker.Sfixed[ibin] - S3max;
         double Stm = sittm[ibin] - S2max;
         double S = Sadd[ibin] - Smax;
         double h = hnorm*walker.h[ibin];
         out << walker.window.unbin(ibin) << " " << S << " " << S/V 
             << " " << walker.h[ibin] << " " << h << " " << ibin 
             << " " << Stm << " " << vttt[ibin] 
             << " " << Sfx << " " << Swl << std::endl;
      }
   }
}

template<typename Walker>
void MC_WangLandau::WriteWLBoltzmann(const Walker& walker, int iupdate)
{
   const bool global = ( walker.window.iwindow<0 );
   if( !global ) return;
   if( mp_window.pool.iproc!=0) return;
   // Write the data
   char fname[100];
   if( walker.wlgamma==0 && iupdate<1 ) iupdate = 1;                          // save at least one WL stage
   sprintf(fname,"WLBoltzmann-%02d.csv",iupdate);                             // Global output
   std::ofstream fout(fname);
   fout << "# Wang-Landau thermodynamic properties calculated via Boltzmann sum" << std::endl;
   fout << "# V = NSite = " << walker.now.V << std::endl;
   fout << "# Wang-Landau gamma = " << walker.wlgamma << std::endl;
   fout << "# iupdate = " << iupdate << std::endl;
   fout << "# Column 1: Temperature in energy units of original data" << std::endl;
   fout << "# Column 2: Specific heat from <E^2>-<E>^2" << std::endl;
   fout << "# Column 3: Entropy S(E)" << std::endl;
   fout << "# Column 4: microcanonical energy E" << std::endl;
   fout << "# Column 5: Helmholtz free energy A = E-kT*S" << std::endl;
   fout << "# Column 6: beta = 1/kT" << std::endl;
   fout << "# Column 7: gamma = -(kT)^2*C" << std::endl;
   fout.precision(15);
   const int nKTpt = 500;
   const double dkT = (kThi-kTlo)/static_cast<double>(nKTpt-1);
   const int npt = walker.S.size();
   if(npt<3) return;
   std::vector<double> E(npt);
   for(int i=0; i<npt; i++) E[i] = walker.window.unbin(i);
   std::vector<double> lng(npt);
   for(int i=0; i<npt; i++) lng[i] = walker.S[i] + walker.Sfixed[i];
   std::vector<double> A(npt);
   for(int ikt=0; ikt<nKTpt; ikt++)
   {
      double kT = kTlo + ikt*dkT;
      for(int i=0; i<npt; i++) A[i] = E[i] - kT*lng[i]; 
      double Amin = std::numeric_limits<double>::max();
      for(int i=0; i<npt; i++) Amin = std::min(Amin,A[i]);
      for(int i=0; i<npt; i++) A[i] -= Amin;
      double Z = 0;
      double E1 = 0;
      double E2 = 0;
      double A1 = 0;
      double beta = 1./kT;
      for(int i=0; i<npt; i++)
      {
         double boltz = std::exp(-A[i]/kT);
         Z  += boltz;
         E1 += boltz*E[i];
         E2 += boltz*E[i]*E[i];
         A1 += boltz*A[i];
      }
      E1 /= Z;
      E2 /= Z;
      A1 /= Z;
      double C = E2-E1*E1;
      double gamma = -C*kT*kT;
      fout << std::setw(20) << kT << " "
           << std::setw(20) << C << " "
           << std::setw(20) << lng[ikt] << " "
           << std::setw(20) << E1 << " "
           << std::setw(20) << A1 << " "
           << std::setw(20) << beta << " "
           << std::setw(20) << gamma << std::endl;
   }
}


template<typename Walker>
void MC_WangLandau::WriteWalkers(const std::vector<Walker>& walkerpool, int iupdate)
{
   if( mp_window.pool.iproc>10 && mp_window.gang.iproc!=0 ) return;
   char fname[100];
   sprintf(fname,"MC_WangLandauConfig-%02d-%02d.txt",mp_window.pool.iproc,iupdate);
   std::ofstream fout(fname);
   fout << "# Configurations from processor " << mp_window.pool.iproc << " at update " << iupdate << std::endl;
   fout << "# Energy followed by values of microconfiguration sigma" << std::endl;
   for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
   {
      fout << walkerpool[iwalk].now.E;
      for(int i=0; i<walkerpool[iwalk].sigma.size(); i++)
         fout << " " << walkerpool[iwalk].sigma[i];
      fout << std::endl;
   }
}



////////////////////////////////////////////////////////////////////////////////
// Use moving walls to force walkers throughout window
// Currently walkers can get trapped outside wall and not find their way in
template<typename Model, typename Walker>
void MC_WangLandau::DoConstrictorSample(Model& model, std::vector<Walker>& walkerpool)
{
   const int NWalker = walkerpool.size();
   const int NWindow_Local = beginwin.size()-1;
   // Clear histogram
   for(int iwalk=0; iwalk<NWalker; iwalk++) 
   {
      Walker& walker(walkerpool[iwalk]);
      if( walker.h.size()!=walker.S.size() ) { walker.h.resize(walker.S.size()); }
      std::fill(walker.h.begin(),walker.h.end(),0); 
      std::fill(walker.S.begin(),walker.S.end(),0);            // Enforce S==0 outside stage window
      walker.wlgamma = 1;
   }
   long NChunkDo = NChunk*walkerpool[0].now.V;
   // Do the actual convergence
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   std::vector<bool> converged(NWindow_Local,false);
   bool all_converged = true;
   for(int i=0; all_converged && i<converged.size(); i++) all_converged = all_converged && converged[i];
   long long totalstep = 0;
   while( !all_converged )
   {
      for(int iwin=0; iwin<NWindow_Local; iwin++) 
      {
         if( !converged[iwin] )
         {
            for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++) 
            {
               DoWangLandau(model,walkerpool[iwalk],NChunkDo);                  // Do the steps using Wang-Landau algorithm
               totalstep += NChunkDo;
            }
         }
         // adjust walls
         Walker average;
         average.copy(walkerpool[beginwin[iwin]]);
         average.clear_histogram();
         for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++)
            for(int j=0; j<average.h.size(); j++) average.h[j] += walkerpool[iwalk].h[j];
         average.window = walkerpool[beginwin[iwin]].window;
         average.WEmin = walkerpool[beginwin[iwin]].WEmin;
         average.WEmax = walkerpool[beginwin[iwin]].WEmax;
#        ifdef USE_MPI
         if( mp_window.ngang>1 )
         {
            std::vector<double> mpi_buffer_d(average.h.size(),0);
            if(mp_window.gang.in()) PI_Allreduce(&(average.h[0]),&(mpi_buffer_d[0]),average.h.size(),MPI_DOUBLE,MPI_SUM,mp_window.gang.comm);
            for(int j=0; j<average.h.size(); j++) average.h[j] = mpi_buffer_d[j];
         }
#        endif
         const int mincts = 50000;
         int ilft = average.window.bin(average.WEmin);
         int irgt = average.window.bin(average.WEmax);
         while( ilft<irgt && average.h[ilft]>mincts ) ilft++;
         while( irgt>ilft && average.h[irgt]>mincts ) irgt--;
         if( irgt<=ilft )
         {
            converged[iwin] = true;
         }
         else
         {
            average.WEmin = average.window.unbin(ilft);
            average.WEmax = average.window.unbin(irgt);
            for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++)
            {
               // all walkers get the same window
               walkerpool[iwalk].WEmin = average.WEmin;                              
               walkerpool[iwalk].WEmax = average.WEmax;
            }
         }
         if( mp_window.gang.iproc==0 ) std::cout << "Constrictor iproc=" << mp_window.pool.iproc << " iwindow=" << iwin << " ilft=" << ilft << " irgt=" << irgt << " converge=" << converged[iwin] << std::endl;
      }
      all_converged = true;
      for(int i=0; all_converged && i<converged.size(); i++) all_converged = all_converged && converged[i];
#     ifdef USE_MPI
      bool local_converged = all_converged;
      if(mp_window.pool.in()) MPI_Allreduce(&local_converged,&all_converged,1,MPI_C_BOOL,MPI_LAND,mp_window.pool.comm);
#     endif
   }
   Walker global_walker;
   double WLQ = DoAnalyzeMPI(walkerpool,global_walker);
   global_walker.wlgamma = walkerpool[0].wlgamma;
   int iupdate = 0;
   WriteWLDOS(global_walker,iupdate);
   WriteWLBoltzman(global_walker,iupdate);
   for(int iwalk=0; iwalk<NWalker; iwalk++) 
   {
      std::fill(walkerpool[iwalk].h.begin(),walkerpool[iwalk].h.end(),0);
      std::fill(walkerpool[iwalk].S.begin(),walkerpool[iwalk].S.end(),0);
      walkerpool[iwalk].wlgamma = 0;
      walkerpool[iwalk].WEmin = walkerpool[iwalk].window.Elo;
      walkerpool[iwalk].WEmax = walkerpool[iwalk].window.Ehi;
      int iwin = walkerpool[iwalk].window.iwindow;
      int istart = global_walker.window.bin(Ewin[1*iwin+0]);
      for(int ibin=0; ibin<walkerpool[iwalk].S.size(); ibin++)
         walkerpool[iwalk].Sfixed[ibin] = global_walker.Sittm[istart+ibin];
   }
}

////////////////////////////////////////////////////////////////////////////////
// Converge the DOS/Entropy to self-consistency by decreasing wlgamma
template<typename Model, typename WLWalker>
void MC_WangLandau::DoWLConverge(Model& model, std::vector<WLWalker>& walkerpool)
{
   wlgamma_start = 1;
   wleta = 0;
   DoConverge(model,walkerpool);
}

////////////////////////////////////////////////////////////////////////////////
// Converge the DOS/Entropy using the Infinite Temperature Transiton Matrix 
// as the DOS result, but a Wang-Landau approach for sampling energy space
template<typename Model, typename WLWalker>
void MC_WangLandau::DoTMConverge(Model& model, std::vector<WLWalker>& walkerpool)
{
#  if 1
   wleta = 1;
   wlgamma_start = 1;
#  else
   DoConstrictorSample(model,walkerpool);
   wleta = 1;
   wlgamma_start = 0;
#  endif
   DoConverge(model,walkerpool);
}


////////////////////////////////////////////////////////////////////////////////
// Converge the DOS
template<typename Model, typename WLWalker>
void MC_WangLandau::DoConverge(Model& model, std::vector<WLWalker>& walkerpool)
{
   Null_Measure measure;
   this->DoConverge(model,walkerpool,measure);
}

////////////////////////////////////////////////////////////////////////////////
// Converge the DOS
template<typename Model, typename WLWalker, typename Measure>
void MC_WangLandau::DoConverge(Model& model, std::vector<WLWalker>& walkerpool, Measure& measure)
{
   //bool verbose = (this->verbose_debug) && (this->verbose) && (mp_window.pool.iproc==0);
   const int NWalker = walkerpool.size();
   long NChunkDo = NChunk*walkerpool[0].now.V;
   if(verbose) std::cout << "wlgamma_start = " << wlgamma_start << std::endl;
   for(int iwalk=0; iwalk<NWalker; iwalk++)
   {
      walkerpool[iwalk].imcs = 0;
      walkerpool[iwalk].MCSS = 0;
      walkerpool[iwalk].wlgamma = wlgamma_start;
      walkerpool[iwalk].clear_histogram();
      for(int ibin=0; ibin<walkerpool[iwalk].S.size(); ibin++) 
      {
         walkerpool[iwalk].S[ibin] = 0;
         walkerpool[iwalk].h[ibin] = 0;
      }
      if( walkerpool[iwalk].Sfixed.size()==0 ) 
         walkerpool[iwalk].Sfixed = walkerpool[iwalk].S;  // zeros from above
   }
   //std::cout << model.header();
   double qthresh = 1;
   double WLQold = 0;
   double fvisit_old = 0;
   double fvisit_target = 1;
   int    fvisit_same = 0;
   unsigned long int istep = 0;
   unsigned long int iloop = 0;
   unsigned long int iupdate = 0;
   WLWalker global_walker;
   while( iupdate<MaxUpdate )
   {
      if(verbose && iupdate>5) 
      {
         verbose = false;
         std::cout << __FILE__ << ":" << __LINE__ << " iupdate=" << iupdate << ". Going silent on DoConverge." << std::endl;
      }
      // Sampling
      for(int iwalk=0; iwalk<NWalker; iwalk++)
         DoWangLandau(model,walkerpool[iwalk],NChunkDo,measure);
      istep += NChunkDo;
      // Analysis
      measure.write();
      double WLQ = DoAnalyzeMPI(walkerpool,global_walker);
      global_walker.wlgamma = walkerpool[0].wlgamma;
      WriteWLDOS(global_walker,iupdate);
      WriteWLBoltzmann(global_walker,iupdate);
      if( output_configs ) WriteWalkers(walkerpool,iupdate);
      if( mp_window.ngang>1 ) DoReplicaExchange(model,walkerpool);
      // Decide about advancing WL parameter
      int ivisit = 0;
      for(int ibin=0; ibin<global_walker.h.size(); ibin++)
         ivisit += (global_walker.h[ibin]>0);
      double fvisit = static_cast<double>(ivisit)/static_cast<double>(global_walker.h.size());
      bool allvisit = (fvisit>=fvisit_target);
      if ( fvisit==fvisit_old )
      {
        bool started = allvisit;
        fvisit_same++;
        if( fvisit_same>=10   && fvisit>0.99 ) allvisit = true;  // Allow for some empty bins
        if( fvisit_same>=50   && fvisit>0.95 ) allvisit = true;  // Allow for some empty bins
        if( fvisit_same>=1000 && fvisit>0.90 ) allvisit = true;   // Allow for some empty bins
        if( global_walker.wlgamma==1 && fvisit_same>=20  ) allvisit = true;  // Allow for many empty
        if( !started && allvisit ) 
        {
           fvisit_target = fvisit;
           if( mp_window.pool.iproc==0 ) std::cout << "New target for fraction of bins visited = " << fvisit_target << std::endl;
        }
      }
      else
      {
         fvisit_same = 0;
      }
      bool qconverged = (std::fabs(WLQ-WLQold)<qthresh);
      if( mp_window.pool.iproc==0 )
      {
         double psecs = static_cast<double>(clock())/static_cast<double>(CLOCKS_PER_SEC);
         std::cout << "loop= " << iloop << " level=" << iupdate << " WLgamma=" << walkerpool[0].wlgamma << " WLQ=" << WLQ << " istep=" << istep << " visit=" << fvisit << " time = " << psecs << " secs" << std::endl;
      }
      // Work on convergence
      if( global_walker.wlgamma>0 && WLQ<Qquit && allvisit   )
      {
         for(int iwalk=0; iwalk<NWalker; iwalk++)
         {
            int iwin = walkerpool[iwalk].window.iwindow;
            int istart = global_walker.window.bin(Ewin[1*iwin+0]);
            for(int ibin=0; ibin<walkerpool[iwalk].S.size(); ibin++)
               walkerpool[iwalk].S[ibin]  = (1.-wleta)*global_walker.S[istart+ibin]+wleta*(global_walker.Sittm[istart+ibin]-walkerpool[iwalk].Sfixed[ibin]);
            walkerpool[iwalk].wlgamma /= 2;
         }
         fvisit = 0;
         WLQold = 0;
         for(int iwalk=0; iwalk<NWalker; iwalk++)
               std::fill(walkerpool[iwalk].h.begin(),walkerpool[iwalk].h.end(),0);
         istep = 0;
         iupdate++;
      }
      if( global_walker.wlgamma==0 && this->convh && qconverged && allvisit )
      {
         double global_min = 1e-20;
         for(int ibin=0; ibin<global_walker.h.size(); ibin++) 
            if( global_walker.h[ibin]>0 && global_walker.h[ibin]<global_min ) global_min = global_walker.h[ibin];
         global_min = std::log(global_min);
         double ittm_min = 0;
         for(int ibin=0; ibin<global_walker.h.size(); ibin++) 
            if( global_walker.h[ibin]>0 && global_walker.Sittm[ibin]<ittm_min ) ittm_min = global_walker.Sittm[ibin];
         for(int iwalk=0; iwalk<NWalker; iwalk++)
         {
            int iwin = walkerpool[iwalk].window.iwindow;
            int istart = global_walker.window.bin(Ewin[1*iwin+0]);
            for(int ibin=0; ibin<walkerpool[iwalk].S.size(); ibin++)
            {
               double lng = global_min;
               if( global_walker.h[istart+ibin]>0 ) lng = std::log(global_walker.h[istart+ibin]);
               double ittm = walkerpool[iwalk].S[ibin];
               if( global_walker.Sittm[istart+ibin]!=0 ) ittm = global_walker.Sittm[istart+ibin]-walkerpool[iwalk].Sfixed[ibin];
               walkerpool[iwalk].S[ibin]  = (1.-wleta)*(walkerpool[iwalk].S[ibin]+lng) + wleta*ittm;
            }
         }
         if( allvisit )
         {  // Only call it an update if visiting all bins
            iupdate++;
            istep = 0;
            qthresh *= 0.1;
         }
         fvisit = 0;
         WLQold = 0;
         for(int iwalk=0; iwalk<NWalker; iwalk++)
               std::fill(walkerpool[iwalk].h.begin(),walkerpool[iwalk].h.end(),0);
      }
      WLQold = WLQ;
      fvisit_old = fvisit;
      iloop++;
   }
   // Copy global walker histogram into walks
   for(int iwalk=0; iwalk<NWalker; iwalk++)
   {
      int iwin = walkerpool[iwalk].window.iwindow;
      int istart = global_walker.window.bin(Ewin[1*iwin+0]);
      for(int ibin=0; ibin<walkerpool[iwalk].h.size(); ibin++)
         walkerpool[iwalk].h[ibin] = global_walker.h[istart+ibin];
   }
}


////////////////////////////////////////////////////////////////////////////////
// Flat Histogram sampling using given DOS
template<typename Model, typename Walker, typename Measure>
void MC_WangLandau::DoSample(Model& model, std::vector<Walker>& walkerpool, Measure& measure)
{
   wlgamma_start = 0;
   wleta = 0;
   return this->DoConverge(model,walkerpool,measure);
}


// Upon compeletion, the walker "average" will have Sfixed used during the sampling,
// S the new ITTM estimate of the DOS, h the total visit histogram, and the total
// transition matrices
template<typename WLWalker>
double MC_WangLandau::DoAnalyze(std::vector<WLWalker>& walkerpool, WLWalker& average, int iwin)
{
   // Analyze flatness individually
   double WLQ = -1;
   for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++)
      WLQ = std::max(WLQ,WLAnalyze(walkerpool[iwalk]));                        // Normalized flatness from UGA13 Eq (14)
#  ifdef USE_MPI
   double WLQ_local = WLQ; 
   MPI_Allreduce(&WLQ_local,&WLQ,1,MPI_DOUBLE,MPI_MAX,mp_window.gang.comm);
#  endif
   average.copy( walkerpool[beginwin[iwin]] );
   average.EStepSum = 0;
   average.EStepCnt = 0;
   average.clear_histogram();
   average.clear_entropy();
   average.Sfixed = walkerpool[beginwin[iwin]].Sfixed;
   for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++)
   {
      for(int j=0; j<average.S.size(); j++) average.S[j] += walkerpool[iwalk].S[j];
      for(int j=0; j<average.h.size(); j++) average.h[j] += walkerpool[iwalk].h[j];
   }
#  ifdef USE_MPI
   // Reduce across MPI processes
   if( mp_window.ngang>1 )
   {
      // Should not reduce h(E) for regular WL, where need to test for flatness before averaging
      std::vector<double> mpi_buffer_d(average.h.size(),0);
      MPI_Allreduce(&(average.h[0]),&(mpi_buffer_d[0]),average.h.size(),MPI_DOUBLE,MPI_SUM,mp_window.gang.comm);
      for(int j=0; j<average.h.size(); j++) average.h[j] = mpi_buffer_d[j];
      mpi_buffer_d.resize(average.S.size(),0);
      MPI_Allreduce(&(average.S[0]),&(mpi_buffer_d[0]),average.S.size(),MPI_DOUBLE,MPI_SUM,mp_window.gang.comm);
      for(int j=0; j<average.S.size(); j++) average.S[j] = mpi_buffer_d[j];
   }
#  endif
   double NWT = static_cast<double>(mp_window.gang.nproc*NWalkPerProcess); 
   if( mp_window.ngang==1 ) NWT = beginwin[iwin+1] - beginwin[iwin];
   for(int j=0; j<average.S.size(); j++) average.S[j] /= NWT;
   long long htot = 0; for(int i=0; i<average.h.size(); i++) htot += average.h[i];
   double Q = TMAnalyze(average);                                                                     // TTT Violation
   double RMSQ = std::sqrt(static_cast<double>(average.h.size())/static_cast<double>(htot));          // sqrt(N) fluctuations
   average.Sfixed = walkerpool[beginwin[iwin]].Sfixed;
   return WLQ;
}  

// Outer layer of analysis and reduction
template<typename WLWalker>
double MC_WangLandau::DoAnalyzeMPI(std::vector<WLWalker>& walkerpool, WLWalker& global_walker)
{
   double WLQ_max = 0;
   std::vector<WLWalker> averagepool(NWindow,walkerpool[0]);        // Copy to get arrays right size. Needs all windows same NBin :(
   const int NWindow_Local = beginwin.size()-1;
   for(int iwin=0; iwin<NWindow_Local; iwin++)
   {
      WLWalker average;
      double WLQ = DoAnalyze(walkerpool,average,iwin);              // Returns maximum WLQ from individual walkers
      WLQ_max = std::max(WLQ,WLQ_max);
      int iwindow = average.window.iwindow;
      averagepool[ average.window.iwindow ].copy(average);
   }
#  ifdef USE_MPI
   if( mp_window.ngang>1 )
   {
      size_t MBin = walkerpool[0].S.size();
      for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
         if(MBin!=walkerpool[iwalk].S.size()) std::cout << __FILE__ << ":" << __LINE__ << " Bad MBin" << std::endl;;
      int buff_size = 5*MBin+1;                                            // The plus one allows passing number of bins
      int tmp = buff_size;
      if(mp_window.pool.in()) MPI_Allreduce(&tmp,&buff_size,1,MPI_INT,MPI_MAX,mp_window.pool.comm);
      std::vector<double> mpi_buffer_d(buff_size);
      // synchronize averagepool across all processes
      // Side effect of DoAnalyze is S=<S> and h=sum(h)
      for(int iwin=0; /*mp_window.ngang>1 &&*/ iwin<NWindow; iwin++)
      {
         // Assumes all windows have same NBin. 
         int root = iwin*mp_window.gang.nproc;                                                            // broadcaster is first process of the window-gang
         int NBin = averagepool[iwin].S.size();
         // Can we avoid transmitting these?
         averagepool[iwin].Sittm.resize(NBin);                                                            // may not have been set yet
         averagepool[iwin].Vittm.resize(NBin); 
         for(int ibin=0;ibin<NBin; ibin++) mpi_buffer_d[0*NBin+ibin] = averagepool[iwin].S[ibin];
         for(int ibin=0;ibin<NBin; ibin++) mpi_buffer_d[1*NBin+ibin] = averagepool[iwin].Sittm[ibin];
         for(int ibin=0;ibin<NBin; ibin++) mpi_buffer_d[2*NBin+ibin] = averagepool[iwin].Vittm[ibin];
         for(int ibin=0;ibin<NBin; ibin++) mpi_buffer_d[3*NBin+ibin] = averagepool[iwin].Sfixed[ibin];
         for(int ibin=0;ibin<NBin; ibin++) mpi_buffer_d[4*NBin+ibin] = averagepool[iwin].h[ibin];
         if(mp_window.pool.in()) MPI_Bcast(&(mpi_buffer_d[0]),buff_size,MPI_DOUBLE,root,mp_window.pool.comm);
         for(int ibin=0;ibin<NBin; ibin++) averagepool[iwin].S[ibin]      = mpi_buffer_d[0*NBin+ibin];
         for(int ibin=0;ibin<NBin; ibin++) averagepool[iwin].Sittm[ibin]  = mpi_buffer_d[1*NBin+ibin];
         for(int ibin=0;ibin<NBin; ibin++) averagepool[iwin].Vittm[ibin]  = mpi_buffer_d[2*NBin+ibin];
         for(int ibin=0;ibin<NBin; ibin++) averagepool[iwin].Sfixed[ibin] = mpi_buffer_d[3*NBin+ibin];
         for(int ibin=0;ibin<NBin; ibin++) averagepool[iwin].h[ibin]      = mpi_buffer_d[4*NBin+ibin];
      }
      double WLQ = WLQ_max;
      if(mp_window.pool.in()) MPI_Allreduce(&WLQ,&WLQ_max,1,MPI_DOUBLE,MPI_MAX,mp_window.pool.comm);
   }
#  endif
   global_walker.sigma = walkerpool[0].sigma;  // Set NSite
   global_walker.now.V = walkerpool[0].now.V;
   CombinedDOS(averagepool,global_walker);
   ITTM global_ittm;
   ittm.get_global(global_ittm);
   global_ittm.calc_SITTM(global_walker.Sittm);
   global_ittm.calc_VITTM(global_walker.Vittm);
   global_walker.EStepSum = global_ittm.EStepSum;
   global_walker.EStepCnt = global_ittm.EStepCnt;
   return WLQ_max;
}


////////////////////////////////////////////////////////////////////////////////
// Combine the DOS across windows and processes
////////////////////////////////////////////////////////////////////////////////

struct MCWL_WinRec 
{ 
public: 
   std::vector<double> E; 
   std::vector<double> lng; 
   std::vector<double> slope;
public:
   MCWL_WinRec(const std::vector<double>& Eval, const std::vector<double>& lngval) { E=Eval; lng=lngval; find_slope(); }
   void find_slope()
   {
      int npts = E.size();
      slope.resize(npts);
      for(int i=1; i<npts; i++)
         slope[i] = (lng[i]-lng[i-1])/(E[i]-E[i-1]);
      slope[0] = slope[1];
   }
   void smooth_slope(int npt=3)
   {
      std::vector<double> smooth(slope.size());
      for(int i=npt; i<(slope.size()-npt); i++)
      {
         smooth[i] = 0;
         for(int j=i-npt; j<i+npt; j++)
            smooth[i] += slope[j];
         smooth[i] /= static_cast<double>(2*npt+1);
      }
      for(int i=npt; i<(slope.size()-npt); i++)
         slope[i] = smooth[i]; 
   }
   int find_energy(double Eval, int istart=0)
   {
      while( istart<E.size() && E[istart]<Eval ) istart++; 
      return istart;
   }
};


template<typename Walker>
void MC_WangLandau::CombinedDOS(const std::vector<Walker>& walkerpool, Walker& combined)
{
   // Set the global window. iwin<0 indicates global window by convention
   combined.window.set_delta(Elo,Ehi,-1,Ebin); 
   int NBinC = combined.window.NBin;
   combined.S.resize(NBinC);
   combined.Sfixed.resize(NBinC);
   combined.h.resize(NBinC);
   int NWin = walkerpool.size();
// for(int iwin=0; iwin<NWin; iwin++) 
// { 
//    combined.EStepSum += walkerpool[iwin].EStepSum; 
//    combined.EStepCnt += walkerpool[iwin].EStepCnt; 
// }
   // Do a simple obliteration approach (for debugging?)
   int ifar = 0;
   for(int iwin=0; iwin<NWin; iwin++)
   {
      int NBinW = walkerpool[iwin].window.NBin;
      int istart = combined.window.bin(Ewin[2*iwin+0]);
      double offS = 0;
      int imatch = 0;
      // find match point
      if( iwin>0 )
      {
         int NGuard = static_cast<double>(NBinW)*fwinover/10.;
         imatch = NGuard;
         int imatch_max = ifar - istart - 1;
         double slope_lo = combined.S[istart+imatch+1] - combined.S[istart+imatch];
         double slope_hi = walkerpool[iwin].S[imatch+1] - walkerpool[iwin].S[imatch];
         bool foo = !(slope_hi>slope_lo);
         while( imatch<imatch_max && (foo^(slope_hi>slope_lo)) )
         {
            slope_lo = combined.S[istart+imatch+1] - combined.S[istart+imatch];
            slope_hi = walkerpool[iwin].S[imatch+1] - walkerpool[iwin].S[imatch];
            imatch++;
         }
         offS = combined.S[istart+imatch] - walkerpool[iwin].S[imatch];
      }
      for(int ibin=imatch; ibin<NBinW; ibin++)
      {
         int ibinc = istart + ibin;
         combined.S[ibinc]      = walkerpool[iwin].S[ibin] + offS;
         combined.Sfixed[ibinc] = walkerpool[iwin].Sfixed[ibin];
         combined.h[ibinc]      = walkerpool[iwin].h[ibin];
      }
      ifar = istart+NBinW;
   }
   double mval = *(std::max_element(combined.S.begin(),combined.S.end()));
   for(int ibin=0; ibin<NBinC; ibin++) combined.S[ibin] -= mval;
   mval = *(std::max_element(combined.Sfixed.begin(),combined.Sfixed.end()));
   for(int ibin=0; ibin<NBinC; ibin++) combined.Sfixed[ibin] -= mval;
}

void MC_WangLandau::positive_temps(std::vector<double>& energy, std::vector<double>& est_lng)
{
   int nel = energy.size();
   int imax = 0;
   double lng_max = est_lng[imax];
   for(int i=0; i<nel; i++)
      if( est_lng[i]>lng_max ) { imax=i; lng_max=est_lng[i]; }
   if( imax>1 )
   {
      energy.resize(imax);
      est_lng.resize(imax);
   }
}


// Equal width windows
void MC_WangLandau::partition_windows()
{
   if( !mp_window.pool.in() ) return;
   int iproc= mp_window.pool.iproc;
   if( iproc==0 ) 
      std::cout << __FILE__ << ":" << __LINE__ <<  "(" << iproc << ") Global Window = [" << Elo << "," << Ehi << "] Ebin=" << Ebin << std::endl;
   if(Elo>Ehi) std::swap(Elo,Ehi);
   float deltawin = (Ehi-Elo)/(static_cast<float>(NWindow)*(1-fwinover)+fwinover);
   int nbin = static_cast<int>( (deltawin/Ebin) + 0.9 );
   Ehi = std::max(Ehi,Elo+Ebin*nbin);
   deltawin = nbin*Ebin;
   if( nbin<10 ) std::cout << __FILE__ << ":" << __LINE__ << " Too few bins per window. nbin=" << nbin << std::endl;
   Ewin.resize(2*NWindow);
   for(int iwin=0; iwin<NWindow; iwin++)
   {
      Ewin[2*iwin+0] = Elo + iwin*deltawin*(1-fwinover);   
      Ewin[2*iwin+1] = Ewin[2*iwin]+deltawin;
   }
   Ewin[2*(NWindow-1)+0] = Ehi - deltawin;                                 // Make last window stretch
   Ewin[2*(NWindow-1)+1] = Ehi;                                            // From end back
   Ewin[0] = Elo;                                                          // When Nwindow==1
   if(iproc==0)
   {
      std::cout << __FILE__ << ":" << __LINE__ << " Windows Ebin=" << Ebin << " numbin=" << nbin << std::endl;
      for(int iwin=0; iwin<NWindow; iwin++) std::cout << iwin << " " << Ewin[2*iwin] << " " << Ewin[2*iwin+1] << std::endl;
   }
}


// Equal delta(lng) windows
void MC_WangLandau::partition_windows(std::vector<double>& energy, std::vector<double>& est_lng)
{
   // Ebin is a better parameter to specify since the optimal binsize is related to the 
   // size of the Monte Carlo step. The number of bins should grow with system size, not the
   // bin size. However, the assumption that nbin is constant across windows is deeply 
   // rooted in the code. That means all the windows need to be the same width until that
   // assumption is broken. Until then, just use the simple window paritioning.
   // MAY HAVE ALREADY REFACTORED THIS->NBIN OUT
   return this->partition_windows();
#  if 0
   if( !fout.is_open() ) fout.open("MC_WangLandau.txt");
   int nel = energy.size();
   if( nel<2 ) 
   {
      fout << __FILE__ << ":" << __LINE__ << " partition_windows nel=" << nel << std::endl;
      return partition_windows();
   }
   if( est_lng.size()!=nel ) return;
   Elo = energy.front();
   Ehi = energy.back();
   Ewin.resize(2*NWindow);
   int imax = 0;
   double lng_max = est_lng[imax];
   for(int i=0; i<nel; i++)
      if( est_lng[i]>lng_max ) { imax=i; lng_max=est_lng[i]; }
   if( imax>(nel-5) )
   {
      // Increasing function
      double width = est_lng.back() - est_lng.front();
      width /= static_cast<double>(NWindow);
      int iwin = 0;
      int ipt = 0;
      while( iwin<NWindow )
      {
         Ewin[2*iwin] = energy[ipt];
         double next = est_lng[ipt] + width;
         while( ipt<nel && est_lng[ipt]<next ) ipt++;
         iwin++;
      }
      for(iwin=0; iwin<(NWindow-1); iwin++)
      {
         Ewin[2*iwin+1] = Ewin[2*(iwin+1)] + fwinover*(Ewin[2*(iwin+2)]-Ewin[2*(iwin+1)]);
         if( Ewin[2*iwin+1] > Ehi ) Ewin[2*iwin+1] = Ehi;   // last window can't go outside range
         if( Ewin[2*iwin] >= Ewin[2*iwin+1] ) Ewin[2*iwin] = Ewin[2*iwin-1] - fwinover*(Ewin[2*iwin-1]-Ewin[2*iwin-2]);
      }   
      Ewin[2*(NWindow-1)+0] = (Ewin[2*(NWindow-2)+1]+Ewin[2*(NWindow-2)])/2.; // Make last window stretch
      Ewin[2*(NWindow-1)+1] = Ehi;                                            // From middle of neighbor to end
      Ewin[0] = Elo;                                                          // When Nwindow==1
   }
   else
   {
      // Maximum in-between
      fout << __FILE__ << ":" << __LINE__ << " not implemented" << std::endl;
      partition_windows();
   }
   return;
#  endif
}


template<typename Walker>
void MC_WangLandau::init_pool(std::vector<Walker>& walkerpool)
{
   bool verbose = false;
   // Set up the windows (here NWindow is global number of windows)
   if( this->NWalkPerProcess<1 ) this->NWalkPerProcess=1;
   int NWalker = 1;
   if(NWindow<1) NWindow = 1;
   mp_window.init(NWindow);
   re_iter = 0;
   // Build walker information
   if( mp_window.pool.nproc==1 )
   {
      // serial
      if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " Serial Window allocation" << std::endl;
      beginwin.resize(NWindow+1);
      int NWalkPerWindow = NWalkPerProcess/NWindow;
      if( NWalkPerWindow<1 )
      {
         std::cout << __FILE__ << ":" << __LINE__ << " Fewer walkers than windows" << std::endl;
         NWalkPerWindow = 1;
         NWalkPerProcess = NWindow;
      }
      for(int i=0; i<NWindow; i++) beginwin[i] = i*NWalkPerWindow;
      beginwin[NWindow] = NWalkPerProcess;
      NWalker = NWalkPerProcess;
   }
   else
   {
      // parallel, all walkers in one window
      beginwin.resize(2);
      beginwin[0] = 0;
      beginwin[1] = this->NWalkPerProcess;
      NWalker = this->NWalkPerProcess;
   }
   // Resize the walker pool
   walkerpool.resize(NWalker,walkerpool[0]);
   // Default windows is uniform length in energy
   if( Ewin.size()<2 )
   {
      if( (Elo==+std::numeric_limits<double>::max()) || (Ehi==-std::numeric_limits<double>::max()) )
      {
         if(mp_window.pool.iproc==0) std::cout << __FILE__ << ":" << __LINE__ << " Must set Elo and Ehi" << std::endl;
         return;
      }
      if(Elo==Ehi)
      {
         if(mp_window.pool.iproc==0) std::cout << __FILE__ << ":" << __LINE__ << " Elo==Ehi" << std::endl;
         return;
      }
      if(verbose && mp_window.pool.iproc==0) std::cout << __FILE__ << ":" << __LINE__ << " Global Window = [" << Elo << "," << Ehi << "]" << std::endl;
      partition_windows();
   }
   if(false) std::cout << __FILE__ << ":" << __LINE__ << "(" << mp_window.pool.iproc << ")" << std::endl;
   ittm.mp = mp_window.pool;
   ittm.init(Elo,Ehi,Ebin);
   if(false) std::cout << __FILE__ << ":" << __LINE__ << "(" << mp_window.pool.iproc << ")" << std::endl;
   if( mp_window.pool.nproc==1 )
   {
      // Serial
      const int NWindow_Local = beginwin.size()-1;
      for(int iwin=0; iwin<NWindow_Local; iwin++)
      {
         for(int iwalk=beginwin[iwin]; iwalk<beginwin[iwin+1]; iwalk++) 
         {
            walkerpool[iwalk].window.set_delta(Ewin[2*iwin],Ewin[2*iwin+1],iwin,Ebin);
            walkerpool[iwalk].WEmin = walkerpool[iwalk].window.Elo; // Ewin[2*iwin];
            walkerpool[iwalk].WEmax = walkerpool[iwalk].window.Ehi; // Ewin[2*iwin+1];
            walkerpool[iwalk].clear_histogram();
            walkerpool[iwalk].clear_entropy();
            walkerpool[iwalk].iwalk_global = iwalk;
         }
         if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iwin=" << iwin << " iwalk=" << beginwin[iwin] << " " 
                               << " binwindow= " << walkerpool[beginwin[iwin]].window.Elo << "," << walkerpool[beginwin[iwin]].window.Ehi << " " 
                               << " bracket= " << walkerpool[beginwin[iwin]].WEmin << "," <<walkerpool[beginwin[iwin]].WEmax << std::endl;
         if(walkerpool[beginwin[iwin]].S.size()<10) 
            std::cout << "iwalk=" << beginwin[iwin] << " S.size=" << walkerpool[beginwin[iwin]].S.size() << std::endl;
         if(walkerpool[beginwin[iwin]].h.size()<10) 
            std::cout << "iwalk=" << beginwin[iwin] << " h.size=" << walkerpool[beginwin[iwin]].h.size() << std::endl;
      }
   }
   else
   {
      // Parallel
      if(verbose) std::cout << __FILE__ << ":" << __LINE__ <<  " Parallel windows NWindow=" << NWindow << " Global Window = [" << Elo << "," << Ehi << "]" << std::endl;
      int iwin = mp_window.igang;
      for(int iwalk=0; iwalk<NWalker; iwalk++)
      {
         // walkerpool[iwalk].window.set(Ewin[2*iwin],Ewin[2*iwin+1],iwin,nbin);
         walkerpool[iwalk].window.set_delta(Ewin[2*iwin],Ewin[2*iwin+1],iwin,Ebin);
         walkerpool[iwalk].WEmin = walkerpool[iwalk].window.Elo; // Ewin[2*iwin];
         walkerpool[iwalk].WEmax = walkerpool[iwalk].window.Ehi; // Ewin[2*iwin+1];
         walkerpool[iwalk].clear_histogram();
         walkerpool[iwalk].clear_entropy();
         walkerpool[iwalk].iwalk_global = NWalker*(mp_window.pool.iproc)+iwalk;
      }
      if(mp_window.gang.iproc==0)
      {
         std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << " iwalk=" << walkerpool[0].iwalk_global << " " << walkerpool[0].window.Elo << "," << walkerpool[0].window.Ehi << " bin=" << walkerpool[0].window.Ebin << "/" << walkerpool[0].window.NBin << " iwin=" << walkerpool[0].window.iwindow << " WE=" << walkerpool[0].WEmin << "," << walkerpool[0].WEmax << std::endl;
      }
   } 
   for(int iwalk=0; iwalk<NWalker; iwalk++)
   {
      walkerpool[iwalk].clear_histogram();
      walkerpool[iwalk].clear_entropy();
   }
}


template<typename Model, typename Walker>
void MC_WangLandau::DoReplicaExchange_serial(Model& model, std::vector<Walker>& walkerpool)
{
   bool verbose = false;
   if( mp_window.pool.nproc!=1 )
   {
      std::cout << __FILE__ << ":" << __LINE__ << " Serial implentation called in parallel environment" << std::endl;
      return;
   }
   // alternate between communication to left and right 
   bool send_right = re_iter % 2;                     // 0,1,0,1,0,1,0,1,0,1,0,1,...
   // Serial
   const int NWindow_Local = beginwin.size()-1;
   for(int iwin=0; iwin<NWindow_Local; iwin++)
   {
      bool when_right = (iwin+1) % 2;
      if( send_right )
      { 
         if( when_right && iwin>1 )
            DoReplicaExchange_serial(model,walkerpool,beginwin[iwin-1],beginwin[iwin],beginwin[iwin+1]);
      }
      else
      {
         if( !when_right && iwin<(NWindow_Local-1) ) 
            DoReplicaExchange_serial(model,walkerpool,beginwin[iwin-1],beginwin[iwin],beginwin[iwin+1]);
      }
   }
}

template<typename Model, typename Walker>
void MC_WangLandau::DoReplicaExchange_serial(Model& model, std::vector<Walker>& walkerpool, int iwalkL, int iwalkM, int iwalkR)
{
   int nwalk = std::min(iwalkM-iwalkL,iwalkR-iwalkM);
   if( nwalk<1 ) 
   {
      std::cout << __FILE__ << ":" << __LINE__ << " error in expected order of " << iwalkL << " < " << iwalkM << " < " << iwalkR << std::endl;
      return;
   }
   int comm_shift = (re_iter/2) % nwalk;
   for(int iwalk=0; iwalk<nwalk; iwalk++)
   {
      // Thomas Vogel, Ying Wai Li, Thomas Wust, and David P. Landau
      // Generic, Hierarchical Framework for Massively Parallel Wang-Landau Sampling
      // Phys. Rev. Lett. 110, 210603 (2013)
      // "Movtivated by the detailed balance condition the acceptance probability
      // for the exchange of configurations X and Y between walkers i and j is
      // Pacc = min[1, g_i(E(X))/g_i(E(Y)) * g_j(E(Y))/g_j(E(X)) ]"
      // Note: X and Y are the configurations of walkers i and j, respectively
      // i=X=right (this branch)   j=Y=left (other branch)
      // Pacc = min[1, g_R(E_R)/g_R(E_L) * g_L(E_L)/g_L(E_R) ]
      // Remember add lndos to multiply dos (every adder can multiply on a log table)
      int ileft = iwalkL + iwalk;
      int irgt = iwalkM + ( (iwalk+comm_shift) % nwalk );
      double E_left = walkerpool[ileft].now.E;
      double E_rght = walkerpool[irgt].now.E;
      double Pacc = walkerpool[irgt].get_lndos(E_rght,LinearInterp) - walkerpool[irgt].get_lndos(E_left,LinearInterp)
                  + walkerpool[ileft].get_lndos(E_left,LinearInterp) - walkerpool[ileft].get_lndos(E_rght,LinearInterp);
      if( Pacc>0 ) Pacc = 0;                                                                   // Pacc = min[ ln(1), ln(ratio_R*ratio_L) ]
      bool accept = static_cast<int>( urng()<std::exp(Pacc) );                                 // convert probability into an action
      if( accept )
      {
         // swap microstate configurations
         Walker tmp(walkerpool[ileft]);
         walkerpool[ileft].sigma = walkerpool[irgt].sigma;
         model.calc_observable(walkerpool[ileft].sigma,walkerpool[ileft].now);
         walkerpool[irgt].sigma = tmp.sigma;
         model.calc_observable(walkerpool[irgt].sigma,walkerpool[irgt].now);
         
      }
   }
}


 
template<typename Model, typename Walker>
void MC_WangLandau::DoReplicaExchange(Model& model, std::vector<Walker>& walkerpool)
{
   bool verbose = false;
#  ifndef USE_MPI
   DoReplicaExchange_serial(model,walkerpool);
   return;
#  else
   if(!mp_window.pool.in()) return;
   int nwalk = walkerpool.size();
   if(nwalk<1) return;
   int nspin = walkerpool[0].sigma.size();
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << std::endl;
   // alternate between communication to left and right 
   bool send_right = re_iter % 2;                     // 0,1,0,1,0,1,0,1,0,1,0,1,...
   bool when_right = (mp_window.igang+1) % 2;
   // offset in communication pattern
   int npergang = mp_window.gang.nproc;
   int comm_shift = (re_iter/2) % npergang;
   // The actual communication
   MPI_Status status;
   int jproc_gang = mp_window.gang.nproc;
   int jproc_pool = mp_window.pool.nproc;
   MPI_Datatype MPIConfigType = MPITypeTraits<typename Walker::ConfigType::value_type>::mpitype;
   int twonwalk = 2*nwalk;
   double dbuffer[4*nwalk];                           // (E_left,E_right,lng_right,lng_left/lng_right)
   double *E_left    = dbuffer + 0*nwalk;
   double *E_right   = dbuffer + 1*nwalk;
   double *lng_right = dbuffer + 2*nwalk;
   double *lng_ratio = dbuffer + 3*nwalk;
   int    ibuffer[nwalk];
   int    *accept    = ibuffer;
   if( send_right==when_right )                       // true,true or false,false
   { 
      // Branch for sender that has lower iproc_pool, left
      if( mp_window.igang<(NWindow-1) )
      {
         // look for jproc to right on the number line        
         jproc_gang = (mp_window.gang.iproc + comm_shift) % npergang;
         jproc_pool = npergang*(mp_window.igang+1) + jproc_gang;
      }
      if( 0<=jproc_pool && jproc_pool<mp_window.pool.nproc )
      {  
         if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " RE Left iproc=" << mp_window.pool.iproc << " jproc=" << jproc_pool << std::endl;
         // STEP 1: Swap E_left and (E_right,lng_right)
         for(int iwalk=0; iwalk<nwalk; iwalk++)
            E_left[iwalk] = walkerpool[iwalk].now.E;
         MPI_Send(E_left,nwalk,MPI_DOUBLE,jproc_pool,0,mp_window.pool.comm);
         MPI_Recv(E_right,twonwalk,MPI_DOUBLE,jproc_pool,1,mp_window.pool.comm,&status);
         // STEP 2: Return lng_left/lng_right
         MPI_Recv(lng_ratio,nwalk,MPI_DOUBLE,jproc_pool,2,mp_window.pool.comm,&status);
         // STEP 3: Send Accept
         // Thomas Vogel, Ying Wai Li, Thomas Wust, and David P. Landau
         // Generic, Hierarchical Framework for Massively Parallel Wang-Landau Sampling
         // Phys. Rev. Lett. 110, 210603 (2013)
         // "Movtivated by the detailed balance condition the acceptance probability
         // for the exchange of configurations X and Y between walkers i and j is
         // Pacc = min[1, g_i(E(X))/g_i(E(Y)) * g_j(E(Y))/g_j(E(X)) ]"
         // Note: X and Y are the configurations of walkers i and j, respectively
         // i=X=right (this branch)   j=Y=left (other branch)
         // Pacc = min[1, g_R(E_R)/g_R(E_L) * g_L(E_L)/g_L(E_R) ]
         // Remember add lndos to multiply dos (every adder can multiply on a log table)
         for(int iwalk=0;iwalk<nwalk; iwalk++)
         {
            // from below: lng_ratio[iwalk] = walkerpool[iwalk].get_lndos(E_left[iwalk]) - lng_right[iwalk];
            //             this is g_j(E(Y))/g_j(E(X)) = g_L(E_L)/g_L(E_R) 
            double Pacc = walkerpool[iwalk].get_lndos(E_right[iwalk],LinearInterp) - walkerpool[iwalk].get_lndos(E_left[iwalk],LinearInterp)
                        + lng_ratio[iwalk];
            if( Pacc>0 ) Pacc = 0;                                       // Pacc = min[ ln(1), ln(ratio_R*ratio_L) ]
            accept[iwalk] = static_cast<int>( urng()<std::exp(Pacc) );   // convert probability into an action
         }
         MPI_Send(accept,nwalk,MPI_INT,jproc_pool,3,mp_window.pool.comm);
         // STEP 4: If Accept then swap configuration
         for(int iwalk=0; iwalk<nwalk; iwalk++)
         {
            if( accept[iwalk] )
            {
               if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << " swapping " << iwalk << " with jproc=" << jproc_pool << std::endl;
               MPI_Send(&(walkerpool[iwalk].sigma[0]),nspin,MPIConfigType,jproc_pool,iwalk+10000,mp_window.pool.comm);
               MPI_Recv(&(walkerpool[iwalk].sigma[0]),nspin,MPIConfigType,jproc_pool,iwalk+20000,mp_window.pool.comm,&status);
               model.calc_observable(walkerpool[iwalk].sigma,walkerpool[iwalk].now);
               if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << " processed " << iwalk << std::endl;
            }
         }
         if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << std::endl;
      }
   }
   else
   {
      // Branch for reciever that has higher iproc_pool, right
      if( mp_window.igang>0 )
      {
         // look for jproc to left on the number line
         jproc_gang = (mp_window.gang.iproc + (npergang-comm_shift)) % npergang;
         jproc_pool = npergang*(mp_window.igang-1) + jproc_gang;
      }
      if( 0<=jproc_pool && jproc_pool<mp_window.pool.nproc )
      {  
         if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " RE Right iproc=" << mp_window.pool.iproc << " jproc=" << jproc_pool << std::endl;
         // STEP 1: Swap E_left and (E_right,lng_right)
         for(int iwalk=0; iwalk<nwalk; iwalk++)
         {
            E_right[iwalk]   = walkerpool[iwalk].now.E;
            lng_right[iwalk] = walkerpool[iwalk].get_lndos(E_right[iwalk],LinearInterp);
         }
         MPI_Recv(E_left,nwalk,MPI_DOUBLE,jproc_pool,0,mp_window.pool.comm,&status);
         MPI_Send(E_right,twonwalk,MPI_DOUBLE,jproc_pool,1,mp_window.pool.comm);
         // STEP 2: Return lng_left/lng_right
         for(int iwalk=0; iwalk<nwalk; iwalk++)
         {
            lng_ratio[iwalk] = walkerpool[iwalk].get_lndos(E_left[iwalk],LinearInterp) - lng_right[iwalk];
         }
         MPI_Send(lng_ratio,nwalk,MPI_DOUBLE,jproc_pool,2,mp_window.pool.comm);
         // STEP 3: Send Accept
         MPI_Recv(accept,nwalk,MPI_INT,jproc_pool,3,mp_window.pool.comm,&status);
         // STEP 4: If Accept then swap configuration
         for(int iwalk=0; iwalk<nwalk; iwalk++)
         {
            if( accept[iwalk] )
            {
               if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << " swapping " << iwalk << " with jproc=" << jproc_pool << std::endl;
               typename Walker::ConfigType config_temp(walkerpool[iwalk].sigma.size());
               //typename Walker::ConfigType config_temp = walkerpool[iwalk].sigma;
               MPI_Recv(&(config_temp[0]),nspin,MPIConfigType,jproc_pool,iwalk+10000,mp_window.pool.comm,&status);
               MPI_Send(&(walkerpool[iwalk].sigma[0]),nspin,MPIConfigType,jproc_pool,iwalk+20000,mp_window.pool.comm);
               walkerpool[iwalk].sigma = config_temp;
               model.calc_observable(walkerpool[iwalk].sigma,walkerpool[iwalk].now);
               if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << " processed " << iwalk << std::endl;
            }
         }
         if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << mp_window.pool.iproc << std::endl;
      }
   }
   // increment the iteration counter
   re_iter++;
#  endif
}

#endif
