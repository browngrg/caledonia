#ifndef MFIA_HAMILTONIAN_HPP
#define	MFIA_HAMILTONIAN_HPP

// Implementation of the "MFIA" Model -- Ising Antiferromagnet with a mean-field interaction
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// Sept 17, 2014:      Testing JSIGN vs -JSIGN
//
// History:
//
// Created on August 21, 2013, 2:38 PM
// February 28, 2014:  Changed to 2d binning
// March 31, 2014:     Changed to binning on N1A and N1B, sublattice occupations
// June 28, 2014:      Adpated to calendonia MC interface
// July 14, 2014:      Implementation of the "MFIA" Model -- Ising Antiferromagnet with a mean-field interaction
// Sept 29, 2014:      Added MfiaMF Hamiltonian for MF exchange interaction
//

#include<vector>
#include<map>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>


struct Mfia_Observables {
public:
    double E,X;   // Required: Energy of the walker and order parameter
    int M;        // Magnetization = sum over spins
    int N;        // Sum over nearest-neighbor spins
    int P;        // Number of spins aligned with staggered pattern
    int R;        // Auxillary variable (row)
    int C;        // Number of up spins (column)
    int U;        // Number of ungerade pairs
    int V;        // Number of spins
    double E_N;
    double E_H;
    double E_A;
public:
    Mfia_Observables() { E=0; X=0; } 
    Mfia_Observables(const Mfia_Observables& w) { this->copy(w); }
    void operator=(const Mfia_Observables& w) { this->copy(w); }
    void copy(const Mfia_Observables& w)
    {
        E = w.E;
        X = w.X;
        M = w.M;
        N = w.N;
        P = w.P;
        R = w.R;
        C = w.C;
        U = w.U;
        V = w.V;
        E_N = w.E_N;
        E_H = w.E_H;
        E_A = w.E_A;
    }
    bool operator==(const Mfia_Observables& w) const
    {
        return ( std::fabs(E-w.E)<1.e-8 && X==w.X && M==w.M && N==w.N && P==w.P && R==w.R && U==w.U );
    }
    bool operator!=(const Mfia_Observables& w) const { return !((*this)==w); }
};


// Stirling's formula for the factorial (approximation for large x)
double lnfac(double x) 
{ 
   if(x<10) 
   {
      static double fac[11] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };
      return log(fac[static_cast<int>(x)]); 
   }
   else 
   {
       return x*log(x)-x; 
   }
}

// Return the log of the sum of two numbers given as logs
double add_logs(double log1, double log2)
{
   using namespace std;
   if(log1<log2) swap(log1,log2);
   double negarg = log2-log1;
   double factor = 1 + exp(negarg);
   return log1 + log(factor);
}



class Mfia_Hamiltonian                  // and if is AFM or FM model
{
public:

    typedef Mfia_Observables     Observables;
    typedef std::vector<short>   Config;

public:

    // Basic Ising-model properties
    int D;
    int L;                       // Side of the lattice
    int V;                       // Volume = Total number of spins
    std::vector<int> jspin;      // List of neighbors
    int JSIGN;                   // FM=+1, AFM=-1
    std::vector<int> stagmask;   // Weights for staggered magnetization

    // The values of the Fields
    double H;                    // The external field
    double A;                    // The mean-field
    double AhbV;                 // A/2/V
    
public:
    
   Mfia_Hamiltonian(int DIM=2, int Lval=8, bool AFM=true);

   ~Mfia_Hamiltonian()
   {
   }

   std::string header() const;
    
   void init(bool verbose=false);
    
public:
    
   // Some code will want this, Mfia ignores it
   double step_size;

   // Calculate all macroscopic variables based on this sigma
   void calc_observable(Config& sigma, Mfia_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   void change(Config& sigma, Mfia_Observables& macro, int ispin, short new_sigma);
    
   // Make an attempted Monte Carlo step
   enum { SPINFLIP = true };
   template<typename MCWalker, typename URNG> void mc_step(MCWalker& walker, URNG& urng);
   template<typename MCWalker, typename URNG> void spin_flip(MCWalker& walker, URNG& urng);
   template<typename MCWalker, typename URNG> void spin_exch(MCWalker& walker, URNG& urng);


public:

   // Create mixed ferromagnetic / antiferromagnetic configuration
   void initial_mixed(Config& sigma) const;

   void initial_ferro(Config& sigma) const;

public:

   void estimate_lng(std::vector<double>& energy, std::vector<double>& lng, int nbin=200) const;

};


 
Mfia_Hamiltonian::Mfia_Hamiltonian(int DIM, int Lval, bool AFM)
{
   D = DIM;
   L = Lval;
   JSIGN = -1;
   if(!AFM) JSIGN = +1;
   A = H = 0;
   this->init();
}


std::string Mfia_Hamiltonian::header() const
{
   char buffer[1024];
   return std::string(buffer);
}

 
void Mfia_Hamiltonian::init(bool verbose)
{
   // Extent of system
   if( D<1 || D>4 ) {
      std::cerr << "Mfia_Hamiltonian::set_size():" << std::endl;
      std::cerr << "Invalid dimension D=" << D<< std::endl;
      D = 2;
      std::cerr << "Using value D=" << D<< std::endl;
   }
   if( (L%2)!=0 )
   {
      std::cerr << "Mfia_Hamiltonian::set_size() assumes L even" << std::endl;
      std::cerr << "specified value is L=" << L << std::endl;
      L--;
      std::cerr << "using L=" << L << std::endl;
   }
   int lda[D+1];
   lda[0] = 1;
   for(int idim=0; idim<D; idim++) lda[idim+1] = L*lda[idim];
   V = lda[D];
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " NSpin=" << V << std::endl;
   AhbV = A/2./static_cast<double>(V);   // A/2/V
   // Neighbor table 
   jspin.resize(2*D*V);
   stagmask.resize(V);
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " resizes successful" << std::endl;
   for(int ispin=0; ispin<V; ispin++) 
   {
      // coordinates on lattice
      int icoord[D];
      int isite = ispin;
      for(int idim=D-1; idim>=0; idim--) 
      {
         icoord[idim] = isite/lda[idim];
         isite -= icoord[idim]*lda[idim];
      }
      // forward neighbors
      for(int idim=0; idim<D; idim++) 
      {
         int jcoord[D];
         for(int k=0; k<D; k++) jcoord[k]=icoord[k];
         jcoord[idim]++;
         if( jcoord[idim]>=L ) jcoord[idim]=0;
         int jsite = 0;
         for(int k=0; k<D; k++) jsite += jcoord[k]*lda[k];
         jspin[2*D*ispin+idim] = jsite;
      }
      // backward neighbors
      for(int idim=0; idim<D; idim++) 
      {
         int jcoord[D];
         for(int k=0; k<D; k++) jcoord[k]=icoord[k]; jcoord[idim]--; if( jcoord[idim]<0 ) jcoord[idim]=L-1; int jsite = 0; for(int k=0; k<D; k++) jsite += jcoord[k]*lda[k];
         jspin[2*D*ispin+D+idim] = jsite;
      }
      // staggered order parameter weights
      int csum = 0;
      for(int idim=0; idim<D; idim++) csum += icoord[idim];
      int odd = csum%2;
      if( odd==1 )
         stagmask[ispin] = -1;
      else
         stagmask[ispin] = +1;
   }
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " init done" << std::endl;
}


void Mfia_Hamiltonian::calc_observable(Config& sigma, Mfia_Observables& macro)
{
   short* S(&(sigma[0]));
   int M = 0;
   int N = 0;
   int P = 0;
   for(int ispin=0; ispin<V; ispin++) 
   {
       M += S[ispin];
       P += stagmask[ispin]*S[ispin];
       for(int k=0; k<D; k++)
           N += S[ispin]*S[jspin[2*D*ispin+k]];
   }
   macro.M = M;
   macro.N = N;
   macro.P = P;
   macro.C = (V+M)/2;
   macro.U = (D*V-N)/2;
   macro.R = macro.U/2;
   macro.E = -JSIGN*N-H*M-AhbV*M*M; // JSIGN*N-H*M-AhbV*M*M; 
   macro.E_N = -JSIGN*N;            // JSIGN*N;
   macro.E_H = -H*M;
   macro.E_A = -AhbV*M*M;
   macro.X = P;
   macro.V = V;
}


void Mfia_Hamiltonian::change(Mfia_Hamiltonian::Config& sigma, Mfia_Observables& macro, int ispin, short new_sigma)
{
   short* S(&(sigma[0]));
   // calculate changes
   int Ni = 0;
   for(int k=0; k<2*D; k++)
       Ni += S[ispin]*S[jspin[2*D*ispin+k]];
   int dN = -2*Ni;
   int dR = -S[ispin];
   int dM = 2*dR;
   int dP = 2*stagmask[ispin]*dR;
   int dM2 = dM*dM + 2*macro.M*dM;        // (M+dM)^2-M^2
   // Assign new values
   macro.E += -JSIGN*dN-H*dM-AhbV*dM2;    // JSIGN*dN-H*dM-AhbV*dM2;
   macro.M += dM;
   macro.N += dN;
   macro.P += dP;
   macro.C -= S[ispin];
   macro.U -= dN/2;
   macro.R = macro.U/2;
   macro.X = macro.P;
   macro.E_N = -JSIGN*macro.N;            // JSIGN*macro.N;
   macro.E_H = -H*macro.M;
   macro.E_A = -AhbV*macro.M*macro.M;
   macro.V   = V;
   // Flip the spin
   S[ispin] = new_sigma;
}
 

template<typename MCWalker, typename URNG> 
void Mfia_Hamiltonian::mc_step(MCWalker& walker, URNG& urng)
{
   if( SPINFLIP )
      this->spin_flip(walker,urng);
   else
      this->spin_exch(walker,urng);
}


template<typename MCWalker, typename URNG> 
void Mfia_Hamiltonian::spin_exch(MCWalker& walker, URNG& urng)
{
   // Randomly choose a move
   int NSite = walker.sigma.size();
   int isite = static_cast<int>(NSite*urng());
   while( isite<0 || isite>=NSite )
   {
      //std::cerr << __FILE__ << ":" << __LINE__ << " markov_propose bad isite=" << isite << "/" << NSite << std::endl;
      isite = static_cast<int>(NSite*urng());
   }
   // Random choose a neighbor
   int nnum = 2*D;
   int jsite = static_cast<int>(nnum*urng());
   while( jsite<0 || jsite>=nnum )
   {
      //std::cerr << __FILE__ << ":" << __LINE__ << " markov_propose bad jsite=" << jsite << "/" << nnum << std::endl;
      jsite = static_cast<int>(nnum*urng());
   }
   jsite = jspin[nnum*isite+jsite];  // convert local neighbor number to site index
   // Do the move
   Mfia_Observables before_exchg(walker.now);   // save macrostate before composing invidual moves
   walker.old.copy(before_exchg);               // This is the old macrostate for the Monte Carlo step
   walker.save_initial();
   walker.add_change(isite);
   walker.add_change(jsite);
   char s_i = walker.sigma[isite];
   char s_j = walker.sigma[jsite];
   this->change(walker.sigma,walker.now,isite,s_j);      // spin_flip updates walker.old as build new state
   this->change(walker.sigma,walker.now,jsite,s_i);
}

template<typename MCWalker, typename URNG> 
void Mfia_Hamiltonian::spin_flip(MCWalker& walker, URNG& urng)
{
   // Randomly choose a move
   int NSite = walker.sigma.size();
   int isite = static_cast<int>(NSite*urng());
   while( isite<0 || isite>=NSite )
   {
      //std::cerr << __FILE__ << ":" << __LINE__ << " mc_step bad isite=" << isite << "/" << NSite << std::endl;
      isite = static_cast<int>(NSite*urng());
   }
   char sold =  walker.sigma[isite];
   char snow = -walker.sigma[isite];
   walker.save_initial();
   walker.add_change(isite);
   this->change(walker.sigma,walker.now,isite,snow);
   if( false )
   {
      // Check that the method gets the right total energy
      double Ealgorithm = walker.now.E;
      this->calc_observable(walker.sigma,walker.now);
      static int ialert = 0;
      if( std::fabs(Ealgorithm-walker.now.E)>1.e-5 && ialert<10 )
      {
         ialert++;
         std::cout << "Mfia_Hamiltonian::spin_flip gets wrong energy" << std::endl;
         std::cout << "Algorithm energy = " << Ealgorithm << std::endl;
         std::cout << "full_calc energy = " << walker.now.E << std::endl;
      } 
   }
}


void Mfia_Hamiltonian::initial_mixed(Mfia_Hamiltonian::Config& sigma) const
{
   this->initial_ferro(sigma);
#  if 0
   // Slab
   for(int ispin=0; ispin<V; ispin++) 
   {
      if( (ispin%L)>(L/2) )
         sigma[ispin] *= stagmask[ispin];
   }
#  else
   // Circle
   int LX = L;
   int LY = (D>1)? L : 1;
   int LZ = (D>2)? L : 1;
   long rad2 = (L/2)*(L/2)-1;
   if( D>3 ) { std::cout << __FILE__ << ":" << __LINE__ << " assumes D<=3" << std::endl; }
   for(int iz=0; iz<LZ; iz++)
   {
      long iz2 = (iz-LZ/2)*(iz-LZ/2);
      for(int iy=0; iy<LY; iy++)
      {
         long iy2 = (iy-LY/2)*(iy-LY/2);
         for(int ix=0; ix<LX; ix++)
         {
            long ix2 = (ix-LX/2)*(ix-LX/2);
            long ri2 = ix2+iy2+iz2;
            if( ri2<rad2 )
            {
               int ispin = ix + iy*LX + iz*LX*LY;
               sigma[ispin] *= stagmask[ispin];
            }
         }
      }
   }
#  endif
}


void Mfia_Hamiltonian::initial_ferro(Mfia_Hamiltonian::Config& sigma) const
{
   sigma.resize(V);
   for(int ispin=0; ispin<V; ispin++) sigma[ispin] = +1;
}


void Mfia_Hamiltonian::estimate_lng(std::vector<double>& energy, std::vector<double>& lng, int nbin) const
{
   std::map<int,double> lngmap;
   int LL = 64;
   int LN = LL*LL;
   if( D==1 ) { LN=LL=128; }
   if( D==3 ) { LL=16; LN=LL*LL*LL; }
   if( D==4 ) { LL=8; LN=LL*LL*LL*LL; }
   if( D>4 )  { std::cerr << __FILE__ << ":" << __LINE__ << " DIM<=4" << std::endl; return; }
   long halfN = LN/2;
   double ahalfN = A/static_cast<double>(2*LN);
   double fourN = 4./static_cast<double>(LN);
   std::vector<double> lngk(halfN+1,0);
   for(int k=0; k<=halfN; k++)
   {
      // lngk[k] = lng_Ising(halfN,k);
      // Returns the entropy of a lattice of N spins with k spins up
      // g(k) =  N! / [ k! * (N-k)! ]
      double lnfacN  = lnfac(halfN);
      double lnfack  = lnfac(k);
      double lnfacNk = lnfac(halfN-k);
      lngk[k] = lnfacN - lnfack - lnfacNk;
   }
   // Tabulate lng(E) for this system
   for(long thetaA=0; thetaA<=halfN; thetaA++)
   {
      long mA = 2*thetaA - halfN;
      for(long thetaB=0; thetaB<=halfN; thetaB++)
      {
         double lng_AB = lngk[thetaA]+lngk[thetaB];
         long mB = 2*thetaB - halfN;
         long M = -LN + 2*(thetaA+thetaB);
         double E = fourN*D*mA*mB - ahalfN*M*M-H*M;
         int iE = static_cast<int>(E);
         std::map<int,double>::iterator itr = lngmap.find(iE);
         if( itr==lngmap.end() )
            lngmap[iE] = lng_AB;
         else
            lngmap[iE] = add_logs(lngmap[iE],lng_AB);
      }
   }
   // Construct the estimate DOS
   if(nbin<50) nbin=50;
   if(nbin>10000) nbin=10000;
   double Elo = lngmap.begin()->first;
   double Ehi = lngmap.rbegin()->first;
   double Ebin = (Ehi-Elo)/static_cast<double>(nbin);
   energy.resize(nbin);
   lng.resize(nbin);
   for(int i=0; i<nbin; i++) 
   {
      energy[i] = Elo + i*Ebin;
      lng[i]=0;
   }
   for(std::map<int,double>::const_iterator itr=lngmap.begin(); itr!=lngmap.end(); itr++)
   {
      int ibin = ((itr->first)-Elo)/Ebin;
      if( ibin<0 ) ibin=0;
      if( ibin>=nbin) ibin=nbin;
      double lng_AB = itr->second;
      if( lng[ibin]<=0 )
         lng[ibin] = lng_AB;
      else
         lng[ibin] = add_logs(lng[ibin],lng_AB);
   }
   // Simple volume scaling
   double scale = static_cast<double>(this->V)/static_cast<double>(LN);
   for(int i=0; i<nbin; i++)
   {
      energy[i] *= scale;
      lng[i] *= scale;
   }
}


////////////////////////////////////////////////////////////////////////////////
// MfiaMF2: Mean-field exchange
//          Configuration only tracks the sublattice coverage
//          Note that MfiaMF model uses all spins, but only MF exchange
////////////////////////////////////////////////////////////////////////////////

class MfiaMF2_Hamiltonian                  // and if is AFM or FM model
{
public:

    typedef Mfia_Observables     Observables;
    typedef std::vector<int>     Config;

public:

    // Basic Ising-model properties
    int D;
    int L;                       // Side of the lattice
    int V;                       // Volume = Total number of spins

    // The values of the Fields
    double H;                    // The external field
    double A;                    // The mean-field
    int    JSIGN;               
    double AhbV;                 // A/2/V
    int    halfV;                // V/2    
    double fourV;                // 4/V;

public:
    
   MfiaMF2_Hamiltonian(int DIM=2, int Lval=8, bool AFM=true);

   std::string header() const;
    
   void init(bool verbose=false);
    
public:
    
   // Some code will want this, Mfia ignores it
   double step_size;

   // Calculate all macroscopic variables based on this sigma
   void calc_observable(Config& sigma, Mfia_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   void change(Config& sigma, Mfia_Observables& macro, int ispin, int new_sigma);
    
   // Make an attempted Monte Carlo step
   template<typename MCWalker, typename URNG> void mc_step(MCWalker& walker, URNG& urng);

public:

   // MC move unique to MfiaMF (since exchange is only between average spins)
   template<typename MCWalker> void swap_sublattices(MCWalker& walker);

   float swapf;   // fractional chance of swapping sublatties each attempted MC

public:

   // Create mixed ferromagnetic / antiferromagnetic configuration
   void initial_mixed(Config& sigma) const;

   void initial_ferro(Config& sigma) const;

};

std::string MfiaMF2_Hamiltonian::header() const
{
   std::string result;
   char buffer[1024];
   sprintf(buffer,"# MfiaMF2 Model J=%d H=%lf A=%lf\n",JSIGN,H,A);
   result += std::string(buffer);
   sprintf(buffer,"# swapf=%f\n",swapf);
   result += std::string(buffer);
   sprintf(buffer,"# Dim=%d L=%d N=%d\n",D,L,V);
   result += std::string(buffer);
   return result;
}

MfiaMF2_Hamiltonian::MfiaMF2_Hamiltonian(int DIM, int Lval, bool AFM)
{
   D = DIM;
   L = Lval;
   JSIGN = -1;
   if(!AFM) JSIGN = +1;
   A = H = 0;
   swapf = 0;
   this->init();
}

void MfiaMF2_Hamiltonian::init(bool verbose)
{
   // Extent of system
   if( D<1 || D>4 ) {
      std::cerr << "MfiaMF2_Hamiltonian::set_size():" << std::endl;
      std::cerr << "Invalid dimension D=" << D<< std::endl;
      D = 2;
      std::cerr << "Using value D=" << D<< std::endl;
   }
   if( (L%2)!=0 )
   {
      std::cerr << "MfiaMF2_Hamiltonian::set_size() assumes L even" << std::endl;
      std::cerr << "specified value is L=" << L << std::endl;
      L--;
      std::cerr << "using L=" << L << std::endl;
   }
   int lda[D+1];
   lda[0] = 1;
   for(int idim=0; idim<D; idim++) lda[idim+1] = L*lda[idim];
   V = lda[D];
   if(verbose) std::cout << __FILE__ << ":" << __LINE__ << " NSpin=" << V << std::endl;
   AhbV = A/2./static_cast<double>(V);   // A/2/V
   halfV = V/2;                          // V/2    
   fourV = 4/static_cast<double>(V);     // 4/V;
   JSIGN = -1;
}

void MfiaMF2_Hamiltonian::calc_observable(Config& sigma, Mfia_Observables& macro)
{
      long thetaA = sigma[0];
      long thetaB = sigma[1];
      long mA = 2*thetaA - halfV;
      long mB = 2*thetaB - halfV;
      long N = fourV*D*mA*mB;
      long M = -V + 2*(thetaA+thetaB);
      long P = 2*(thetaA-thetaB);
      float E = fourV*D*mA*mB - AhbV*M*M-H*M;
      macro.M = M;
      macro.N = N;
      macro.P = P;
      macro.C = (V+M)/2;
      macro.U = (D*V-N)/2;
      macro.R = macro.U/2;
      macro.E = -JSIGN*N-H*M-AhbV*M*M;
      macro.E_N = -JSIGN*N;           
      macro.E_H = -H*M;
      macro.E_A = -AhbV*M*M;
      macro.X = P;
      macro.V = V;
      // debugging
      if( false && ( M<-V || M>V ) )
      {
         std::cout << __FILE__ << ":" << __LINE__ << " bad macro:" <<std::endl
                   << "sigma  = " << sigma[0] << " " << sigma[1] << std::endl
                   << "thetaA = " << thetaA << std::endl
                   << "thetaB = " << thetaB << std::endl
                   << "mA     = " << mA     << std::endl
                   << "mB     = " << mB     << std::endl
                   << "M      = " << M      << std::endl
                   << "Phi    = " << P      << std::endl
                   << "N      = " << E      << std::endl;
      }
}

void MfiaMF2_Hamiltonian::change(Config& sigma, Mfia_Observables& macro, int ispin, int new_sigma)
{
   sigma[ispin] = new_sigma;
   calc_observable(sigma,macro);
}

template<typename MCWalker>
void MfiaMF2_Hamiltonian::swap_sublattices(MCWalker& walker)
{
   long tmp;
   tmp = walker.sigma[0];
   walker.sigma[0] = walker.sigma[1];
   walker.sigma[1] = tmp;
   calc_observable(walker.sigma,walker.now);
}
  
// Make an attempted Monte Carlo step
template<typename MCWalker, typename URNG> 
void MfiaMF2_Hamiltonian::mc_step(MCWalker& walker, URNG& urng)
{
   // Swap sublattices to flip staggered order parameter
   if( swapf>0 && urng()<swapf ) this->swap_sublattices(walker);
   // Randomly choose a move
   int NSite = walker.sigma.size();
   int isite = (urng()<0.5);
   int sold = walker.sigma[isite];
   bool upspin = ((halfV*urng())<sold);
   int snow;
   if( upspin )
      snow = sold - 1;
   else
      snow = sold + 1;
   if( snow<0 ) snow = 0;
   if( snow>halfV ) snow = halfV;
   walker.save_initial();
   walker.add_change(isite);
   this->change(walker.sigma,walker.now,isite,snow);
   // debugging
   if( false && ( ((snow-sold)*(snow-sold)>1) || snow<0 || snow>halfV || snow!=walker.sigma[isite] ) )
   {
      std::cerr << __FILE__ << ":" << __LINE__ << " more than one spin-flip" << std::endl
                << "isite = " << isite << "/" << NSite << std::endl
                << "sold = " << sold << std::endl
                << "snow = " << snow << std::endl
                << "sigma = " << walker.sigma[0] << " " << walker.sigma[1] << std::endl;
   }
   // std::cout << isite << " " << increase << " " << sold << " " << snow << " " << walker.sigma[isite] << std::endl;
}

void MfiaMF2_Hamiltonian::initial_mixed(Config& sigma) const
{
   sigma.resize(2); 
   sigma[0] = halfV;
   sigma[1] = 0;
}

void MfiaMF2_Hamiltonian::initial_ferro(Config& sigma) const
{
   sigma.resize(2);
   sigma[0] = halfV;
   sigma[1] = halfV;
}




#endif	// MFIA_HAMILTONIAN_HPP

