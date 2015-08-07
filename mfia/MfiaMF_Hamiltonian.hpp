#ifndef MFIAMF_HAMILTONIAN_HPP
#define	MFIAMF_HAMILTONIAN_HPP

// Implementation of the completely mean-field "MFIA" Model -- Ising Antiferromagnet with a mean-field interaction
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// Sept 17, 2014:      Branched Mean-Field-Exchange version from original MFIA model
//
// History:
//
// Created on August 21, 2013, 2:38 PM
// February 28, 2014:  Changed to 2d binning
// March 31, 2014:     Changed to binning on N1A and N1B, sublattice occupations
// June 28, 2014:      Adpated to calendonia MC interface
// July 14, 2014:      Major revision for MFIA model


#include<vector>
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
        E_N = w.E_N;
        E_H = w.E_H;
        E_A = w.E_A;
    }
    bool operator==(const Mfia_Observables& w) const
    {
        return ( E==w.E && X==w.X && M==w.M && N==w.N && P==w.P && R==w.R && U==w.U );
    }
    bool operator!=(const Mfia_Observables& w) const { return !((*this)==w); }
};

class MfiaMF_Hamiltonian                  // and if is AFM or FM model
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
    
   MfiaMF_Hamiltonian(int DIM=2, int Lval=8, bool AFM=true);

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
   template<typename MCWalker, typename URNG> void mc_step(MCWalker& walker, URNG& urng);

public:

   // Create mixed ferromagnetic / antiferromagnetic configuration
   void initial_mixed(Config& sigma) const;

   void initial_ferro(Config& sigma) const;

};


    
   MfiaMF_Hamiltonian::MfiaMF_Hamiltonian(int DIM, int Lval, bool AFM)
   {
      D = DIM;
      L = Lval;
      JSIGN = -1;
      if(!AFM) JSIGN = +1;
      A = H = 0;
      this->init();
   }


   std::string MfiaMF_Hamiltonian::header() const
   {
      char buffer[1024];
      sprintf(buffer,"# MfiaMF Model J=%d H=%lf A=%lf\n",JSIGN,H,A);
      sprintf(buffer,"# Dim=%d L=%d N=%d\n",D,L,V);
      return std::string(buffer);
   }

    
   void MfiaMF_Hamiltonian::init(bool verbose)
   {
      // Extent of system
      if( D<1 || D>4 ) {
         std::cerr << "MfiaMF_Hamiltonian::set_size():" << std::endl;
         std::cerr << "Invalid dimension D=" << D<< std::endl;
         D = 2;
         std::cerr << "Using value D=" << D<< std::endl;
      }
      if( (L%2)!=0 )
      {
         std::cerr << "MfiaMF_Hamiltonian::set_size() assumes L even" << std::endl;
         std::cerr << "specified value is L=" << L << std::endl;
         L--;
         std::cerr << "using L=" << L << std::endl;
      }
      if(verbose) std::cout << "MfiaMF: Completely mean-field -- including exchange" << std::endl;
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


   void MfiaMF_Hamiltonian::calc_observable(Config& sigma, Mfia_Observables& macro)
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
      // MfiaMeanFIeldDOS.cpp: float E = fourN*DIM*mA*mB - ahalfN*M*M-H*M;
      // Mfia_Hamiltonian.hpp: macro.E = -JSIGN*N-H*M-AhbV*M*M; 
      macro.E   = -JSIGN*D*(M*M-P*P)/static_cast<float>(V)-H*M-AhbV*M*M;
      macro.E_N = -JSIGN*D*(M*M-P*P)/static_cast<float>(V);
      macro.E_H = -H*M;
      macro.E_A = -AhbV*M*M;
      macro.X = P;
  }


   void MfiaMF_Hamiltonian::change(MfiaMF_Hamiltonian::Config& sigma, Mfia_Observables& macro, int ispin, short new_sigma)
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
      // macro.E += JSIGN*dN-H*dM-AhbV*dM2;
      macro.M += dM;
      macro.N += dN;
      macro.P += dP;
      macro.C -= S[ispin];
      macro.U -= dN/2;
      macro.R = macro.U/2;
      // MfiaMeanFIeldDOS.cpp: float E = fourN*DIM*mA*mB - ahalfN*M*M-H*M;
      // Mfia_Hamiltonian.hpp: macro.E = JSIGN*N-H*M-AhbV*M*M; 
      // Above macro.E   = -JSIGN*D*(M*M-P*P)/static_cast<float>(V)-H*M-AhbV*M*M;
      macro.X = macro.P;
      macro.E_N = -JSIGN*D*(macro.M*macro.M-macro.P*macro.P)/static_cast<float>(V);
      macro.E_H = -H*macro.M;
      macro.E_A = -AhbV*macro.M*macro.M;
      macro.E   = macro.E_N + macro.E_H + macro.E_A;
      // Flip the spin
      S[ispin] = new_sigma;
   }
    

   template<typename MCWalker, typename URNG> 
   void MfiaMF_Hamiltonian::mc_step(MCWalker& walker, URNG& urng)
   {
      // Randomly choose a move
      int NSite = walker.sigma.size();
      int isite = static_cast<int>(NSite*urng());
      while( isite<0 || isite>=NSite )
      {
         std::cerr << __FILE__ << ":" << __LINE__ << " mc_step bad isite=" << isite << "/" << NSite << std::endl;
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
            std::cout << "MfiaMF_Hamiltonian::spin_flip gets wrong energy" << std::endl;
            std::cout << "Algorithm energy = " << Ealgorithm << std::endl;
            std::cout << "full_calc energy = " << walker.now.E << std::endl;
         } 
      }
   }


void MfiaMF_Hamiltonian::initial_mixed(MfiaMF_Hamiltonian::Config& sigma) const
{
   sigma.resize(V);
   for(int ispin=0; ispin<V; ispin++) 
   {
      sigma[ispin] = +1;
      if( (ispin%L)>(L/2) )
         sigma[ispin] *= stagmask[ispin];
   }
}


void MfiaMF_Hamiltonian::initial_ferro(MfiaMF_Hamiltonian::Config& sigma) const
{
   sigma.resize(V);
   for(int ispin=0; ispin<V; ispin++) sigma[ispin] = +1;
}

#endif	// MFIA_HAMILTONIAN_HPP

