#ifndef HEISENBERG_HAMILTONIAN_HPP_
#define HEISENBERG_HAMILTONIAN_HPP_

//#include"R1279.hpp"
//#include"SeedFromClock.hpp"

#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#include<map>
#include<algorithm>
#include<limits>

// Implementation of the extended Heisenberg model
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// August 5, 2015 changed M[3] -> M3d[3], made M scalar
//                This is for compatibility with Ising/Lattice Gas codes
// August 15, 2014 Initial Version

struct Heisenberg_Observables {
public:
    double E,X;   // Required: Energy of the walker and order parameter
    double E_N;   // Nearest neighbor energy
    double E_K;   // Anisotropy energy
    double E_D;   // Demagnetization energy
    double E_H;   // Zeeman energy
    double M3d[3];// Global magnetization
    double M;     // Scalar magnetization
    double Mmag2;  
    double Mstag; // Staggered magnetization
    double V;     // Volume = Number of spins
public:
    Heisenberg_Observables() { E=0; X=0; E_N=E_K=E_D=E_H=0; V=0; } 
    Heisenberg_Observables(const Heisenberg_Observables& w) { this->copy(w); }
    void operator=(const Heisenberg_Observables& w) { this->copy(w); }
    void copy(const Heisenberg_Observables& w)
    {
       E = w.E;
       E_N = w.E_N;
       E_K = w.E_K;
       E_D = w.E_D;
       E_H = w.E_H;
       X = w.X;
       M = w.M;
       Mmag2 = w.Mmag2;
       Mstag = w.Mstag;
       V = w.V;
       for(int i=0; i<3; i++) M3d[i] = w.M3d[i];
    }
    bool operator==(const Heisenberg_Observables& w) const
    {
#      if 0
          return ( E==w.E && X==w.X && Mstag==w.Mstag );
#      else
          using std::fabs;
          const double eps = 1.e-6;
          //bool pass = ( fabs(E-w.E)<eps && fabs(X-w.X)<eps && fabs(Mstag-w.Mstag)<eps );
          bool pass = ( fabs((E-w.E)/E)<eps && fabs((X-w.X)/X)<eps && fabs((Mstag-w.Mstag)/Mstag)<eps );
          if( false && !pass )
          {
             std::cout << __FILE__ << ":" << __LINE__ << " Observables not equal" << std::endl;
             std::cout << "pass = " << pass << std::endl;
             std::cout << "E: " << ( fabs((E-w.E)/E)<eps ) << " " << E << " " << w.E << " " << fabs((E-w.E)/E) << std::endl;
             std::cout << "X: " << ( fabs((X-w.X)/X)<eps ) << " " << X << " " << w.X << " " << fabs((X-w.X)/X) << std::endl;
             std::cout << "Stag: " << ( fabs((Mstag-w.Mstag)/Mstag)<eps ) << " " << Mstag << " " << w.Mstag << " " << fabs((Mstag-w.Mstag)/Mstag) << std::endl;
          }
          return pass;
#         endif
    }
    bool operator!=(const Heisenberg_Observables& w) const { return !((*this)==w); }
};



class Heisenberg_Hamiltonian 
{
public:

   typedef Heisenberg_Observables   Observables;
   typedef std::vector<double>      Config;

public:

   enum { DIM = 3 };               // Lattice dimension
   int L;                          // Side of the lattice
   int nspin;                      // Number of spins
   std::vector<int> jspin;         // List of neighbors
   std::vector<double> stagmask;   // Encodes sublattice for staggered magnetization
 
   double H;                       // The Zeeman field

public:

   // Part of the MC Hamiltonian interface
   double step_size;               // 0<step_size<=2
   double MPI;                     // Value of pi used in step_along spin

   template<typename URNG>
   void step_along_spin(double cartesian_spin[3], double step_size, URNG& urng) const;

   // Calculate all macroscopic variables based on this sigma
   virtual void calc_observable(Config& sigma, Heisenberg_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   virtual void change(Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[4]);
    
   // Make an attempted Monte Carlo step
   template<typename MCWalker, typename URNG> void mc_step(MCWalker& walker, URNG& urng);

public:

   virtual std::string header() const 
   { 
      char buffer[1024];
      sprintf(buffer,"# Simple Heisenberg Model H=%lf\n# Dim=%d L=%d N=%d\n",H,DIM,L,nspin);
      return std::string(buffer);
   }

   // Default constructor
   Heisenberg_Hamiltonian(int Lval=8) 
   {
      step_size = 1;
      H = 0;
      set_size(Lval);
   }
   

   // Construct a density of states g(E) on range [Emin,Emax] using polynomial coefficients
   void set_size(int Lval=8) 
   {
      // Extent in each dimension
      L = Lval;
      this->init();
   }

   virtual void init(bool verbose=false);

   void initial_mixed(Config& sigma) const;
   void initial_ferro(Config& sigma) const;
   template<typename URNG> void initial_random(Config& sigma, URNG& urng) const;

};



void Heisenberg_Hamiltonian::init(bool verbose) 
{
   MPI = 4.*std::atan(1.);
   int lda[DIM+1];
   lda[0] = 1;
   for(int idim=0; idim<DIM; idim++)    lda[idim+1] = L*lda[idim];
   nspin = lda[DIM];
   stagmask.resize(nspin);
   // neighbor table
   jspin.resize(2*DIM*nspin);
   for(int ispin=0; ispin<nspin; ispin++) {
      // coordinates on lattice
      int icoord[DIM];
      int isite = ispin;
      for(int idim=DIM-1; idim>=0; idim--) {
         icoord[idim] = isite/lda[idim];
         isite -= icoord[idim]*lda[idim];
      }
      // forward neighbors
      for(int idim=0; idim<DIM; idim++) {
         int jcoord[DIM];
         for(int k=0; k<DIM; k++) jcoord[k]=icoord[k];
         jcoord[idim]++;
         if( jcoord[idim]>=L ) jcoord[idim]=0;
         int jsite = 0;
         for(int k=0; k<DIM; k++) jsite += jcoord[k]*lda[k];
         jspin[2*DIM*ispin+idim] = jsite;
      }
      // backward neighbors
      for(int idim=0; idim<DIM; idim++) {
         int jcoord[DIM];
         for(int k=0; k<DIM; k++) jcoord[k]=icoord[k];
         jcoord[idim]--;
         if( jcoord[idim]<0 ) jcoord[idim]=L-1;
         int jsite = 0;
         for(int k=0; k<DIM; k++) jsite += jcoord[k]*lda[k];
         jspin[2*DIM*ispin+DIM+idim] = jsite;
      }
      // staggered order parameter weights
      int csum = 0;
      for(int idim=0; idim<DIM; idim++) csum += icoord[idim];
      int odd = csum%2;
      if( odd==1 )
         stagmask[ispin] = -1;
      else
         stagmask[ispin] = +1;
   } 
}

 

void Heisenberg_Hamiltonian::calc_observable(Heisenberg_Hamiltonian::Config& sigma, Heisenberg_Observables& macro)
{
   macro.V = nspin;
   double* S(&sigma[0]);            // Hiesenberg specific
   if( false )
   {
      // Check validity of sigma
      if( (sigma.size()/3) != nspin)
         std::cout << __FILE__ << ":" << __LINE__ << " Bad sigma size=" << sigma.size() << "," << sigma.size()/3 << " nspin=" << nspin << std::endl;
      for(int ispin=0; ispin<nspin; ispin++)
      {
         double mag2 = S[3*ispin+0]*S[3*ispin+0]+S[3*ispin+1]*S[3*ispin+1]+S[3*ispin+2]*S[3*ispin+2];
         if( std::fabs(mag2-1.0)>1.e-2 )
             std::cout << __FILE__ << ":" << __LINE__ << " Bad spin ispin=" << ispin << " mag2=" << mag2 << " (" << S[0] << "," << S[1] << "," << S[2] << ")" << std::endl;
      }
   }
   double net = 0;
   double mag[3] = { 0, 0, 0 };
   for(int ispin=0; ispin<nspin; ispin++) {
      for(int j=0; j<3; j++) mag[j] += S[3*ispin+j];
      for(int k=0; k<DIM; k++) 
      {
            net += S[3*(ispin)+0]*S[3*jspin[2*DIM*ispin+k]+0] 
                 + S[3*(ispin)+1]*S[3*jspin[2*DIM*ispin+k]+1] 
                 + S[3*(ispin)+2]*S[3*jspin[2*DIM*ispin+k]+2];
      }
   }
   double stag = 0;
   for(int ispin=0; ispin<nspin; ispin++) 
      stag += stagmask[ispin]*S[3*ispin+2]; // stagmask[ispin]*( mag[0]*S[3*ispin+0] + mag[1]*S[3*ispin+1] + mag[2]*S[3*ispin+2] );
   macro.E_N = -net;
   macro.E_H = -mag[2]*H;
   for(int k=0; k<3; k++) macro.M3d[k] = mag[k];
   macro.Mmag2 = macro.M3d[0]*macro.M3d[0] + macro.M3d[1]*macro.M3d[1] + macro.M3d[2]*macro.M3d[2];
   macro.Mstag = stag;
   macro.E = macro.E_N + macro.E_H;
   macro.X = std::sqrt(macro.Mmag2);
   macro.M = (H==0)? macro.X : macro.M3d[2];
}


template<typename URNG>
void Heisenberg_Hamiltonian::step_along_spin(double cartesian_spin[3], double step_size, URNG& urng) const
{
   if( false && (step_size>1 || step_size<0) )
   {
      std::cerr << "step_along_spin: invalid step_size = " << step_size << std::endl;
   }
   double tmp_spin[3]; 
   for(int k=0; k<3; k++) 
      tmp_spin[k] = cartesian_spin[k]; 
   // anonymous scope so can declare z,r below
   { 
      double  z = (1.0 - 2.0*step_size*urng());
      double  t = 2.0 * MPI * urng();
      double  r = std::sqrt(1.-z*z);
      cartesian_spin[0] = r*std::cos(t);
      cartesian_spin[1] = r*std::sin(t);
      cartesian_spin[2] = z;
   }
   double r = tmp_spin[0]*tmp_spin[0]+tmp_spin[1]*tmp_spin[1];
   if(r<std::numeric_limits<double>::min())   // 1e-20)
   {
     if (tmp_spin[2]<0)
     {
        cartesian_spin[0] = -cartesian_spin[0];
        cartesian_spin[1] = -cartesian_spin[1];
        cartesian_spin[2] = -cartesian_spin[2];
     }
       return;
   }
   r = std::sqrt(r);
   double x,y,z;
   x = tmp_spin[0];
   y = tmp_spin[1];
   z = tmp_spin[2];
   double RI[3][3];
   RI[0][0] = x*z/r;
   RI[1][0] = -y/r;
   RI[2][0] = x;
   RI[0][1] = y*z/r;
   RI[1][1] = x/r;
   RI[2][1] = y;
   RI[0][2] = -r;
   RI[1][2] = 0;
   RI[2][2] = z;
   // Change this to Mat-Vec multiply
   double tmp[3];
   for(int i=0;i<3;i++) {
      tmp[i] = 0;
      for(int j=0;j<3;j++) {
         tmp[i] += RI[j][i]*cartesian_spin[j];
      }
   }
   for(int k=0; k<3; k++) cartesian_spin[k] = tmp[k];
   // check conservation of spin length
   if( false ) {
      double rho2 = cartesian_spin[0]*cartesian_spin[0]
                  + cartesian_spin[1]*cartesian_spin[1]
                  + cartesian_spin[2]*cartesian_spin[2];
      double diff = std::fabs(rho2-1);
      if( diff>1.e-6) {
         std::cerr << "step_along_spin non-unit spin" << std::endl
                   << cartesian_spin[0] << ", "
                   << cartesian_spin[1] << ", "
                   << cartesian_spin[2] << std::endl
                   << "norm = " << rho2 << " diff = " << diff << std::endl;
      }
   }
}


// Make an attempted Monte Carlo step
template<typename MCWalker, typename URNG>
void Heisenberg_Hamiltonian::mc_step(MCWalker& walker, URNG& urng)
{
   // Randomly choose a move
   int NSite = walker.sigma.size()/3;
   int isite = static_cast<int>(NSite*urng());
   while( isite<0 || isite>=NSite )
   {
      //std::cerr << __FILE__ << ":" << __LINE__ << " markov_propose bad isite=" << isite << "/" << NSite << std::endl;
      isite = static_cast<int>(NSite*urng());
   }
   double new_sigma[3]; for(int i=0; i<3; i++) new_sigma[i] = walker.sigma[3*isite+i];
   step_along_spin(new_sigma,this->step_size,urng);
   // Do the move
   walker.save_initial();
   for(int k=0; k<3; k++) walker.add_change(3*isite+k);
   this->change(walker.sigma,walker.now,isite,new_sigma);
}


void Heisenberg_Hamiltonian::change(Heisenberg_Hamiltonian::Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[3])
{
   double* S(&sigma[0]);
   // Calculate before
   double E0 = 0;
   for(int k=0; k<2*DIM; k++)
      E0 += -S[3*ispin+0]*S[3*jspin[2*DIM*ispin+k]+0] 
          + -S[3*ispin+1]*S[3*jspin[2*DIM*ispin+k]+1] 
          + -S[3*ispin+2]*S[3*jspin[2*DIM*ispin+k]+2];
   macro.Mstag -= stagmask[ispin]*S[3*ispin+2]; // stagmask[ispin]*( macro.M[0]*S[3*ispin+0] + macro.M[1]*S[3*ispin+1] + macro.M[2]*S[3*ispin+2] );
   for(int k=0; k<3; k++)  macro.M3d[k] -= S[3*ispin+k];
   // Make move
   for(int k=0; k<3; k++) S[3*ispin+k] = new_sigma[k];
   // Calculate after
   for(int k=0; k<3; k++)  macro.M3d[k] += S[3*ispin+k];
   macro.Mmag2 = macro.M3d[0]*macro.M3d[0] + macro.M3d[1]*macro.M3d[1] + macro.M3d[2]*macro.M3d[2];
   macro.Mstag += stagmask[ispin]*S[3*ispin+2]; // stagmask[ispin]*( macro.M[0]*S[3*ispin+0] + macro.M[1]*S[3*ispin+1] + macro.M[2]*S[3*ispin+2] );
   double E1 = 0;
   for(int k=0; k<2*DIM; k++)
      E1 += -S[3*ispin+0]*S[3*jspin[2*DIM*ispin+k]+0] 
          + -S[3*ispin+1]*S[3*jspin[2*DIM*ispin+k]+1] 
          + -S[3*ispin+2]*S[3*jspin[2*DIM*ispin+k]+2];
   // update
   macro.E_N = macro.E_N + (E1-E0);
   macro.E_H = -macro.M3d[2]*H;
   macro.E   = macro.E_N + macro.E_H;
   macro.X   = std::sqrt(macro.Mmag2);
   macro.M   = (H==0)? macro.X : macro.M3d[2];
   // check energy
   if( false ) {
      static int icount = 0;
      Heisenberg_Observables direct(macro);
      calc_observable(sigma,direct);
      if( macro!=direct )
      {
         if( icount<100 ) 
         {
            icount++;
            std::cout << __FILE__ << ":" << __LINE__ << " change E != direct E (icount=" << icount << ")" << std::endl;
            std::cout << "E     = " << macro.E << ", " << direct.E << " diff=" << macro.E - direct.E << std::endl;
            std::cout << "E_N   = " << macro.E_N << ", " << direct.E_N << " diff=" << macro.E_N - direct.E_N << std::endl;
            std::cout << "E_H   = " << macro.E_H << ", " << direct.E_H << " diff=" << macro.E_H - direct.E_H << std::endl;
            std::cout << "Mmag2 = " << macro.Mmag2 << ", " << direct.Mmag2 << " diff=" << macro.Mmag2 - direct.Mmag2 << std::endl;
            std::cout << "Mstag = " << macro.Mstag << ", " << direct.Mstag << " diff=" << macro.Mstag - direct.Mstag << std::endl;
            std::cout << "M = (" << macro.M3d[0] << "," << macro.M3d[1] << "," << macro.M3d[2] << "),"
                      << "    (" << direct.M3d[0] << "," << direct.M3d[1] << "," << direct.M3d[2] << ")" << std::endl;
         }
         macro = direct;
      }
   }
}


void Heisenberg_Hamiltonian::initial_mixed(Heisenberg_Hamiltonian::Config& sigma) const
{
   this->initial_ferro(sigma);
   for(int ispin=0; ispin<nspin; ispin++)
   {
      if( (ispin%L)>(L/2) )
         sigma[3*ispin+2] *= stagmask[ispin];
   }
}


void Heisenberg_Hamiltonian::initial_ferro(Heisenberg_Hamiltonian::Config& sigma) const
{
   sigma.resize(3*nspin);
   for(int ispin=0; ispin<nspin; ispin++) 
   {
      sigma[3*ispin+0] = 0;
      sigma[3*ispin+1] = 0;
      sigma[3*ispin+2] = 1;
   }
}

template<typename URNG>
void Heisenberg_Hamiltonian::initial_random(Heisenberg_Hamiltonian::Config& sigma, URNG& urng) const
{
   sigma.resize(3*nspin);
   for(int ispin=0; ispin<nspin; ispin++) 
   {
      sigma[3*ispin+0] = 0;
      sigma[3*ispin+1] = 0;
      sigma[3*ispin+2] = 1;
      step_along_spin(&(sigma[3*ispin]),1,urng);
   }
}

#endif
