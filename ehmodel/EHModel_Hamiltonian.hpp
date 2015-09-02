#ifndef EHMODEL_HAMILTONIAN_HPP_
#define EHMODEL_HAMILTONIAN_HPP_

// Implementation of the extended Heisenberg model
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// August 5, 2015 -- Modified M[]->M3d to mirror Heisenberg_Observables
// September 25, 2014 -- Initial version


// include Heisenberg_Observables
#include"Heisenberg_Hamiltonian.hpp"

#include<iostream>
#include<string>
#include<cstring>
#include<cmath>
#include<vector>
#include<map>
#include<algorithm>
#include<limits>
#include<unistd.h>


class EHModel_Hamiltonian : public Heisenberg_Hamiltonian
{
public:

   typedef Heisenberg_Hamiltonian       BASE;
   typedef typename BASE::Observables   Observables;
   typedef typename BASE::Config        Config;


public:

   int nngbr;                      // Number of neighbors per spin
   std::vector<double> KiJij;      // Array of exchange constants and anisotropy information
   std::vector<double> moment;     // Magnetic moment of each atom

   // Calculate all macroscopic variables based on this sigma
   void calc_observable(Config& sigma, Heisenberg_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   void change(Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[3]);

public:

   std::string header() const 
   { 
      char buffer[1024];
      sprintf(buffer,"# Extended Heisenberg Model H=%lf\n# Dim=%d L=%d N=%d\n",H,DIM,L,nspin);
      return std::string(buffer);
   }

   EHModel_Hamiltonian(int Lval=8);
   
   void init(bool verbose=false);

   char model_fn[100];
   char model_ft[100];
   std::vector<std::string> filetypes;
   void read_pairs(std::string filename);
   void write_KiJij(std::string filename);

};


EHModel_Hamiltonian::EHModel_Hamiltonian(int Lval)
{
   L = Lval; 
   std::strcpy(model_fn,"EHModel-input.txt");
   filetypes.resize(2);
   filetypes[0] = "EHModelPairs";
   filetypes[1] = "EHModelKiJij";
   std::strcpy(model_ft,filetypes[0].c_str());
}




void EHModel_Hamiltonian::init(bool verbose) 
{
   // A simple Heisenberg model
   MPI = 4.*std::atan(1.);
   int lda[DIM+1];
   lda[0] = 1;
   for(int idim=0; idim<DIM; idim++)    lda[idim+1] = L*lda[idim];
   nspin = lda[DIM];
   stagmask.resize(nspin); 
   for(int ispin=0;  ispin<nspin; ispin++) stagmask[ispin]=1;
   bool haveread = false;
   if(true)
   {
      if( access( model_fn, F_OK ) != -1 )
      {
         int icase = 0; while(icase<filetypes.size() && filetypes[icase].compare(model_ft)!=0) icase++;
         haveread = true;
         switch (icase)
         {
         case 0:
            read_pairs(model_fn);
            break;
         case 1:
            break;
         default:
            haveread = false;
         }
      }
      else
      {
         std::cout << "Cannot open file \"" << model_fn << "\"" << std::endl;
      }
   }
   if( !haveread )
   {
      // Backup is simple Heisenberg model
      std::cout << "Using simple Heisenberg model for Hamiltonian" << std::endl;
      this->BASE::init(verbose);
      nngbr = 2*BASE::DIM;
      KiJij.resize( nspin*(nngbr+4) );
      for(int ispin=0; ispin<nspin; ispin++)
      {
         KiJij[ispin*(nngbr+4)+0] = 0;  // Khat_x
         KiJij[ispin*(nngbr+4)+1] = 0;
         KiJij[ispin*(nngbr+4)+2] = 0;
         KiJij[ispin*(nngbr+4)+3] = 0;  // |Ki|
         for(int j=0; j<nngbr; j++)
            KiJij[ispin*(nngbr+4)+4+j] = 1;
      }
   } 
}


void EHModel_Hamiltonian::calc_observable(Config& sigma, Heisenberg_Observables& macro)
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
   double E_N = 0;
   double E_K = 0;
   double mag[3] = { 0, 0, 0 };
   double stag = 0;
   int iptr = 0;
   for(int ispin=0; ispin<nspin; ispin++) 
   {
      // increment global magnetization
      for(int j=0; j<3; j++) mag[j] += moment[ispin]*S[3*ispin+j];
      stag += stagmask[ispin]*moment[ispin]*S[3*ispin+2];
      // calculate anisotropy energy KiJij[kx,ky,kz,Kmag]
      double sdotk = 0;
      for(int j=0; j<3; j++) sdotk += S[3*ispin+j]*KiJij[iptr++];
      E_K += KiJij[iptr++]*sdotk;
      // calculate exchange energy
      for(int j=0; j<nngbr; j++) 
      {
            E_N += KiJij[iptr++]
                * ( S[3*(ispin)+0]*S[3*jspin[nngbr*ispin+j]+0] 
                +   S[3*(ispin)+1]*S[3*jspin[nngbr*ispin+j]+1] 
                +   S[3*(ispin)+2]*S[3*jspin[nngbr*ispin+j]+2] );
      }
   }
   E_N /= 2.;  // Double counted two-body interactions
   macro.E_N = -E_N;
   macro.E_K = -E_K;
   macro.E_H = -mag[2]*H;
   for(int k=0; k<3; k++) macro.M3d[k] = mag[k];
   macro.Mmag2 = macro.M3d[0]*macro.M3d[0] + macro.M3d[1]*macro.M3d[1] + macro.M3d[2]*macro.M3d[2];
   macro.Mstag = stag;
   macro.E = macro.E_N + macro.E_K + macro.E_H;
   macro.X = std::sqrt(macro.Mmag2);
   macro.M = (H==0)? macro.X : macro.M3d[2];
}



void EHModel_Hamiltonian::change(EHModel_Hamiltonian::Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[3])
{
   double* S(&sigma[0]);
   // Calculate before
   int iptr = (nngbr+4)*ispin;
   double sdotk = 0;
   for(int j=0; j<3; j++) sdotk += S[3*ispin+j]*KiJij[iptr++];
   double E_K0 = KiJij[iptr++]*sdotk;
   double E_N0 = 0;
   for(int j=0; j<nngbr; j++) 
   {
         E_N0 += KiJij[iptr++]
             * ( S[3*(ispin)+0]*S[3*jspin[nngbr*ispin+j]+0] 
             +   S[3*(ispin)+1]*S[3*jspin[nngbr*ispin+j]+1] 
             +   S[3*(ispin)+2]*S[3*jspin[nngbr*ispin+j]+2] );
   }
   macro.Mstag -= stagmask[ispin]*moment[ispin]*S[3*ispin+2]; 
   for(int k=0; k<3; k++)  macro.M3d[k] -= moment[ispin]*S[3*ispin+k];
   // Make move
   for(int k=0; k<3; k++) S[3*ispin+k] = new_sigma[k];
   // Calculate after
   for(int k=0; k<3; k++)  macro.M3d[k] += moment[ispin]*S[3*ispin+k];
   macro.Mmag2 = macro.M3d[0]*macro.M3d[0] + macro.M3d[1]*macro.M3d[1] + macro.M3d[2]*macro.M3d[2];
   iptr = (nngbr+4)*ispin;
   sdotk = 0;
   for(int j=0; j<3; j++) sdotk += S[3*ispin+j]*KiJij[iptr++];
   double E_K1 = KiJij[iptr++]*sdotk;
   double E_N1 = 0;
   for(int j=0; j<nngbr; j++) 
   {
         E_N1 += KiJij[iptr++]
             * ( S[3*(ispin)+0]*S[3*jspin[nngbr*ispin+j]+0] 
             +   S[3*(ispin)+1]*S[3*jspin[nngbr*ispin+j]+1] 
             +   S[3*(ispin)+2]*S[3*jspin[nngbr*ispin+j]+2] );
   }
   macro.Mstag += stagmask[ispin]*S[3*ispin+2];
   // update
   macro.E_N = macro.E_N + (-1)*(E_N1-E_N0);   // minus sign in Hamiltonian
   macro.E_K = macro.E_K + (-1)*(E_K1-E_K0);
   macro.E_H = -macro.M3d[2]*H;
   macro.E   = macro.E_N + macro.E_K + macro.E_H;
   macro.X   = std::sqrt(macro.Mmag2);
   macro.M   = (H==0)? macro.X : macro.M3d[2];
   // check energy
   if( true ) {
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
            std::cout << "E_K   = " << macro.E_N << ", " << direct.E_N << " diff=" << macro.E_N - direct.E_N << std::endl;
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


void EHModel_Hamiltonian::read_pairs(std::string filename)
{
   std::ifstream fin(filename.c_str());
   std::string buffer;
   const int nheader = 3;
   for(int i=0; i<nheader; i++)
      std::getline(fin,buffer); 
   fin >> nspin;
   if( nspin<1 )
   {
      std::cerr << __FILE__ << ":" << __LINE__ << " EHModel_Hamiltonian::read_pairs" << std::endl
               << "invalid number of spins. nspin=" << nspin << std::endl;
      return;
   }
   for(int ispin=0; ispin<nspin; ispin++)
   {
      int index, npair;
      fin >> index >> npair;
      if( ispin==0 )
      {
         this->nngbr = npair;
         if( nngbr<0 )
         {
            std::cerr << __FILE__ << ":" << __LINE__ << " EHModel_Hamiltonian::read_pairs" << std::endl
                      << "invalid number of pairs for first spin. npair=" << npair << std::endl;
            return;
         }
         KiJij.resize( nspin*(nngbr+4) ); std::fill(KiJij.begin(),KiJij.end(),0.);
         jspin.resize( nspin*nngbr );     std::fill(jspin.begin(),jspin.end(),-1);
      }
      if( npair<0 || npair>nngbr )
      {
         std::cerr << __FILE__ << ":" << __LINE__ << " EHModel_Hamiltonian::read_pairs" << std::endl
                   << "assumes first spin has the maximum number of neighbors." << std::endl
                   << "ispin=" << ispin << " npair=" << npair << " EHModel.nngbr=" << nngbr << std::endl;
         return;
      }
      for(int ipair=0; ipair<npair; ipair++)
      {
         int rank;
         fin >> jspin[index*nngbr+ipair] >> rank >> KiJij[index*(nngbr+4)+4+ipair];
         if( rank!=1 )
         {
            std::cerr << __FILE__ << ":" << __LINE__ << " EHModel_Hamiltonian::read_pairs" << std::endl
                      << "Assumes rank=1 for all Jij" << std::endl;
            return;
         }
      }
   }
   moment.resize(nspin);
   for(int ispin=0; ispin<nspin; ispin++) moment[ispin] = 1;
   if( true ) write_KiJij("EHModel-KiJij.txt");
}


void EHModel_Hamiltonian::write_KiJij(std::string filename)
{
   std::ofstream fout(filename.c_str());
   fout << "# EHModel KiJij output" << std::endl;
   fout << "# nspin numneighbor" << std::endl;
   fout << "# Kx, Ky Kz, |K|; jspin Jij" << std::endl;
   fout << nspin << " " << nngbr << std::endl;
   for(int ispin=0; ispin<nspin; ispin++)
   {
      fout << KiJij[ispin*(nngbr+4)+0] << " " << KiJij[ispin*(nngbr+4)+1] << " " << KiJij[ispin*(nngbr+4)+2] << " " << KiJij[ispin*(nngbr+4)+3] << std::endl;
      for(int ingbr=0; ingbr<nngbr; ingbr++)
         fout << jspin[ispin*nngbr+ingbr] << " " << KiJij[ispin*(nngbr+4)+4+ingbr] << std::endl;
      
   }
}


#endif
