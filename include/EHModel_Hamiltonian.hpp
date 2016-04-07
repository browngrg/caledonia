#ifndef EHMODEL_HAMILTONIAN_HPP_
#define EHMODEL_HAMILTONIAN_HPP_

// Implementation of the extended Heisenberg model
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// August 5, 2015 -- Modified M[]->M3d to mirror Heisenberg_Observables
// September 25, 2014 -- Initial version


// include Heisenberg_Observables
#include"Heisenberg_Hamiltonian.hpp"

#include<iostream>
#include<fstream>
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

   // jspin[offset+nindex] comes from the base class
   enum { KRECSIZE = 9 };          // Storage space of anisotropy information (symmetry, K, K', dir1, dir2)
   int nngbr;                      // Number of neighbors per spin
   std::vector<double> KiJij;      // Array of exchange constants and anisotropy information
   std::vector<double> moment;     // Magnetic moment of each atom

   // Calculate all macroscopic variables based on this sigma
   void calc_observable(Config& sigma, Heisenberg_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   void change(Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[3]);

   void change_field(Heisenberg_Observables& macro, float Hnew) { BASE::change_field(macro,Hnew); }

public:

   std::string header() const 
   { 
      char buffer[1024];
      sprintf(buffer,"# Extended Heisenberg Model H=%lf\n# Dim=%d L=%d N=%d\n",H,DIM,L,nspin);
      return std::string(buffer);
   }

   EHModel_Hamiltonian(int Lval=8);
  
   template<typename OPTIONS> void add_options(OPTIONS& options);

   void init(bool verbose=false);

   char model_fn[100];
   char model_ft[100];
   std::vector<std::string> filetypes;
   void read_pairs(std::string filename);
   void write_KiJij(std::string filename);

   double calc_Ki(Config& sigma, int ispin);

   void map_energy(std::string filename="EHModelEnergyMap.csv");

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

template<typename OPTIONS>
void EHModel_Hamiltonian::add_options(OPTIONS& options)
{
   options.add_option("H",        "magnitude of magnetic field",    'H', &H );
   options.add_option("ehfile",   "EHModel input file",             ' ', model_fn );
   options.add_option("ehfilet",  "EHModel input file type",        ' ', model_ft );
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
      KiJij.resize( nspin*(nngbr+KRECSIZE) );
      for(int ispin=0; ispin<nspin; ispin++)
      {
         KiJij[ispin*(nngbr+KRECSIZE)+0] = 0;  // type = none
         KiJij[ispin*(nngbr+KRECSIZE)+1] = 0;  // K1
         KiJij[ispin*(nngbr+KRECSIZE)+2] = 0;  // K2
         KiJij[ispin*(nngbr+KRECSIZE)+3] = 0;  // Kdir1
         KiJij[ispin*(nngbr+KRECSIZE)+4] = 0;
         KiJij[ispin*(nngbr+KRECSIZE)+5] = 0;
         KiJij[ispin*(nngbr+KRECSIZE)+6] = 0;  // Kdir2
         KiJij[ispin*(nngbr+KRECSIZE)+7] = 0;
         KiJij[ispin*(nngbr+KRECSIZE)+8] = 0;
         for(int j=0; j<nngbr; j++)
            KiJij[ispin*(nngbr+KRECSIZE)+KRECSIZE+j] = 1;
      }
   } 
}

/* *  Calculate the magnetic anisotropy energy of site i
   *
   *  The energy for uniaxial anisotropy is given by
   *  \f[ E_{{\rm u},i} = K_i (1 - (\hat{e}_i \cdot \hat{s}_i)^2) \f]
   *  where <I>K<SUB>i</SUB></I> is the strength of the anisotropy,
   *  \f$\hat{e}_i\f$ is the axis of the anisotropy, and 
   *  \f$\hat{s}_i\f$ is the direction of the <I>i</I>-th spin.
   *
   *  The energy for square anisotropy is given by
   *  \f[ E_{{\rm s},i} = K_i ( (\hat{e}_{1,i} \cdot \hat{s}_i)^2
   *                           +(\hat{e}_{2,i} \cdot \hat{s}_i)^2 )
   *                     + K'_i (\hat{e}_{1,i} \cdot \hat{s}_i)^2
   *                            (\hat{e}_{2,i} \cdot \hat{s}_i)^2
   *  \f]
   *  where <I>K<SUB>i</SUB></I> and <I>K'<SUB>i</SUB></I> together 
   *  describe the strength of the anisotropy and 
   *  \f$hat{e}_{1,i}\f$ and \f$\hat{e}_{2,i}\f$ describe its axes. 
   *
   *  The energy for cubic anisotropy is given by
   *  \f[ E_{{\rm c},i} = K_i ( (\hat{e}_{1,i} \cdot \hat{s}_i)^2
   *                           +(\hat{e}_{2,i} \cdot \hat{s}_i)^2 
   *                           +(\hat{e}_{3,i} \cdot \hat{s}_i)^2 )
   *                     + K'_i (\hat{e}_{1,i} \cdot \hat{s}_i)^2
   *                            (\hat{e}_{2,i} \cdot \hat{s}_i)^2
   *                            (\hat{e}_{3,i} \cdot \hat{s}_i)^2
   *  \f]
   *  where <I>K<SUB>i</SUB></I> and <I>K'<SUB>i</SUB></I> together 
   *  describe the strength of the anisotropy and 
   *  \f$hat{e}_{1,i}\f$, \f$\hat{e}_{2,i}\f$, and \f$\hat{e}_{3,i}\f$ 
   *  describe its axes. It is assumed that the three axes are mutually
   *  orthogonal, so only two of them actually need to be specified. 
   *  The direction of the third can be calculated.
   */
double EHModel_Hamiltonian::calc_Ki(Config& sigma, int ispin)
{
   double Ek = 0;
   double a1sqr,a2sqr,a3sqr;
   double* Krec = &(KiJij[ispin*(nngbr+KRECSIZE)]);   // atype, K1, K2, Kdir1, Kdir2
   int atype = static_cast<int>(Krec[0]);
   switch(atype) {
   case 0:
      // no anisotropy
      break;
   case 1:
      // uniaxial anisotropy
      a1sqr = 0;
      for(int j=0; j<3; j++) a1sqr += sigma[3*ispin+j]*Krec[3+j];  // S*Kdir1
      a1sqr *= a1sqr;                                              // (S*Kdir1)^2
      Ek = Krec[1]*(1.-a1sqr);                                     // K*[1-(S*Kdir1)^2]
      break;
   case 2:
      // square anisotropy
      a1sqr = 0;
      for(int j=0; j<3; j++) a1sqr += sigma[3*ispin+j]*Krec[3+j];  // S*Kdir1
      a1sqr *= a1sqr;
      a2sqr = 0;
      for(int j=0; j<3; j++) a2sqr += sigma[3*ispin+j]*Krec[6+j];  // S*Kdir2
      a2sqr *= a2sqr;
      Ek = Krec[1]*(a1sqr+a2sqr) + Krec[2]*a1sqr*a2sqr;            // K1*[ (S*K1)^2 + (S*K2)^2 ] + K2*(S*K1)^2*(S*K2)^2
      break;
   case 3:
      // cubic anisotropy
      a1sqr = 0;
      for(int j=0; j<3; j++) a1sqr += sigma[3*ispin+j]*Krec[3+j];  // S*Kdir1
      a1sqr *= a1sqr;
      a2sqr = 0;
      for(int j=0; j<3; j++) a2sqr += sigma[3*ispin+j]*Krec[6+j];  // S*Kdir2
      a2sqr *= a2sqr;
      a3sqr = 1. - a1sqr - a2sqr;                                  // dot product of two unit vectors
      Ek = Krec[1]*(a1sqr*a2sqr+a2sqr*a3sqr+a3sqr*a1sqr) + Krec[2]*a1sqr*a2sqr*a3sqr;
      break;
   default:
      std::cout << __FILE__ << ":" << __LINE__ << " Error anisotropy type " << atype << " not recognized" << std::endl;
   }
   return Ek;
}


void EHModel_Hamiltonian::calc_observable(Config& sigma, Heisenberg_Observables& macro)
{
   if( !use_walker_H ) macro.H = H;
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
   for(int ispin=0; ispin<nspin; ispin++) 
   {
      // increment global magnetization
      for(int j=0; j<3; j++) mag[j] += moment[ispin]*S[3*ispin+j];
      stag += stagmask[ispin]*moment[ispin]*S[3*ispin+2];
      E_K += calc_Ki(sigma,ispin);
      // calculate exchange energy
      int iptr = (nngbr+KRECSIZE)*ispin + KRECSIZE;
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
   macro.E_K =  E_K;
   macro.E_H = -mag[2]*macro.H;
   for(int k=0; k<3; k++) macro.M3d[k] = mag[k];
   macro.Mmag2 = macro.M3d[0]*macro.M3d[0] + macro.M3d[1]*macro.M3d[1] + macro.M3d[2]*macro.M3d[2];
   macro.Mstag = stag;
   macro.E = macro.E_N + macro.E_K + macro.E_H;
   macro.X = std::sqrt(macro.Mmag2);
   macro.M = macro.M3d[2];
}


void EHModel_Hamiltonian::change(EHModel_Hamiltonian::Config& sigma, Heisenberg_Observables& macro, int ispin, double new_sigma[3])
{
   if( !use_walker_H ) macro.H = H;
   double* S(&sigma[0]);
   // Calculate before
   double E_K0 = calc_Ki(sigma,ispin);
   double E_N0 = 0;
   int iptr = (nngbr+KRECSIZE)*ispin + KRECSIZE;
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
   iptr = (nngbr+KRECSIZE)*ispin;
   double E_K1 = calc_Ki(sigma,ispin);
   double E_N1 = 0;
   iptr = (nngbr+KRECSIZE)*ispin + KRECSIZE;
   for(int j=0; j<nngbr; j++) 
   {
         E_N1 += KiJij[iptr++]
             * ( S[3*(ispin)+0]*S[3*jspin[nngbr*ispin+j]+0] 
             +   S[3*(ispin)+1]*S[3*jspin[nngbr*ispin+j]+1] 
             +   S[3*(ispin)+2]*S[3*jspin[nngbr*ispin+j]+2] );
   }
   macro.Mstag += stagmask[ispin]*moment[ispin]*S[3*ispin+2];
   // update
   macro.E_N = macro.E_N + (-1)*(E_N1-E_N0);   // minus sign in Hamiltonian
   macro.E_K = macro.E_K +      (E_K1-E_K0);
   macro.E_H = -macro.M3d[2]*macro.H;
   macro.E   = macro.E_N + macro.E_K + macro.E_H;
   macro.X   = std::sqrt(macro.Mmag2);
   macro.M   = macro.M3d[2];
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
            std::cout << "E_K   = " << macro.E_K << ", " << direct.E_K << " diff=" << macro.E_K - direct.E_K << std::endl;
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
         KiJij.resize( nspin*(nngbr+KRECSIZE) ); std::fill(KiJij.begin(),KiJij.end(),0.);
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
         fin >> jspin[index*nngbr+ipair] >> rank >> KiJij[index*(nngbr+KRECSIZE)+KRECSIZE+ipair];
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
   fout << "# moment; Ktype, K1, K2, K1x, K1y K1z, K2x, K2y, K2z; jspin Jij" << std::endl;
   fout << nspin << " " << nngbr << std::endl;
   for(int ispin=0; ispin<nspin; ispin++)
   {
      fout << moment[ispin] << std::endl;
      fout << KiJij[ispin*(nngbr+KRECSIZE)+0];
      for(int i=1; i<KRECSIZE; i++) fout << " " << KiJij[ispin*(nngbr+KRECSIZE)+i];
      fout << std::endl;
      for(int ingbr=0; ingbr<nngbr; ingbr++)
         fout << jspin[ispin*nngbr+ingbr] << " " << KiJij[ispin*(nngbr+KRECSIZE)+KRECSIZE+ingbr] << std::endl;
   }
}


void EHModel_Hamiltonian::map_energy(std::string filename)
{
   const int nz = 128+1;             // number of pts along z axis for each phi, needs to be odd to get z=0
   double pi = 4.*atan(1.);
   FILE* fout = fopen(filename.c_str(),"w");
   fprintf(fout,"# Energy for uniformly magnetized particle\n");
   fprintf(fout,"# Column 1: phi\n");
   fprintf(fout,"# Column 2: z = cos(theta)\n");
   fprintf(fout,"# Column 3: total energy\n");
   fprintf(fout,"# Column 4: Exchange energy\n");
   fprintf(fout,"# Column 5: Anisotropy energy\n");
   fprintf(fout,"# Column 6: Dipole-dipole energy\n");
   fprintf(fout,"# Column 7: Zeman energy\n");
   fprintf(fout,"# Column 8: global Mx\n");
   fprintf(fout,"# Column 9: global My\n");
   fprintf(fout,"# Column 10: global Mz\n");
   Config sigma; this->initial(sigma);
   Observables macro;
   for(int iphi=0; iphi<5; iphi++)
   {
      double phi = iphi*pi/4.;
      double cosp = std::cos(phi);
      double sinp = std::sin(phi);
      for(int iz=0; iz<nz; iz++)
      {
         double z = 1.-2.*static_cast<double>(iz)/static_cast<double>(nz-1);
         double r = std::sqrt(1.-z*z);
         double x = r*cosp;
         double y = r*sinp;
         for(int ispin=0; ispin<nspin; ispin++)
         {
            sigma[3*ispin+0] = x;
            sigma[3*ispin+1] = y;
            sigma[3*ispin+2] = z;
         }
         this->calc_observable(sigma,macro);
         fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",phi,z,macro.E,macro.E_N,macro.E_K,macro.E_D,macro.E_H,macro.M3d[0],macro.M3d[1],macro.M3d[2]);
      }
   }
   fclose(fout);
}

#endif
