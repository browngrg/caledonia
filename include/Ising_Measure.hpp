#ifndef ISING_MEASURE_HPP
#define ISING_MEASURE_HPP


#include<string>
#include<stdio.h>
#include<vector>
#include<map>
#include"MPI_Struct.hpp"
#include"StatRec.hpp"


class Ising_Measure
{
public:

   int         icall;    // number of calls
   bool        written;  // write since last sample?
   std::string header;   // Canned header information from other objects
   std::string filename;
   MPI_Struct  mp;       // MPI Information

   std::vector<float> kTlist;

private:

   enum { MAXK = 5 };
   std::vector<double> statbuff;
   double *emom,*mmom,*xmom,*mMagmom,*xMagmom;
   double V;

public:

   Ising_Measure();

   ~Ising_Measure();

   // Use the number of temperatures to initialize storage
   void init(const std::vector<float>& kT);

   // Add a string to the header information
   void add_header(std::string new_lines);

   // Add a sample state of walker to statistics (called by sampling object)
   template<typename MCWalker>
   void add_sample(const MCWalker& walker, bool at_equilibrium=true);

   // Clear accummulated data and start over
   void clear();

   // Write current statistics (could be a check-point)
   void write();

};



Ising_Measure::Ising_Measure() 
{ 
   written=true; 
   icall = 0;
   // Any automatic set-up actions go here
   filename = "";
   header="";
   V = 1;
}

Ising_Measure::~Ising_Measure() 
{ 
   if( !written ) this->write(); 
   // Any clean-up actions go here:
}

void Ising_Measure::add_header(std::string new_lines)
{
   header = header + new_lines;
}

template<typename MCWalker>
void Ising_Measure::add_sample(const MCWalker& walker, bool at_equilibrium)
{
   icall++;
   written = false;
   // Actions that happen regardless of equilibrium:
   V = walker.now.V;
   // Equilibrium Measurements Come Here:
   if( !at_equilibrium) return;
   int ibin = walker.kTi;            // Defined for RE_Walker
   double xval = walker.now.E;
   double xk = 1;
   for(int k=0; k<MAXK; k++) 
   {
      emom[ibin*MAXK+k] += xk;
      xk *= xval;
   }
   xval = walker.now.M;
   xk = 1;
   for(int k=0; k<MAXK; k++) 
   {
      mmom[ibin*MAXK+k] += xk;
      xk *= xval;
   }
   if(xval<0) xval*=-1;
   mMagmom[ibin] += xval;
   xval = walker.now.X;
   xk = 1;
   for(int k=0; k<MAXK; k++) 
   {
      xmom[ibin*MAXK+k] += xk;
      xk *= xval;
   }
   if(xval<0) xval*=-1;
   xMagmom[ibin] += xval;
}

void Ising_Measure::clear()
{
   written = false;
   icall = 0;
   // Actions to clear accumulated statistics go here:
   for(int i=0; i<statbuff.size(); i++) statbuff[i]=0;
}

void Ising_Measure::write()
{
   // mpi_gather()
   int nbuff = statbuff.size();
   std::vector<double> stat_global(nbuff,0);
#  ifdef USE_MPI
   if(mp.in()) MPI_Allreduce(&(statbuff[0]),&(stat_global[0]),nbuff,MPI_DOUBLE,MPI_SUM,mp.comm);
#  else
   mp.iproc=0;
   stat_global = statbuff;
#  endif
   if( mp.iproc==0 )
   {
      int NBin = kTlist.size();
      double *emom,*mmom,*xmom,*mMagmom,*xMagmom;
      emom = &(stat_global[0*MAXK*NBin]);
      mmom = &(stat_global[1*MAXK*NBin]);
      xmom = &(stat_global[2*MAXK*NBin]);
      mMagmom = &(stat_global[3*MAXK*NBin+0*NBin]);
      xMagmom = &(stat_global[3*MAXK*NBin+1*NBin]);
      std::ofstream fout("Ising.csv");
      fout << "# Ising statistics" << std::endl;
      fout << header;
      fout << "# Column 1: T temperature in energy units" << std::endl;
      fout << "# Column 2: <|M|>/N" << std::endl;
      fout << "# Column 3: <M>/N" << std::endl;
      fout << "# Column 4: <M^2>/N" << std::endl;
      fout << "# Column 5: <M^2>-<M>^2)" << std::endl;
      fout << "# Column 6: chi_m/N = (<M^2>-<|M|>^2)/T/N" << std::endl;
      fout << "# Column 7: Binder cummulant M" << std::endl;
      fout << "# Column 8: <E>/N Energy per site" << std::endl;
      fout << "# Column 9: <E^2>-<E>^2" << std::endl;
      fout << "# Column 10: C_H/N = (<E^2>-<E>^2)/T/T/N" << std::endl;
      fout << "# Column 11: <|X|>/N  (Order Parameter defined by Hamiltonian/walker)" << std::endl;
      fout << "# Column 12: <X>/N" << std::endl;
      fout << "# Column 13: <X^2>/N" << std::endl;
      fout << "# Column 14: <X^2>-<X>^2" << std::endl;
      fout << "# Column 15: Binder cummulant X" << std::endl;
      fout << "# Column 16: nsamp" << std::endl;
      std::vector<double> momv(MAXK);
      for(int ibin=0; ibin<NBin; ibin++)
      {
         double kT = kTlist[ibin];
         int nsamp = emom[ibin*MAXK];
         for(int i=0; i<MAXK; i++) momv[i]=emom[ibin*MAXK+i];
         if( momv[0]!=nsamp ) std::cout << __FILE__ << ":" << __LINE__ << " sampling error" << std::endl;
         StatRec stats(momv);
         double Emean = stats.mu;
         double Esq = stats.mup_2;                                 // read mu_2 prime, the noncentral moment
         double Efluc = stats.mu_2;                                // mu_2, the central moment
         double C_H = Efluc/kT/kT/V;
         for(int i=0; i<MAXK; i++) momv[i]=mmom[ibin*MAXK+i];
         if( momv[0]!=nsamp ) std::cout << __FILE__ << ":" << __LINE__ << " sampling error" << std::endl;
         stats(momv);
         double Mmag = mMagmom[ibin]/mmom[ibin*MAXK];
         if( momv[0]!=nsamp ) std::cout << __FILE__ << ":" << __LINE__ << " sampling error" << std::endl;
         double Mmean = stats.mu;
         double Msq = stats.mup_2;
         double Mfluc = stats.mu_2;
         double chim = (Msq-Mmag*Mmag)/kT/V;
         double Mcumm = stats.binder_4;
         for(int i=0; i<MAXK; i++) momv[i]=xmom[ibin*MAXK+i];
         if( momv[0]!=nsamp ) std::cout << __FILE__ << ":" << __LINE__ << " sampling error" << std::endl;
         stats(momv);
         double Xmag = xMagmom[ibin]/xmom[ibin*MAXK];
         double Xmean = stats.mu;
         double Xsq = stats.mup_2;
         double Xfluc = stats.mu_2;
         double Xcumm = stats.binder_4;
         fout << kT << " " << Mmag/V  << " " << Mmean/V << " " << Msq/V << " " << Mfluc << " " << chim  << " " << Mcumm << " "
                           << Emean/V << " " << Efluc   << " " << C_H   << " "
                           << Xmag/V  << " " << Xmean/V << " " << Xsq/V << " " << Xfluc << " " << Xcumm << " " << emom[ibin*MAXK] << std::endl;
      }
   }
   written = true;
}


void Ising_Measure::init(const std::vector<float>& kT)
{
   kTlist = kT;
   int NBin = kTlist.size();
   statbuff.resize( 3*MAXK*NBin+2*NBin);
   emom = &(statbuff[0*MAXK*NBin]);
   mmom = &(statbuff[1*MAXK*NBin]);
   xmom = &(statbuff[2*MAXK*NBin]);
   mMagmom = &(statbuff[3*MAXK*NBin+0*NBin]);
   xMagmom = &(statbuff[3*MAXK*NBin+1*NBin]);
   this->clear();
}



#endif
