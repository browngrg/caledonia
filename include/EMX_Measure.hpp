#ifndef EMX_MEASURE_HPP
#define EMX_MEASURE_HPP


#include<string>
#include<stdio.h>
#include<vector>
#include"MPI_Gang.hpp"
#include"StatRec.hpp"


class EMX_Measure
{
public:

   int    icall;         // number of calls
   bool   written;       // write since last sample?
   std::string header;   // Canned header information from other objects
   std::string filename;
   MPI_Gang    mp_window;
 
   float Elo,Ehi,Ebin;
   float Mlo,Mhi,Mbin;
   int NBinE,NBinM;

   bool make_map;

public:

   EMX_Measure();

   ~EMX_Measure();

   // Add a string to the header information
   void add_header(std::string new_lines);

   // Add a sample state of walker to statistics (called by sampling object)
   template<typename MCWalker>
   void add_sample(const MCWalker& walker, bool at_equilibrium=true);

   // Clear accummulated data and start over
   void clear();

   // Write current statistics (could be a check-point)
   void write();

   void init(float Emin, float Emax, float Ebin);

private:

   int bin(float E) const { int ibin=static_cast<int>((E-Elo)/Ebin); if(ibin<0) ibin=0; if(ibin>=NBinE) ibin=NBinE-1; return ibin; }
   int binM(float M) const { int ibin=static_cast<int>((M-Mlo)/Mbin); if(ibin<0) ibin=0; if(ibin>=NBinM) ibin=NBinM-1; return ibin; }

private:

   enum { MAXK = 5 };
   std::vector<double> statbuff;
   double *emom,*mmom,*xmom,*mMagmom,*xMagmom;
   std::vector<long long> em_map;

   
};



EMX_Measure::EMX_Measure() 
{ 
   written=true; 
   icall = 0;
   make_map = true;
   // Any automatic set-up actions go here
   filename = "";
}

EMX_Measure::~EMX_Measure() 
{ 
   if( !written ) this->write(); 
   // Any clean-up actions go here:
}

void EMX_Measure::add_header(std::string new_lines)
{
   header = header + new_lines;
}

template<typename MCWalker>
void EMX_Measure::add_sample(const MCWalker& walker, bool at_equilibrium)
{
   icall++;
   written = false;
   // Actions that happen regardless of equilibrium:
   int ibin = bin(walker.now.E);
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
   // Equilibrium Measurements Come Here:
   if( !at_equilibrium) return;
   int jbin = binM(static_cast<float>(walker.now.M)/static_cast<float>(walker.now.V)); 
   em_map[ibin*NBinM+jbin]++;
}

void EMX_Measure::clear()
{
   written = false;
   icall = 0;
   // Actions to clear accumulated statistics go here:
   for(int i=0; i<statbuff.size(); i++) statbuff[i]=0;
   for(int i=0; i<em_map.size(); i++) em_map[i]=0;
}

void EMX_Measure::write()
{
   // mpi_gather()
   int nbuff = statbuff.size();
   std::vector<double> stat_global(nbuff,0);
#  ifdef USE_MPI
   MPI_Allreduce(&(statbuff[0]),&(stat_global[0]),nbuff,MPI_DOUBLE,MPI_SUM,mp_window.comm_pool);
#  else
   em_global = em_map;
#  endif
   if( mp_window.iproc_gang==0 )
   {
      double *emom,*mmom,*xmom,*mMagmom,*xMagmom;
      emom = &(stat_global[0*MAXK*NBinE]);
      mmom = &(stat_global[1*MAXK*NBinE]);
      xmom = &(stat_global[2*MAXK*NBinE]);
      mMagmom = &(stat_global[3*MAXK*NBinE+0*NBinE]);
      xMagmom = &(stat_global[3*MAXK*NBinE+1*NBinE]);
      std::ofstream fout("EMX.csv");
      fout << "# EMX statistics" << std::endl;
      fout << header;
      fout << "# Column 1: Energy" << std::endl;
      fout << "# Column 2: <E^2>-<E>" << std::endl;
      fout << "# Column 3: <|M|>" << std::endl;
      fout << "# Column 4: <M>" << std::endl;
      fout << "# Column 5: <M^2>-<M>" << std::endl;
      fout << "# Column 6: Binder cummulant M" << std::endl;
      fout << "# Column 7: <|X|>  (Order Parameter defined by Hamiltonian/walker)" << std::endl;
      fout << "# Column 8: <X>" << std::endl;
      fout << "# Column 9: <X^2>-<X>" << std::endl;
      fout << "# Column 10: Binder cummulant X" << std::endl;
      std::vector<double> momv(MAXK);
      for(int ibin=0; ibin<NBinE; ibin++)
      {
         for(int i=0; i<MAXK; i++) momv[i]=emom[ibin*MAXK+i];
         StatRec stats(momv);
         double Emean = stats.mu;
         double Efluc = stats.mu_2;
         for(int i=0; i<MAXK; i++) momv[i]=mmom[ibin*MAXK+i];
         stats(momv);
         double Mmag = mMagmom[ibin]/mmom[ibin*MAXK];
         double Mmean = stats.mu;
         double Mfluc = stats.mu_2;
         double Mcumm = stats.binder_4;
         for(int i=0; i<MAXK; i++) momv[i]=xmom[ibin*MAXK+i];
         stats(momv);
         double Xmag = xMagmom[ibin]/xmom[ibin*MAXK];
         double Xmean = stats.mu;
         double Xfluc = stats.mu_2;
         double Xcumm = stats.binder_4;
         fout << Emean << " " << Efluc << " " 
              << Mmag  << " " << Mmean << " " << Mfluc << " " << Mcumm << " "
              << Xmag  << " " << Xmean << " " << Xfluc << " " << Xcumm << std::endl;
      }
   }
   if( em_map.size()>0 )
   {
      std::vector<long long> em_global(em_map.size());
#     ifdef USE_MPI
      if( mp_window.nproc_pool>1 )
      {
         int NBuff = em_map.size();
         MPI_Allreduce(&(em_map[0]),&(em_global[0]),NBuff,MPI_LONG_LONG,MPI_SUM,mp_window.comm_pool);
      }
      else
         em_global = em_map;
#     endif
      //
      if( mp_window.iproc_pool==0 )
      {
         std::ofstream fout("EMMap.csv");
         fout << "# EM Counts" << std::endl;
         fout << header;
         fout << "# Column 1: Energy" << std::endl;
         fout << "# Column 2: Magnetization" << std::endl;
         fout << "# Column 3: Row-normalized frequency" << std::endl;
         fout << "# Column 3: Raw Counts" << std::endl;
         for(int ibin=0; ibin<NBinE; ibin++)
         {
            float E = Elo + ibin*Ebin;
            long long rowN = 0;
            for(int jbin=0; jbin<NBinM; jbin++) rowN += em_global[ibin*NBinM+jbin];
            double rowf = 0; if(rowN>0) rowf = 1./static_cast<double>(rowN);
            for(int jbin=0; jbin<NBinM; jbin++)
            {
               float M = Mlo + jbin*Mbin;
               long long  raw = em_global[ibin*NBinM+jbin];
               double norm = rowf*raw;
               fout << E << " " << M << " " << norm << " " << raw << std::endl;
            }
         }
      }
   }
   written = true;
}

void EMX_Measure::init(float Emin, float Emax, float EBin)
{
   float Erange = Emax - Emin;
   Ebin = EBin;
   NBinE = Erange/Ebin;
   Elo = Emin;
   Ehi = Elo + Ebin*NBinE;
   statbuff.resize( 3*MAXK*NBinE + 2*NBinE );
   emom = &(statbuff[0*MAXK*NBinE]);
   mmom = &(statbuff[1*MAXK*NBinE]);
   xmom = &(statbuff[2*MAXK*NBinE]);
   mMagmom = &(statbuff[3*MAXK*NBinE+0*NBinE]);
   xMagmom = &(statbuff[3*MAXK*NBinE+1*NBinE]);
   Mbin = 0.02;
   NBinM = 100;
   Mlo = -1;
   Mhi = Mlo + Mbin*NBinM;
   if(make_map) em_map.resize(NBinE*NBinM);
   this->clear();
}

#endif
