#ifndef MEAS_TRACE_HPP
#define MEAS_TRACE_HPP


#include "MPI_Struct.hpp"


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>


/// Average time-dependent trajectory over repeated trials
class Meas_Trace
{

public:

   int icall;
   bool written;
   MPI_Struct mp;
   bool cyclic;

public:

   Meas_Trace(); 

   ~Meas_Trace();

   template<typename Walker>
   void add_sample(const Walker& walker, bool at_equilibrium=true);

   void clear();

   void write();

   void init(long NStepCycle, int NMeas); 

private:

   long NStep;
   long NStepCycle;
   int  NMeas;
   enum { MAXK = 3 };
   std::vector<double> statbuff;
   double* emom;
   double* mmom;
   double* H;
   double* t;
   float NSpin;
   std::ofstream tout;
};


Meas_Trace::Meas_Trace()
{
   cyclic = true;
   mp = MPI_Struct::world();
   clear();
}

Meas_Trace::~Meas_Trace()
{
#  ifdef USE_MPI
   int done;
   MPI_Finalized(&done);
   if( done ) return;
#  endif
   if( !written ) write();
}


void Meas_Trace::clear()
{
   icall = 0;
   written = false;
   NSpin = 0;
   for(int i=0; i<statbuff.size(); i++) statbuff[i] = 0;
   if( tout.is_open() ) tout.close();
   tout.open("Trace.csv");
}

void Meas_Trace::init(long _Cycle, int _NMeas)
{
   NMeas  = _NMeas;
   NStepCycle = _Cycle;
   NStep  = NStepCycle/NMeas;
   statbuff.resize( 2*MAXK*NMeas + 2*NMeas);
   emom = &(statbuff[0*MAXK*NMeas]);
   mmom = &(statbuff[1*MAXK*NMeas]);
   H = &(statbuff[2*MAXK*NMeas+0*NMeas]);
   t = &(statbuff[2*MAXK*NMeas+1*NMeas]);
   this->clear();
}

template<typename Walker>
void Meas_Trace::add_sample(const Walker& walker, bool at_equilibrium)
{
   NSpin = walker.now.V;
   long tloop = static_cast<long>(walker.MCSS);
   if( this->cyclic ) tloop = tloop%NStepCycle;
   long iloop = static_cast<long>(walker.MCSS)/NStepCycle;
   int imeas  = tloop/NStep;
   if( imeas>=NMeas ) return;
   H[imeas] = walker.now.H;
   t[imeas] = tloop;
   double Ek = 1;
   for(int k=0; k<MAXK; k++)
   {
      emom[imeas*MAXK+k] += Ek;
      Ek *= walker.now.E;
   }
   double Mk = 1;
   for(int k=0; k<MAXK; k++)
   {
      mmom[imeas*MAXK+k] += Mk;
      Mk *= walker.now.M;
   }
   if( iloop<5 && walker.iwalk_global==0 )
   {
      tout << walker.MCSS << " " << walker.now.M << " " << walker.now.H << std::endl;
   }
   else
   {
      if( tout.is_open() ) tout.close();
   }
   written = false;
}


void Meas_Trace::write()
{
   // mpi_gather
   int nbuff = statbuff.size() - 2*NMeas;       // Exclude H and t;
   std::vector<double> stat_global(nbuff,0);
#  ifdef USE_MPI
   if(mp.in()) MPI_Allreduce(&(statbuff[0]),&(stat_global[0]),nbuff,MPI_DOUBLE,MPI_SUM,mp.comm);
#  else
   stat_global = stat_buff; 
#  endif
   if( mp.iproc==0 && stat_global[0]>0 )    // stat_global[0] is number of energy samples
   {
      double* emom = &(stat_global[0*MAXK*NMeas]);
      double* mmom = &(stat_global[1*MAXK*NMeas]);
      std::ofstream fout("TraceAve.csv");
      fout.precision(15);
      fout << "# Trace of Magnetization for sweeps and hysteresis" << std::endl;
      fout << "# Column 1: Applied Field" << std::endl;
      fout << "# Column 2: Mean magnetization" << std::endl;
      fout << "# Column 3: Std Dev of magnetization" << std::endl;
      fout << "# Column 4: Normalized mean magnetization M/V" << std::endl;
      fout << "# Column 5: Mean Energy" << std::endl;
      fout << "# Column 6: Std Dev of energy" << std::endl;
      fout << "# Column 7: Normalized mean energy E/V" << std::endl;
      fout << "# Column 8: time in MCSS" << std::endl;
      fout << "# Column 9: counts in bin " << std::endl;
      for(int imeas=0; imeas<NMeas; imeas++)
      {
         double mmean = mmom[MAXK*imeas+1]/mmom[MAXK*imeas+0];
         double mstdd = std::sqrt(mmom[MAXK*imeas+2]/mmom[MAXK*imeas+0]-mmean*mmean);
         double emean = emom[MAXK*imeas+1]/emom[MAXK*imeas+0];
         double estdd = std::sqrt(emom[MAXK*imeas+2]/emom[MAXK*imeas+0]-emean*emean);
         fout << H[imeas] << " " << mmean << " " << mstdd << " " << mmean/NSpin << " "
              << emean << " " << estdd << " " << emean/NSpin << " " << t[imeas] << " " << mmom[imeas*MAXK+0] << std::endl;
      }
   }
   written = true;
}



#endif 
