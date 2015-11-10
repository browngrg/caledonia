#ifndef MEAS_DIFFUSION_HPP
#define MEAS_DIFFUSION_HPP


#include<stdio.h>
#include<vector>
#include<map>


class Meas_Diffusion
{

public:

   int icall;
   bool written;

   int delay;

public:

   Meas_Diffusion(float Emin=0, float Emax=1, float _Ebin=1);

   void init(float Emin, float Emax, float Ebin);

   template<typename SIMULATION> void init(SIMULATION& sim);

   void clear() { old.clear(); for(int i=0; i<NBinE; i++) Emom[i]=0; }

   void write();

   template<typename MCWalker>
   void add_sample(const MCWalker& walker, bool at_equilibrium=true);

private:

   float Elo, Ehi, Ebin;
   int NBinE;

   int bin(float E) const { int ibin=static_cast<int>((E-Elo)/Ebin); if(ibin<0) ibin=0; if(ibin>=NBinE) ibin=NBinE-1; return ibin; }

   enum { MAXK = 3 };
   std::vector<double> Emom;
   std::map<int,double> old;
};


Meas_Diffusion::Meas_Diffusion(float Emin, float Emax, float _Ebin)
{ 
   icall = 0;
   written = true;
   delay=1000; 
   init(Emin,Emax,_Ebin); 
}


void Meas_Diffusion::init(float Emin, float Emax, float _Ebin)
{
   float Erange = Emax - Emin;
   Ebin = _Ebin;
   NBinE = Erange/Ebin;
   Elo = Emin;
   Ehi = Elo + Ebin*NBinE;
   Emom.resize( MAXK*NBinE );
   clear();
}


template<typename SIMULATION> void Meas_Diffusion::init(SIMULATION& sim)
{
   this->init(sim.Elo,sim.Ehi,sim.Ehi);
}


template<typename MCWalker>
void Meas_Diffusion::add_sample(const MCWalker& walker, bool at_equilibrium)
{
   if( !at_equilibrium ) return;
   icall++;
   long long imcs = walker.imcs;
   if( (imcs%delay) != 0 ) return;
   int iwalk = walker.iwalk_global;
   double Enow = walker.now.E;
   std::map<int,double>::iterator iEold = old.find(iwalk);
   if( iEold != old.end() )
   {
      written = false;
      double Eold = iEold->second;
      double dE = Enow - Eold;
      int ibin = bin(Eold);
      double Ek = 1;
      for(int k=0; k<MAXK; k++)
      {
         Emom[ ibin*MAXK + k ] += Ek;
         Ek *= dE;
      }
   }
   old[iwalk] = Enow;
}


void Meas_Diffusion::write()
{
   std::vector<double> Emom_global(Emom.size());
#  ifdef USE_MPI
   for(int i=0; i<Emom_global.size(); i++) Emom_global[i] = 0;
   MPI_Allreduce(&(Emom[0]),&(Emom_global[0]),Emom.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   int iproc = 0;
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   if( iproc!=0 ) return;
#  else
   Emom_global = Emom;
#  endif
   FILE* fout = fopen("Meas_Diffusion.csv","w");
   fprintf(fout,"# Diffusion of Walkers vs Energy for delay of %d mcs\n",delay);
   fprintf(fout,"# Column 1: Energy\n");
   fprintf(fout,"# Column 2: Mean bias <E1-E0>\n");
   fprintf(fout,"# Column 3: Mean-squared change <(E1-E0)^2>\n");
   fprintf(fout,"# Column 4: RMS change sqrt(Col 3)\n");
   fprintf(fout,"# Column 5: number of cts\n");
   fprintf(fout,"# Column 6: accumulated sum deltaE\n");
   fprintf(fout,"# Column 7: accumulated sum deltaE^2\n");
   for(int ibin=0; ibin<NBinE; ibin++)
   {
      double N = Emom_global[ibin*MAXK+0];
      if( N>1 ) 
      {
         double E = Elo + (ibin+0.5)*Ebin;
         double mean = Emom_global[ibin*MAXK+1]/N;
         double meansq = Emom_global[ibin*MAXK+2]/N;
         double rms = std::sqrt(meansq);
         fprintf(fout,"%lf %lf %lf %lf",E,mean,meansq,rms);
         for(int k=0; k<MAXK; k++) fprintf(fout," %20.16le",Emom_global[ibin*MAXK+k]);
         fprintf(fout,"\n");
      }
   }
}
#endif 
