#ifndef MFIA_MEASURE_HPP
#define MFIA_MEASURE_HPP


#include<string>
#include<vector>
#include<stdio.h>


class Mfia_Measure
{
public:

   int    icall;         // number of calls
   bool   written;       // write since last sample?
   std::string header;   // Canned header information from other objects

public:

   Mfia_Measure();

   ~Mfia_Measure();

   // Add a string to the header information
   void add_header(std::string new_lines);

   // Add a sample state of walker to statistics (called by sampling object)
   template<typename MCWalker>
   void add_sample(MCWalker& walker, bool at_equilibrium=true);

   // Clear acculated data and start over
   void clear();

   // Write current statistics (could be a check-point)
   void write();

   // Enable writing of trajectory
   void trajectory_open(std::string Efname="Mfia-Energy.csv", 
                        std::string Xfname="Mfia-Magnet.csv");

   // Close trajectory and a final check state (calc_observable==change)
   template<typename MCWalker>
   void trajectory_close(MCWalker& final_state);

   // Close trajectory without final state printed
   void trajectory_close();

private:

   FILE* Eout;           // trace file for energy
   FILE* Xout;           // trace file for observables

   template<typename MCWalker>
   void print_trajectory(MCWalker& walker) { print_energy(walker); print_observable(walker); }

   template<typename MCWalker>
   void print_energy(MCWalker& walker);
   void print_energy_header();

   template<typename MCWalker>
   void print_observable(MCWalker& walker);
   void print_observable_header();

private:

   double V;                        // volume of system
   std::string hist_fname;         // finame for histogram
   std::vector<long long> phist;
   std::vector<long long> mhist;

public:

   void histogram_open(int nbin=100, std::string fname="Mfia-Histograms.csv");
   void print_histogram();
 
};



Mfia_Measure::Mfia_Measure() 
{ 
   Eout = NULL;
   Xout = NULL;
   written=true; 
   icall = 0;
   // Any automatic set-up actions go here:
}

Mfia_Measure::~Mfia_Measure() 
{ 
   if( !written ) this->write(); 
   this->trajectory_close();
   // Any clean-up actions go here:
}

void Mfia_Measure::add_header(std::string new_lines)
{
   header = header + new_lines;
}


template<typename MCWalker>
void Mfia_Measure::add_sample(MCWalker& walker, bool at_equilibrium)
{
   // Actions that happen regardless of equilibrium:
   V = walker.now.V;
   icall++;
   if(icall<100000) this->print_trajectory(walker);
   if( !at_equilibrium ) return;
   written = false;
   // Equilibrium Measurements Come Here:
   int nbin = mhist.size();
   if( nbin>0 )
   {
      int ibin = static_cast<int>(static_cast<double>(walker.now.M+V)/static_cast<double>(2*V)*static_cast<double>(nbin));
      if( ibin>=0 && ibin<nbin ) mhist[ibin] += 1;
      ibin = static_cast<int>(static_cast<double>(walker.now.P+V)/static_cast<double>(2*V)*static_cast<double>(nbin));
      if( ibin>=0 && ibin<nbin ) phist[ibin] += 1;
   }
}

void Mfia_Measure::clear()
{
   written = false;
   icall = 0;
   // Actions to clear accumulated statistics go here:
}

void Mfia_Measure::write()
{
   this->print_histogram();
}

void Mfia_Measure::trajectory_open(std::string Efname, std::string Xfname)
{
   if( Eout==NULL && Efname.compare("")!=0 )
   {
      Eout = fopen(Efname.c_str(),"w");
      fprintf(Eout,"%s",header.c_str());
      this->print_energy_header();
   }
   if( Xout==NULL && Xfname.compare("")!=0 )
   {
      Xout = fopen(Xfname.c_str(),"w");
      fprintf(Xout,"%s",header.c_str());
      this->print_observable_header();
   }
}

template<typename MCWalker>
void Mfia_Measure::trajectory_close(MCWalker& walker)
{
   if( Eout!=NULL )
   {
      fprintf(Eout,"# This should match last line of output\n");
      fprintf(Eout,"#");
      this->print_energy(walker);
      fclose(Eout);
      Eout = NULL;
   }
   if( Xout!=NULL )
   {
      fprintf(Xout,"# This should match last line of output\n");
      fprintf(Xout,"#");
      this->print_observable(walker);
      fclose(Xout);
      Xout = NULL;
   }
}


void Mfia_Measure::trajectory_close()
{
   if( Eout!=NULL )
   {
      fclose(Eout);
      Eout = NULL;
   }
   if( Xout!=NULL )
   {
      fclose(Xout);
      Xout = NULL;
   }
}


void Mfia_Measure::print_energy_header()
{
   if( Eout==NULL ) return;
   fprintf(Eout,"# Column 1: Monte Carlo step (imcs)\n");
   fprintf(Eout,"# Column 2: Monte Carlo time (imcs/NGrid)\n");
   fprintf(Eout,"# Column 3: E, the total energy\n");
   fprintf(Eout,"# Column 4: E_N, the exchange energy\n");
   fprintf(Eout,"# Column 5: E_M, the Zeeman energy\n");
   fprintf(Eout,"# Column 6: E_A, the Mean-field energy\n");
}


template<typename MCWalker>
void Mfia_Measure::print_energy(MCWalker& walker)
{
   if( Eout==NULL ) return;
   fprintf(Eout,"%12lld %12lf %18.8le %18.8le %18.8le %18.8le\n", 
                walker.imcs, walker.MCSS, walker.now.E, walker.now.E_N, walker.now.E_H, walker.now.E_A);
}


void Mfia_Measure::print_observable_header()
{
   if( Xout==NULL ) return;
   fprintf(Xout,"# Column 1: Monte Carlo step (imcs)\n");
   fprintf(Xout,"# Column 2: Monte Carlo time (imcs/NGrid)\n");
   fprintf(Xout,"# Column 3: M, the magnetization = sum_i sigma_i\n");
   fprintf(Xout,"# Column 4: P, the staggered magnetization M_A - M_B\n");
   fprintf(Xout,"# Column 5: (N+M)/2, Number of up spins (Column)\n");
   fprintf(Xout,"# Column 6: (D*N-N)/4 = U/2 (Row)\n");
   fprintf(Xout,"# Column 7: U, number of ungerade pairs\n");
   fprintf(Xout,"# Column 8: N, sum over nearest-neighbor spins\n");
}


template<typename MCWalker>
void Mfia_Measure::print_observable(MCWalker& walker)
{
   if( Xout==NULL ) return;
   fprintf(Xout,"%12lld %12lf %12d %12d %12d %12d %12d %12d\n",
           walker.imcs, walker.MCSS, walker.now.M, walker.now.P, walker.now.R, walker.now.C, walker.now.U, walker.now.N);
}


void Mfia_Measure::histogram_open(int nbin, std::string fname)
{
   if(nbin<1) return;
   mhist.resize(nbin);
   phist.resize(nbin);
   for(int ibin=0; ibin<nbin; ibin++) { phist[ibin]=0; mhist[ibin]=0; }
   hist_fname = fname;
}


void Mfia_Measure::print_histogram()
{
   int nbin = mhist.size();
   if( nbin<1 ) return;
   FILE* fout = fopen(hist_fname.c_str(),"w");
   fprintf(fout,"%s",header.c_str());
   fprintf(fout,"# Column 1: Normalized order parameter\n");
   fprintf(fout,"# Column 2: Order parameter\n");
   fprintf(fout,"# Column 3: Normalized M pdf\n");
   fprintf(fout,"# Column 4: Normalized P pdf\n");
   fprintf(fout,"# Column 5: M counts\n");
   fprintf(fout,"# Column 6: P counts\n");
   long long mcts = 0;
   long long pcts = 0;
   for(int ibin=0; ibin<nbin; ibin++)
   {
      mcts += mhist[ibin];
      pcts += phist[ibin];
   }
   double mnorm = static_cast<double>(nbin)/2./static_cast<double>(mcts);
   double pnorm = static_cast<double>(nbin)/2./static_cast<double>(pcts);
   for(int ibin=0; ibin<nbin; ibin++)
   {
      double x = 2.*(static_cast<double>(ibin)+0.5)/static_cast<double>(nbin)-1;
      double xV = V*x;
      double mpdf = mnorm*mhist[ibin];
      double ppdf = pnorm*phist[ibin];
      fprintf(fout,"%12f %12f %12lf %12lf %12lld %12lld\n",x,xV,mpdf,ppdf,mhist[ibin],phist[ibin]);
   }
   fclose(fout);
} 

#endif
