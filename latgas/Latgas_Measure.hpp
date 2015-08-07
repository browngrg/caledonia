#ifndef LATGAS_MEASURE_HPP
#define LATGAS_MEASURE_HPP


#include<string>
#include<stdio.h>


class Latgas_Measure
{
public:

   int    icall;         // number of calls
   bool   written;       // write since last sample?
   std::string header;   // Canned header information from other objects

public:

   Latgas_Measure();

   ~Latgas_Measure();

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
   void trajectory_open(std::string Efname="Latgas-Energy.csv", 
                        std::string Xfname="Latgas-Coverage.csv");

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

};



Latgas_Measure::Latgas_Measure() 
{ 
   Eout = NULL;
   Xout = NULL;
   written=true; 
   icall = 0;
   // Any automatic set-up actions go here:
}

Latgas_Measure::~Latgas_Measure() 
{ 
   if( !written ) this->write(); 
   this->trajectory_close();
   // Any clean-up actions go here:
}

void Latgas_Measure::add_header(std::string new_lines)
{
   header = header + new_lines;
}


template<typename MCWalker>
void Latgas_Measure::add_sample(MCWalker& walker, bool at_equilibrium)
{
   // Actions that happen regardless of equilibrium:
   icall++;
   this->print_trajectory(walker);
   if( !at_equilibrium ) return;
   written = false;
   // Equilibrium Measurements Come Here:
}

void Latgas_Measure::clear()
{
   written = false;
   icall = 0;
   // Actions to clear accumulated statistics go here:
}

void Latgas_Measure::write()
{
   printf("# Latgas_Measure called %d times\n",icall);
   written = true;
}

void Latgas_Measure::trajectory_open(std::string Efname, std::string Xfname)
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
void Latgas_Measure::trajectory_close(MCWalker& walker)
{
   if( Eout!=NULL )
   {
      fprintf(Eout,"# This should match last line of output\n");
      fprintf(Eout,"#");
      this->print_energy(walker);
      fclose(Eout);
      Eout == NULL;
   }
   if( Xout!=NULL )
   {
      fprintf(Xout,"# This should match last line of output\n");
      fprintf(Xout,"#");
      this->print_observable(walker);
      fclose(Xout);
      Xout == NULL;
   }
}


void Latgas_Measure::trajectory_close()
{
   if( Eout!=NULL )
   {
      fclose(Eout);
      Eout == NULL;
   }
   if( Xout!=NULL )
   {
      fclose(Xout);
      Xout == NULL;
   }
}


void Latgas_Measure::print_energy_header()
{
   if( Eout==NULL ) return;
   fprintf(Eout,"# Column 1: Monte Carlo step (imcs)\n");
   fprintf(Eout,"# Column 2: Monte Carlo time (imcs/NGrid)\n");
   fprintf(Eout,"# Column 3: E, the total energy\n");
   fprintf(Eout,"# Column 4: Zero-body interaction\n");
   fprintf(Eout,"# Column 5: One-body interactions\n");
   fprintf(Eout,"# Column 6: Two-body interactions\n");
   fprintf(Eout,"# Column 7: Three-body interactions\n");
   fprintf(Eout,"# Column 8: Four-body interactions\n");
}

template<typename MCWalker>
void Latgas_Measure::print_energy(MCWalker& walker)
{
   if( Eout==NULL ) return;
   fprintf(Eout,"%12lld %12lf %18.8le %18.8le %18.8le %18.8le %18.8le %18.8le\n",
           walker.imcs, walker.MCSS, walker.now.E, walker.now.ebody[0], walker.now.ebody[1], walker.now.ebody[2], walker.now.ebody[3], walker.now.ebody[4]);
}


void Latgas_Measure::print_observable_header()
{
   if( Xout==NULL ) return;
   fprintf(Xout,"# Column 1: Monte Carlo step (imcs)\n");
   fprintf(Xout,"# Column 2: Monte Carlo time (imcs/NGrid)\n");
   fprintf(Xout,"# Column 3: Theta[0]\n");
   fprintf(Xout,"# Column 4: Theta[1]\n");
}

template<typename MCWalker>
void Latgas_Measure::print_observable(MCWalker& walker)
{
   if( Eout==NULL ) return;
   fprintf(Xout,"%12lld %12lf %18.8le %18.8le\n",
           walker.imcs, walker.MCSS, walker.now.theta[0], walker.now.theta[1]);
}

#endif
