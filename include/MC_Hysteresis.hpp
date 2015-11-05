#ifndef MC_HYSTERESIS_HPP
#define MC_HYSTERESIS_HPP


#include "Random.hpp"
#include "MPI_Struct.hpp"
#include "ProcessTime.hpp"
#include "Null_Measure.hpp"


#include <vector>
#include <cmath>
#include <iostream>


class MC_Hysteresis
{

public:

   MultiplyWithCarry<4> urng;

   float H0;              // The initial field, or upper field value
   float H1;              // The final field, or lower field value
   float kT;              // The temperature of the walkers

   long NStep0;           // Period (in MCSS), or time of first part
   long NStep1;           // time of second parts
   long NCycle;           // Number of cycles to simulate
   long NStep;            // MC Steps per spin (MCSS) between measurements
   MPI_Struct mp;         // Multiprocessing information
   int NWalkPerProcess;   // Number of (embarassingly) parallel      

   ProcessTime ptime;
   double ckpt_freq;
   double ckpt_last;
   double wall_limit;
   bool   at_wall_limit;

public:

   /** The different types of time-dependent field
    * <ul>
    * <li> <i>quench<\i> field held at H0 for NStep0, the swithces to H1 for NStep1. The initial
    * configuration is restored with each cycle.
    * <li> <i>sweep<\i> field is held at H0 for NStep0 to equilbrate, then swept linearly to H1
    * in NStep1 time. The initial configuration is restore with each cycle.
    * <li> <i>sine<\i> A sine-wave field with maximum H0, minimum H1, and a period of NStep1. The
    * is an initial equilibration period of NStep0 for the first cycle, and the initial 
    * configuration is never restored.
    * <li> <i>square<\i> A square field with field at H0 for NStep0, followed by H1 for NStep1.
    * There is an equilibration for NStep0 before the first cycle, but the initial configuration
    * is not restored on subsequent cycles.
    * <li> <i>sawtooth<\i> A triangular field that starts at H0, sweeps to H1 in NStep0, and the
    * sweeps back to H0 in NStep1. Therre is a equilibration of NStep0 before the first cycle,
    * but the initial configuration is not restored on subsequent cycles;
    */
   enum WaveForm { quench, sweep, sine, square, sawtooth } waveform;

public:

   MC_Hysteresis();

   template<typename Walker>
   void InitPool(std::vector<Walker>& walkerpool);

   template<typename Hamiltonian, typename Walker>
   void DoNStep(long NStep, Hamiltonian& hamilton, Walker& walker);

   template<typename Hamiltonian, typename Walker>
   void DoSample(Hamiltonian& hamilton, std::vector<Walker>& walkerpool);

   template<typename Hamiltonian, typename Walker, typename MeasureObj>
   void DoSample(Hamiltonian& hamilton, std::vector<Walker>& walkerpool, MeasureObj& measure);

   float Htime(float t)
   {
      if( t<0 ) return H0;
      long tloop = static_cast<long>(t)%NStepTot;
      switch( waveform )
      {
         case quench:    return (t<NStep0)? H0 : H1;
         case sweep:     return (t<NStep0)? H0 : H0 - sweep0*(t-NStep0); 
         case sine:      return (H0+H1)/2. + (HAmp/2.)*std::cos(twopiT*t);
         case square:    return (tloop<NStep0)? H0 : H1;
         case sawtooth:  return (tloop<NStep0)? H0-sweep0*tloop : H1+sweep1*(tloop-NStep0);
      }
   }

private:

   long NStepTot;
   float sweep0, sweep1;
   float HAmp;
   float twopiT;
};


MC_Hysteresis::MC_Hysteresis()
{
   mp = MPI_Struct::world();         // Multiprocessing information
   NWalkPerProcess = 1;
   H0 = 10;
   H1 =-10;
   kT = 1;
   waveform = sine;
   NStep0 = 0;   
   NStep1 = 0;  
   NCycle = NStep0+NStep1; 
   NStep = 1; 
   ckpt_freq = 10*60;      // 10 minutes
   wall_limit = 24*60*60;  // 1 day
   at_wall_limit = false;
   ckpt_last = 0;
}


template<typename Walker>
void MC_Hysteresis::InitPool(std::vector<Walker>& walkerpool)
{
#  ifdef USE_MPI
   MPI_Comm_rank(mp.comm,&mp.iproc);
   MPI_Comm_size(mp.comm,&mp.nproc);
#  endif
   ptime.mp = mp;
   ptime.start();
   ckpt_last = ptime.elapsed();
   at_wall_limit = false;
   if(this->NWalkPerProcess<1) this->NWalkPerProcess=1;
   walkerpool.resize(NWalkPerProcess);
   for(int iwalk=1; iwalk<walkerpool.size(); iwalk++)
   {
      walkerpool[iwalk] = walkerpool[0];
      walkerpool[iwalk].iwalk_global = NWalkPerProcess*mp.iproc+iwalk;
   }
   NStepTot = NStep0 + NStep1;
   HAmp = H0 - H1;
   twopiT = 8.0*std::atan(1.0)/static_cast<float>(NStepTot);
   sweep0 = HAmp/static_cast<float>(NStep0);
   sweep1 = HAmp/static_cast<float>(NStep1);
}



template<typename Hamiltonian, typename Walker>
void MC_Hysteresis::DoNStep(long NStep, Hamiltonian& hamilton, Walker& walker)
{
   for(long istep=0; istep<NStep; istep++)
   {
      // Change the field
      float Hold = walker.now.H;
      float t = walker.MCSS;
      float Hnew = Htime(t);
      if(Hnew!=Hold) hamilton.change_field(walker.now,Hnew);
      // Propose a Monte Carlo move
      walker.save_initial();                                      // Save state, clear the stack of change information
      hamilton.mc_step(walker,urng);                              // Propose a new state, changes to sigma saved in walker
      // Calculate chance of acceptance
      double delta_E = walker.now.E - walker.old.E;               // Change in energy
      bool accept = (delta_E<0);                                  // always accept changes that decrease energy
      if( !accept )                                               // Use Metropolis criterion to conditionally acecept
      {                                                           // proposals that increase energy
	 double chance = exp(-delta_E/kT);
	 double u = urng();
	 accept = (u<chance);
      }
      // Restore initial state if proposal rejected
      if( !accept ) walker.restore_initial(); 
      walker.imcs++;
      walker.MCSS = static_cast<double>(walker.imcs)/static_cast<double>(walker.now.V);
   }
}


template<typename Hamiltonian, typename Walker>
void MC_Hysteresis::DoSample(Hamiltonian& hamilton, std::vector<Walker>& walkerpool)
{
   // TODO: Replace with a Meas_Trace and Meas_HysteresisLoop depending on waveform
   Null_Measure measure;
   this->DoSample(hamilton,walkerpool,measure);
}


template<typename Hamiltonian, typename Walker, typename MeasureObj>
void MC_Hysteresis::DoSample(Hamiltonian& hamilton, std::vector<Walker>& walkerpool, MeasureObj& measure)
{
   // Get ready
   bool written = false;
   Walker initial_state = walkerpool[0];
   hamilton.change_field(initial_state.now,H0);
   long NSpin = walkerpool[0].now.V;
   // Do equilibration
   long nstep_eq = NStep0*NSpin;
   if( waveform==quench || waveform==sweep ) nstep_eq = 0;
   for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
   {
      walkerpool[iwalk].imcs  = -nstep_eq;
      walkerpool[iwalk].MCSS  = static_cast<double>(walkerpool[iwalk].imcs)/static_cast<double>(NSpin);
      hamilton.change_field(walkerpool[iwalk].now,H0);
      DoNStep(nstep_eq,hamilton,walkerpool[iwalk]);
      measure.add_sample(walkerpool[iwalk]);
   }
   // Do the simulation
   long nmeas_cycle = NStepTot/NStep;                            // number of measurements per cycle
   long nstep_meas  = NStep*NSpin;                               // number of mcs to between measurements
   long nstep_cycle = NStepTot*NSpin;                            // number of mcs for whole cycle
   long nstep_endup = nstep_cycle - nstep_meas*nmeas_cycle;      // potentially some extra steps to finish cycle
   if( nstep_endup>0 ) nmeas_cycle++;
   for(long icycle=0; icycle<NCycle && !at_wall_limit; icycle++)
   {
      if( waveform==quench || waveform==sweep )              // These waveforms are independent repeated trials
      {
         for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
         {
            walkerpool[iwalk] = initial_state;
            walkerpool[iwalk].imcs = 0;
            walkerpool[iwalk].MCSS = 0;
            measure.add_sample(walkerpool[iwalk]);
         }
      }
      for(long imeas=0; imeas<nmeas_cycle && !at_wall_limit; imeas++)
      {
         long nstep_this = nstep_meas;
         if( (imeas+1)==nmeas_cycle && nstep_endup>0 ) nstep_this = nstep_endup;
         for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
            DoNStep(nstep_this,hamilton,walkerpool[iwalk]);
         for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
            measure.add_sample(walkerpool[iwalk]);
         written = false;
         double ptime_sec = ptime.elapsed();
         double ckpt_sec = ptime_sec - ckpt_last;
         if( ckpt_sec>ckpt_freq )
         {
            ckpt_last = ptime_sec;
            measure.write();
            written = true;
         }
         double remain_sec = wall_limit - ptime.elapsed();
         at_wall_limit = (remain_sec<ckpt_freq);
      }
   }
   if( !written )
   {
      measure.write();
      written = true;
   }
}




#endif
