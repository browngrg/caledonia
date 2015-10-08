#ifndef MC_METROPOLIS_HPP
#define MC_METROPOLIS_HPP


#include "Random.hpp"
#include "MPI_Struct.hpp"
#include "MPITypeTraits.hpp"
#include "ProcessTime.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>


class MC_Metropolis
{

public:

   MultiplyWithCarry<4> urng;  // The random number generator

   float kTmin,kTmax;          // The temperature range
   std::vector<float> kTlist;  // Number of different temperature

   long NStep;                 // Number of steps (MCSS) to take between measurements to minimize correlations
   long NMeas;                 // Number of measurements to take
   long NTherm;                // Number of measurement blocks before begin averaging
   long NExchg;                // Number of mesurements blocks between replica exchange attempts
   long re_iter;               // Loop control variable for replica exchange

   MPI_Struct mp;              // Multiprocessor information

   int NWindow;
   int NWalkPerProcess;

//   std::string filename;       // For output of energy statistics

   ProcessTime ptime;
   double ckpt_freq;              // How often to write results (in seconds)
   double ckpt_last;              // When the last ckpt write occurred
   double wall_limit;             // Longest time to run (in seconds)
   bool   at_wall_limit;

public:

   MC_Metropolis()
   {
      this->mp = MPI_Struct::world();
      this->kTmin  = 0.5;
      this->kTmax  = 3;
      this->NStep  = 10;
      this->NMeas  = 10000;
      this->NTherm = 100;
      this->NExchg = 10;
      this->NWindow = 1;
      this->NWalkPerProcess = 10;
      this->re_iter = 0;
      at_wall_limit = false;
      wall_limit = 24*60*60;
      ckpt_freq  = 10*60;
      ckpt_last  = 0;
   }

   ~MC_Metropolis()
   {
   }

   template<typename Hamiltonian, typename MCWalker>
   void DoNStep(long NStep, Hamiltonian& hamilton, MCWalker& walker);

   template<typename Hamiltonian, typename MCWalker, typename MeasureObj>
   void DoSample(Hamiltonian& model, std::vector<MCWalker>& walker, MeasureObj& measure_obj);

   template<typename MCWalker>
   void InitPool(std::vector<MCWalker>& walkerpool);

   template<typename Model, typename Walker>
   void DoReplicaExchange(Model& model, std::vector<Walker>& walkerpool);

   template<typename MCWalker>
   void write(const std::vector<MCWalker>& walkerpool);

   template<typename Model, typename MCWalker>
   void read_config(Model& model, std::vector<MCWalker>& walkerpool);

   template<typename MCWalker> 
   void set_kT(const std::vector<float>& new_kT, std::vector<MCWalker>& walkerpool);

   template<typename MCWalker>
   void CalcOptimal_kT(std::vector<MCWalker>& walker, bool DoAdjust=false);

   std::string header();

private:

   template<typename MCWalker>
   bool check_walker(std::string file, int line, const MCWalker& walker);

};


/// Basic Metropolis sampling engine.
/** This method advances a Monte Carlo walker N Monte Carlo steps (mcs), accepting or
 *  rejecting proposed moves based on the energy of the walker.
 *  @param NStep is the number of simple Monte Carlo steps to attempt.
 *  @param hamilton calculates the energy and other macrovariables while it proposes an attempted move
 *  @param walker is a single Monte Carlo walker with one microconfiguration sigma, and records of
 *         associated macrovariables, <em>now</em> and <em>old</em>. 
 */
template<typename Hamiltonian,typename MCWalker>
void MC_Metropolis::DoNStep(long NStep, Hamiltonian& hamilton, MCWalker& walker)
{
   for(long istep=0; istep<NStep; istep++)
   {
      // Propose a Monte Carlo move
      walker.save_initial();                                      // Save state, clear the stack of change information
      hamilton.mc_step(walker,urng);                              // Propose a new state, changes to sigma saved in walker
      // Calculate chance of acceptance
      double delta_E = walker.now.E - walker.old.E;               // Change in energy
      bool accept = (delta_E<0);                                  // always accept changes that decrease energy
      if( !accept )                                               // Use Metropolis criterion to conditionally acecept
      {                                                           // proposals that increase energy
	 double chance = exp(-delta_E/walker.kT);
	 double u = urng();
	 accept = (u<chance);
      }
      // Restore initial state if proposal rejected
      if( !accept ) walker.restore_initial(); 
      walker.imcs++;
   }
   walker.MCSS = static_cast<double>(walker.imcs)/static_cast<double>(walker.now.V);
}


template<typename Hamiltonian, typename MCWalker, typename MeasureObj>
void MC_Metropolis::DoSample(Hamiltonian& hamilton, std::vector<MCWalker>& walker, MeasureObj& measure_obj)
{
   int NSpin = walker[0].now.V;
   // Do sampling
   for(int iwalk=0; iwalk<walker.size(); iwalk++)
   {
      walker[iwalk].imcs = -NTherm*NStep*NSpin;
      walker[iwalk].MCSS = static_cast<double>(walker[iwalk].imcs)/static_cast<double>(walker[iwalk].now.V);
   }
   long NMeasTot = NMeas + NTherm;
   bool written = true;
   for(long imeas=0; imeas<NMeasTot && !at_wall_limit; imeas++)
   {
      for(int iwalk=0; iwalk<walker.size(); iwalk++)
         this->DoNStep(NSpin*NStep,hamilton,walker[iwalk]);
      bool at_equilibrium = (imeas>=NTherm);
      for(int iwalk=0; iwalk<walker.size(); iwalk++)
      {
         measure_obj.add_sample(walker[iwalk],at_equilibrium);
         if( at_equilibrium ) walker[iwalk].sample();
      }
      written = false;
      double ptime_sec = ptime.elapsed();
      double ckpt_sec = ptime_sec - ckpt_last;
      if( ckpt_sec>ckpt_freq )
      {
         ckpt_last = ptime_sec;
         measure_obj.write();
         this->write(walker);
         written = true;
         this->CalcOptimal_kT(walker,false);
      }
      if( (imeas%NExchg)==0 ) DoReplicaExchange(hamilton,walker);
      double remain_sec = wall_limit - ptime.elapsed();
      at_wall_limit = (remain_sec<ckpt_freq);
   }
   if( !written )
   {
      measure_obj.write();
      this->write(walker);
      this->CalcOptimal_kT(walker,false);
   }
}

template<typename MCWalker>
void MC_Metropolis::InitPool(std::vector<MCWalker>& walkerpool)
{
#  ifdef USE_MPI
   if( mp.in() ) MPI_Comm_rank(mp.comm,&mp.iproc);
   if( mp.in() ) MPI_Comm_size(mp.comm,&mp.nproc);
#  endif
   ptime.mp = mp;
   ptime.start();
   ckpt_last = ptime.elapsed();
   at_wall_limit = false;
   if(this->NWalkPerProcess<1) this->NWalkPerProcess=1;
   walkerpool.resize(this->NWalkPerProcess);
   for(int iwalk=1; iwalk<walkerpool.size(); iwalk++)
      walkerpool[iwalk]=walkerpool[0];
   for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
   {
      walkerpool[iwalk].iwalk_global = NWalkPerProcess*mp.iproc+iwalk;
   }
   // Assign temperatures in a zig-zag pattern over range
   std::vector<float> new_kTlist;
   std::ifstream fin("MC_MetropTemps.csv");
   if( fin && fin.is_open() )
   {
      // Read temps
      while( fin && fin.is_open() && !fin.eof() )
      {
         std::string buff;
         std::getline(fin,buff);
         if( buff.size()>0 && buff[0]!='#' )
         {
            float kT;
            sscanf(buff.c_str(),"%f",&kT);
            new_kTlist.push_back(kT);
         }
      }
      if( mp.iproc==0 ) std::cout << "Read " << new_kTlist.size() << " temperatures from input" << std::endl;
   }
   if( new_kTlist.size()==0 )
   {
      // Default is evenly spaced temperatures
      int kTnum = NWalkPerProcess*mp.nproc;
      float dTemp = 0;
      if(kTnum>1) { dTemp = (kTmax-kTmin)/static_cast<float>(kTnum-1); }
      new_kTlist.resize(kTnum);
      for(int i=0; i<kTnum; i++) new_kTlist[i] = kTmin + dTemp*i;
   }
   this->set_kT(new_kTlist,walkerpool);
   if(false)
   {
      std::ostringstream iproc_str; iproc_str << mp.iproc;
      std::string filename = "MC_MetropTemps_" + iproc_str.str() + ".csv";
      std::ofstream fout(filename.c_str());
      for(int i=0; i<walkerpool.size(); i++)
         fout << i << " " << walkerpool[i].iwalk_global << " " << walkerpool[i].kTi << " " << walkerpool[i].kT << std::endl;
   }
}


template<typename Walker>
void MC_Metropolis::set_kT(const std::vector<float>& new_kT, std::vector<Walker>& walkerpool)
{
   this->kTlist = new_kT;
   int kTnum = kTlist.size();
   int kTmax = NWalkPerProcess*mp.nproc;                 // limited by number of walkers
   if( kTnum>kTmax ) { kTlist.resize(kTmax); kTnum=kTlist.size(); }
   int ii = (NWalkPerProcess*mp.iproc)%(2*kTnum);
   bool up = true;
   if( ii>=kTnum ) 
   {
      ii = (kTnum-1) - (ii-kTnum);
      up =  false;
   }
   for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
   {
      walkerpool[iwalk].kTi = ii;
      walkerpool[iwalk].kT = kTlist[ walkerpool[iwalk].kTi ];
      for(int i=0; i<Walker::MAXK; i++) walkerpool[iwalk].emom[i]=0;
      if(up) ii++; else ii--;
      if(ii>=kTnum) { ii=(kTnum-1)-(ii-kTnum); up=false; }
      if(ii<0)      { ii=-ii; up=true; }
   }
   kTmin = kTlist.front();
   kTmax = kTlist.back();
   if( false && mp.in() )
   {
      char filename[100];
      sprintf(filename,"MC_MetropSetkT_%d",mp.iproc);
      std::ofstream fout(filename);
      fout << "# Column 1: iwalk local" << std::endl;
      fout << "# Column 2: iwalk global" << std::endl;
      fout << "# Column 3: kTi" << std::endl;
      fout << "# Column 4: kT" << std::endl;
      fout << "# Column 5: kT[kTi]" << std::endl;
      for(int iwalk=0; iwalk<walkerpool.size(); iwalk++)
      {
         fout << iwalk << " " << walkerpool[iwalk].iwalk_global << " " << walkerpool[iwalk].kTi << " " << walkerpool[iwalk].kT << " " << kTlist[ walkerpool[iwalk].kTi ] << std::endl;
      }
   }
}


std::string MC_Metropolis::header()
{
   std::string result;
   char buffer[1024];
   sprintf(buffer,"# Metropolis Sampling kTmin=%f kTmax=%f\n",kTmin,kTmax);
   result += std::string(buffer);
   sprintf(buffer,"# NMeas=%ld NTherm=%ld NStep=%ld\n",NMeas,NTherm,NStep);
   result += std::string(buffer);
   return result;
}


template<typename MCWalker>
void MC_Metropolis::CalcOptimal_kT(std::vector<MCWalker>& walker, bool DoAdjust)
{
   const int MAXK = MCWalker::MAXK;
   // mpi_gather
   int kTnum = kTlist.size();
   int buff_size = kTnum*MAXK;
   std::vector<double> emom_buff(buff_size,0);
   for(int iwalk=0; iwalk<walker.size(); iwalk++)
   {
      int kTi = walker[iwalk].kTi*MAXK;
      for(int k=0; k<MAXK; k++)
         emom_buff[kTi+k] += walker[iwalk].emom[k];
   }
   std::vector<double> emom_global(buff_size);
   MPI_Allreduce(&(emom_buff[0]),&(emom_global[0]),buff_size,MPI_DOUBLE,MPI_SUM,mp.comm);
   // Adjust temperatures
   bool bad = false;
   std::vector<float> mean(kTnum);
   std::vector<float> sigma(kTnum);
   for(int kTi=0; kTi<kTnum; kTi++)
   {
      bad = bad || (emom_global[kTi*MAXK+0]==0);
      mean[kTi] = emom_global[kTi*MAXK+1]/emom_global[kTi*MAXK+0];
      float var  = emom_global[kTi*MAXK+2]/emom_global[kTi*MAXK+0] - mean[kTi]*mean[kTi];
      if( var<0 ) var=0;
      sigma[kTi] = std::sqrt(var);
   }
   if( bad ) return;
   std::vector<float> new_kT;
   new_kT.push_back(kTmin);
   float Ei = mean[0];
   float sigi = sigma[0]; 
   int ilast = 0;
   while( new_kT.back()<kTmax && new_kT.size()<100000 )
   {
      int i = ilast+1;
      if( sigi==0 && i<(sigma.size()-1) )
      {
         while( sigma[i]==0 && i<(sigma.size()-1) ) i++;
         ilast = i;
         Ei = mean[i];
         sigi = sigma[i];
         new_kT.push_back(kTlist[i]);
         continue;
      }
      const float f = 1.5;        // overlap factor
      float Ej = Ei + f*sigi;
      float sigj = sigi;
      float x = 0;
      int jj = ilast;
      for(int iter=0; iter<10; iter++)
      {
         jj = ilast;
         while( jj<mean.size() && mean[jj]<Ej ) jj++;
         if(jj>=mean.size()) jj = mean.size()-1;
         if(jj==0) jj=1;
         x = ( Ej - mean[jj-1] ) / ( mean[jj] - mean[jj-1] );
         sigj = sigma[jj-1] + x*(sigma[jj]-sigma[jj-1]);
         Ej = Ei + f*(sigi+sigj)/2.;
      }
      ilast = jj;
      if( jj>=(mean.size()-1) )
      {
         new_kT.push_back(kTmax);
      }
      else
      {   
         new_kT.push_back( std::min( kTmax,  kTlist[jj-1] + x*(kTlist[jj]-kTlist[jj-1]) ) );
         Ei = Ej;
         sigi = sigj;
      }
   } 
   if( mp.iproc==0 )
   {
      std::ofstream fout("MC_MetropOptimal.csv");
      fout << "# Optimal temperatures for MC_Metropolis" << std::endl;
      fout << "# Column 2: kT" << std::endl;
      for(int i=0; i<new_kT.size(); i++) fout << new_kT[i] << std::endl;
   }
   if( DoAdjust ) { this->set_kT(new_kT,walker); }
}


/// Read configurations from previous run at same conditions.
/** This method is not very smart about the configurations it reads. It ignores the temperature
 *  and just copies the configuration and calculates all of the observables from scratch. It will
 *  not change the number of walkers, resize the microconfiguration <em>sigma<em> or anything.
 *  Its there to load the configurations at the end of the last run to allow daisy chaining
 *  runs of systems that are hard to relax. 
 */
template<typename Hamilton, typename MCWalker>
void MC_Metropolis::read_config(Hamilton& hamilton, std::vector<MCWalker>& walker)
{
   char cfilename[100];
   sprintf(cfilename,"MC_MetropConfig_%d",mp.iproc);
   std::ifstream fin(cfilename);
   if( fin && fin.is_open() )
   {
      for(int iwalk=0; iwalk<walker.size(); iwalk++)
      {
         double kT,E,imcs;
         fin >> kT >> E >> imcs;
         if( kT !=walker[iwalk].kT )
            std::cout << __FILE__ << ":" << __LINE__ << " Config temp not match set temp" << std::endl;
         for(int i=0; i<walker[iwalk].sigma.size(); i++) fin >> walker[iwalk].sigma[i];
         hamilton.calc_observable(walker[iwalk].sigma,walker[iwalk].now);
      }
   }
}


template<typename MCWalker>
void MC_Metropolis::write(const std::vector<MCWalker>& walker)
{
   // write configs
   char cfilename[100];
   sprintf(cfilename,"MC_MetropConfig_%d",mp.iproc);
   if( mp.iproc<100 )
   {
      std::ofstream fout(cfilename);
      if( fout && fout.is_open() )
      {
         for(int iwalk=0; iwalk<walker.size(); iwalk++)
         {
            fout << walker[iwalk].kT << " " << walker[iwalk].now.E << " " << walker[iwalk].imcs;
            for(int i=0; i<walker[iwalk].sigma.size(); i++) fout << " " << walker[iwalk].sigma[i];
            fout << std::endl;
         }
      }
   }
   // Output Energy moments
   const int MAXK = MCWalker::MAXK;
   int kTnum = kTlist.size();
   int buff_size = kTnum*MAXK;
   std::vector<double> emom_buff(buff_size,0);
   for(int iwalk=0; iwalk<walker.size(); iwalk++)
   {
      int kTi = walker[iwalk].kTi*MAXK;
      for(int k=0; k<MAXK; k++)
         emom_buff[kTi+k] += walker[iwalk].emom[k];
   }
   std::vector<double> emom_global(buff_size);
   MPI_Allreduce(&(emom_buff[0]),&(emom_global[0]),buff_size,MPI_DOUBLE,MPI_SUM,mp.comm);
   std::vector<long> re_buff(kTnum,0);
   for(int iwalk=0; iwalk<walker.size(); iwalk++)
   {
      re_buff[ walker[iwalk].kTi ] += walker[iwalk].re_swap;
   }
   std::vector<long> re_swap(kTnum);
   MPI_Allreduce(&(re_buff[0]),&(re_swap[0]),kTnum,MPI_LONG,MPI_SUM,mp.comm);
   // Root process writes
   if( mp.iproc==0 )
   {
      std::ofstream fout("MC_Metrop.csv");
      fout << "# Basic Energy statistics for walkers" << std::endl;
      fout << header();
      fout << "# Column 1: kT" << std::endl;
      fout << "# Column 2: <E>" << std::endl;
      fout << "# Column 3: <E^2>-<E>^2" << std::endl;
      fout << "# Column 4: standard deviation" << std::endl;
      fout << "# Column 5: overlap of T_i + T_i+1 = mean(sigma)/delta(E)" << std::endl;
      fout << "# Column 6: number of replica-exchange swaps" << std::endl;
      fout << "# Other columns raw accumulated moments" << std::endl;
      float mean_old;
      float stdd_old;
      for(int kTi=0; kTi<kTnum; kTi++)
      {
         float mean = emom_global[kTi*MAXK+1]/emom_global[kTi*MAXK+0];
         float var  = emom_global[kTi*MAXK+2]/emom_global[kTi*MAXK+0] - mean*mean;
         float stdd = std::sqrt(var);
         float overlap = 0;
         if( kTi>0 )
            overlap = 1.-(stdd+stdd_old)/2./(mean-mean_old);
         fout << std::setw(5) << kTlist[kTi] << " " << mean << " " << var << " " << stdd << " " 
              << overlap << " " << re_swap[kTi];
         for(int i=0; i<MAXK; i++) fout << " " << std::setprecision(15) << emom_global[kTi*MAXK+i];
         fout << std::endl;
         mean_old = mean;
         stdd_old = stdd;
      }
   }
}

// This can be used for debugging
template<typename MCWalker>
bool MC_Metropolis::check_walker(std::string file, int line, const MCWalker& walker)
{
   // this is for ising-like walkers
   int dm = walker.now.M - walker.old.M;
   if( (dm*dm)>4 )
   {
      std::cerr << file << ":" << line << " step flips more than one spin" << std::endl
                << "M: " << walker.old.M << " -> " << walker.now.M << " delta=" << dm << std::endl
                << "X: " << walker.old.X << " -> " << walker.now.X << std::endl
                << "E: " << walker.old.E << " -> " << walker.now.E << std::endl
                << "sigma = " << walker.sigma[0] << " " << walker.sigma[1] << std::endl
                << "imcs=" << walker.imcs << " MCSS=" << walker.MCSS << std::endl;
      return false;
   }
   return true;
}

template<typename Model, typename Walker>
void MC_Metropolis::DoReplicaExchange(Model& model, std::vector<Walker>& walkerpool)
{
   bool send_right = re_iter % 2;
#  ifdef USE_MPI
   if( !mp.in() ) return;
   const int iproc = mp.iproc;
   const int nproc = mp.nproc;
   if( nproc>1 )
   {
      const int irght = (iproc+1)%nproc;
      const int ileft = (iproc-1+nproc)%nproc;
      bool when_right = (iproc+1) % 2;
      MPI_Status status;
      const int jwalk = walkerpool.size()-1;
      // Send E_right, Receive E_left
      int tag = re_iter % 1000;
      std::vector<double> right(2),left(2);
      right[0] = walkerpool[jwalk].now.E; 
      right[1] = walkerpool[jwalk].kT; 
      MPI_Send(&(right[0]),2,MPI_DOUBLE,irght,tag,mp.comm);
      MPI_Recv(&(left[0]), 2,MPI_DOUBLE,ileft,tag,mp.comm,&status);   
      double E_left = left[0];
      double kT_left = left[1];
      // Decide accept/reject
      double delta = (1./kT_left - 1./walkerpool[0].kT) * (E_left - walkerpool[0].now.E);
      bool accept_left = (delta>0);                                  // always accept changes that decrease energy
      if( !accept_left )                            
      {
         if( delta==0 )
            accept_left = urng()<0.5;
         else
            accept_left = (urng()<exp(delta));
      }
      bool accept_right;
      MPI_Send(&accept_left, 1,MPI_C_BOOL,ileft,tag+1000,mp.comm);
      MPI_Recv(&accept_right,1,MPI_C_BOOL,irght,tag+1000,mp.comm,&status);
      // Swap configurations (pattern borrowed from MC_WangLandau.hpp
      int NSPIN = walkerpool[0].sigma.size();
      MPI_Datatype MPIConfigType = MPITypeTraits<typename Walker::ConfigType::value_type>::mpitype;
      std::vector<typename Walker::ConfigType::value_type> buffer(NSPIN);
      if( send_right==when_right )
      {
         if( accept_right )
         {
            MPI_Send(&(walkerpool[jwalk].sigma[0]),NSPIN,MPIConfigType,irght,tag+3000,mp.comm);
            MPI_Recv(&(buffer[0]),NSPIN,MPIConfigType,irght,tag+2000,mp.comm,&status);
            walkerpool[jwalk].re_swap++;
            walkerpool[jwalk].sigma = buffer;
            model.calc_observable(walkerpool[jwalk].sigma,walkerpool[jwalk].now);
         }
         if( accept_left )
         {
            MPI_Recv(&(buffer[0]),NSPIN,MPIConfigType,ileft,tag+3000,mp.comm,&status);
            MPI_Send(&(walkerpool[0].sigma[0]),NSPIN,MPIConfigType,ileft,tag+2000,mp.comm);
            walkerpool[0].re_swap++;
            walkerpool[0].sigma = buffer;
            model.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
         }
      }
      else
      {
         if( accept_left )
         {
            MPI_Recv(&(buffer[0]),NSPIN,MPIConfigType,ileft,tag+3000,mp.comm,&status);
            MPI_Send(&(walkerpool[0].sigma[0]),NSPIN,MPIConfigType,ileft,tag+2000,mp.comm);
            walkerpool[0].re_swap++;
            walkerpool[0].sigma = buffer;
            model.calc_observable(walkerpool[0].sigma,walkerpool[0].now);
         }
         if( accept_right )
         {
            MPI_Send(&(walkerpool[jwalk].sigma[0]),NSPIN,MPIConfigType,irght,tag+3000,mp.comm);
            MPI_Recv(&(buffer[0]),NSPIN,MPIConfigType,irght,tag+2000,mp.comm,&status);
            walkerpool[jwalk].re_swap++;
            walkerpool[jwalk].sigma = buffer;
            model.calc_observable(walkerpool[jwalk].sigma,walkerpool[jwalk].now);
         }
      }
   }
#  endif
   for(int iwalk=(send_right)? 0:1; iwalk<(walkerpool.size()-1); iwalk+=2)
   {
      double betai = 1./walkerpool[iwalk].kT;
      double betaj = 1./walkerpool[iwalk+1].kT;
      double Ei = walkerpool[iwalk].now.E;
      double Ej = walkerpool[iwalk+1].now.E;
      double delta = (betai-betaj)*(Ei-Ej);
      bool accept = (delta>0);                                  // always accept changes that decrease energy
      if( !accept )                            
      {
         if( delta==0 )
            accept = urng()<0.5;
         else
            accept = (urng()<exp(delta));
      }
      if( accept )
      {
         typename Walker::ConfigType sigma_temp = walkerpool[iwalk].sigma;
         walkerpool[iwalk].sigma = walkerpool[iwalk+1].sigma;
         walkerpool[iwalk+1].sigma = sigma_temp;
         typename Walker::ObserveType now_temp = walkerpool[iwalk].now;
         walkerpool[iwalk].now = walkerpool[iwalk+1].now;
         walkerpool[iwalk+1].now = now_temp;
         walkerpool[iwalk].re_swap++;
         walkerpool[iwalk+1].re_swap++;
      }
   }
   re_iter++;
}


#endif 
