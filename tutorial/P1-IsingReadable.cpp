// P1_IsingReadable.cpp -- easily readable Ising Monte Carlo simulation
// Author: Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// Date: June 27, 2016
//
// See "The Hobbyhorse of Magnetic Systems: The Ising Model," 
// E. Ibarra-Garcia-Padilla and F.J. Poveda-Cuevas, arXiv:106.05800
// For a discussion of the Iisng model, and a source for the notation here


#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<ctime>


// This calculates the sum of the values of s[j] for the nearest neighbor spins of k.
// Since this gets calculated in two different places, it makes sense to make a function.
int calc_sumj(int k, const std::vector<short>& s, int L)
{
   int ky = k/L; 
   int kx = k%L;         
   int Lm1 = L-1;
   int sum_j = 0;
   int jx = (kx + Lm1) % L;   // jx = kx-1 = (kx+L-1)%L with periodic boundary conditions, and no negative numbers
   sum_j += s[jx+ky*L];
   jx = (kx+1)%L;             // jx = kx+1 with periodic boundary conditions
   sum_j += s[jx+ky*L];
   int jy = (ky + Lm1) % L;
   sum_j += s[kx+jy*L];
   jy = (ky+1) % L;
   sum_j += s[kx+jy*L];
   return sum_j;
}


int main(int argc, char* argv[])
{

   // These are the variables
   int L = 16;                     // Length along one side of a square lattice
   int N = L*L;                    // Total number of spins in the lattice
   std::vector<short> s;           // The spins of the system, defined in Eq. (1)
   double J = 1.;                  // Exchange energy, defined in Eq. (1)
   double H = 0.;                  // Zemann Field, defined in Eq. (1)
   double T = 2;                   // The temperature of the system
   int    M = 0;                   // Net magnetization, defined in Eq. (2)
   double m = 0;                   // net magnetization per spin, defined in Eq. (9)
   double E = 0;                   // The energy of the system
   long   imcs = 0;                // Total number of Monte Carlo steps taken (per T)
   double MCSS = 0;                // Monte Carlo steps per spin (size-independent MC time)
   long   imcs_tot = 0;            // Total number of Monte Carlo steps taken for program
   int    Nmeas = 10000;           // Number of samples to average together
   int    Nstep = 10;              // Number of MCSS between measurements
   int    Ntherm = 100;            // Number of measurements to throw away while system comes to equilibrium
   int    NumTemp = 100;           // Number of temperture points to simulate
  
   // Set the system size L 
   std::cout << "What is the system size L? " << std::endl;
   std::cin >> L;
   N = L*L;
   std::cout << "Simulating " << L << "x" << L << " Ising Model" << std::endl;

   // This is c++11
   // Note that L cannot be changed after krng declared
   std::random_device rd;                          // Let the machine generate random seeds
   std::mt19937 igen(rd());                        // An integer random number generator 
   std::uniform_real_distribution<> urng;          // Converts integer random numbers to real random numbers [0,1)
   std::uniform_int_distribution<>  krng(0,N-1);   // Randomly chooses a site in the lattice

   std::ofstream fout("Time.csv");
   fout << "#MC Simulation of " << L << "x" << L << " Ising Model" << std::endl;
   fout << "#Column 1: imcs, total Monte Carlo steps" << std::endl;
   fout << "#Column 2: MCSS, Monte Carlo steps per spin" << std::endl;
   fout << "#Column 3: T, Temperature in energy units" << std::endl;
   fout << "#Column 4: E, Energy of the system" << std::endl;
   fout << "#Column 5: M, net magnetization of the system" << std::endl;
   fout << "#Column 6: m=M/N, net magnetization per spin" << std::endl;

   std::ofstream aout("Average.csv");
   aout << "# MC Simulation of " << L << "x" << L << " Ising Model" << std::endl;
   aout << "# Column 1: T temperature in energy units" << std::endl;
   aout << "# Column 2: <|M|>/N" << std::endl;
   aout << "# Column 3: <M^2>/N" << std::endl;
   aout << "# Column 4: chi_m/N, Eq. (57b)" << std::endl;
   aout << "# Column 5: <E>/N" << std::endl;
   aout << "# Column 6: <E^2>/N" << std::endl;
   aout << "# Column 7: C_H/N, Eq. (57a)" << std::endl;
   aout << "# Column 8: number of samples" << std::endl;

   std::ofstream mout("Movie.txt");
   std::ofstream mout2("Movie2.txt");

   // Keep timing information to compare different example codes
   time_t start_time;
   time(&start_time);

   for(int itemp=0; itemp<100; itemp++)
   {
      // Set the temperture
      T = 1.3 + 2.*static_cast<double>(itemp)/static_cast<double>(NumTemp);   // T = [1.3,3.3)
      std::cout << "itemp=" << itemp << " T=" << T << std::endl;

      // Initialize
      s.resize(N);                              // Set the size of the vector of spins
      for(int k=0; k<s.size(); k++)          // Randomly assign up/down with equal probability
         s[k] = (urng(igen)>0.5)? 1 : -1;    // like the paper 
      double sum_kj = 0;
      for(int k=0; k<s.size(); k++)
      {
         int sum_j = calc_sumj(k,s,L);          // sum_j s[j]
         sum_kj += s[k]*sum_j;                  // sum_j s[k]*s[j]
      }                                         // after loop, sum_kj = sum_<kj> s[k]*s[j]
      E = -(J/2.)*sum_kj;                       // Eq. (1)
      M = 0;
      for(int k=0; k<s.size(); k++) M += s[k];               // Eq. (2);
      m = static_cast<double>(M)/static_cast<double>(N);
      imcs = 0;
      MCSS =0;

      double msum[5];                           // Used to accumulate moments of m for averaging
      for(int i=0; i<5; i++) msum[i]=0;
      double esum[5];                           // Used to accumulate moments of e for averaging
      for(int i=0; i<5; i++) esum[i]=0;

      // Do the simulation
      for(int imeas=0; imeas<(Nmeas+Ntherm); imeas++)
      {
         for(int istep=0; istep<(Nstep*N); istep++)
         {
            // 1. Choose randomly one site k in the lattice
            int k = krng(igen); 
            // 2. Calculate the energy difference DeltaEk between the actual energy and
            //    the energy if the spin is flipped
            // DeltaEk = Ek(n) -  Ek(m), where {m}->{n} is the move, see Eq. (49)
            // Ek(m) = -sum_j J*s[k,{m}]*s[j] = -J*(+s[k,{m}])*sum_j s[j]
            // Ek(n) = -sum_j J*s[k,{n}]*s[j] = -J*(-s[k,{m}])*sum_j s[j]
            // DeltaEk = J*s[k]*sum_j - -J*s[k]*sum_j = 2*J*s[k]*sum_j
            int sum_j = calc_sumj(k,s,L);
            double DeltaEk = 2.*J*s[k]*sum_j;
            // 3. If DeltaEk<0 accept the new configuration. Otherwise, accept it with probability
            //    exp(-beta*DeltaEk) where beta = 1/kBT
            bool accept = (DeltaEk<0);
            if( !accept )
            {
               double alpha = std::exp(-DeltaEk/T);    // Eq. (49)
               accept = (urng(igen)<alpha);
            }
            if( accept )
            {
               E += DeltaEk;
               M -= 2*s[k];
               s[k] *= -1;
            }
            imcs++;
            imcs_tot++;
         }
         MCSS = static_cast<double>(imcs)/static_cast<double>(N);
         // Accumulate information for averages
         if( imeas>=Ntherm )
         {
             m = static_cast<double>(M)/static_cast<double>(N);
             fout << imcs << " " << MCSS << " " << T << " " << E << " " << M << " " << m << std::endl;
             double mi = 1;             // Accumulate i-th moments of m in an array
             for(int i=0; i<5; i++)
             {
                msum[i] += mi;
                mi *= std::fabs(M);
             }
             double ei = 1;             // Accumulate i-th moments of e in an array
             for(int i=0; i<5; i++)
             {
                esum[i] += ei;
                ei *= E;
             }
         }
         // Output pictures of relaxation at first temperature
         if(itemp==0)
         {
            mout2 << "t=" << imcs << " mcs, " << MCSS << " MCSS" << std::endl;
            for(int iy=0; iy<L; iy++)
            {
               for(int ix=0; ix<L; ix++)
               {
                  // int ky = k/L; 
                  // int kx = k%L;         
                  int k = ix+L*iy;
                  if(s[k]>0) mout2 << "."; else mout2 << "*";
               }
               mout2 << std::endl;
            }
         }
      }
      // Output the averages for this temperature
      double rN = static_cast<double>(N);
      double mave = msum[1]/msum[0];         //   <m> = sum_i m / sum_i 1
      double m2ave = msum[2]/msum[0];        // <m^2> = sum_i m^2 / sum_i 1
      double chim = (m2ave-mave*mave)/T;     // Eq. (57b)
      double eave = esum[1]/esum[0];
      double e2ave = esum[2]/esum[0];
      double CH = (e2ave-eave*eave)/T/T;     // Eq. (57a)
      aout << T << " " << mave/rN << " " << m2ave/rN << " " << chim/rN << " " 
           << eave/rN << " " << e2ave/rN << " " << CH/rN << " " << msum[0] << std::endl;
      // Print a picture of the configuration
      mout << "T=" << T << " M=" << M << " " << " E=" << E << std::endl;
      for(int iy=0; iy<L; iy++)
      {
         for(int ix=0; ix<L; ix++)
         {
            // int ky = k/L; 
            // int kx = k%L;         
            int k = ix+L*iy;
            if(s[k]>0) mout << "."; else mout << "*";
          }
          mout << std::endl;
       }
    }
    fout.close();
    aout.close();
    mout.close();
 
    // Output information about execution time
    double psec = static_cast<double>(clock())/static_cast<double>(CLOCKS_PER_SEC);   // process clock
    time_t stop_time;
    time(&stop_time);
    double secs = difftime(stop_time,start_time);                                     // wall clock
    double mcsec = static_cast<double>(imcs_tot)/secs;
    if(secs<=0) mcsec = static_cast<double>(imcs)/psec;
    std::cout << "Total process time = " << psec << " seconds" << std::endl;
    std::cout << "Wall time in simulation = " << secs << " seconds" << std::endl;
    std::cout << "Monte Carlo Steps / second = " << mcsec << std::endl;

}
