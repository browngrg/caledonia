#ifndef WHAM_RANDOM_H_
#define WHAM_RANDOM_H_

#include<iostream>
#include<iomanip>
#include<limits>
#include<cmath>
#include<vector>
#include<time.h>

#ifdef USE_MPI
#include<mpi.h>
#endif


// A basic 32-bit LinearCongruential random number generator
class LinearCongruential32
{
public:
   typedef unsigned long INT32;              // long is at least 32 bits
   typedef unsigned long long INT64;         // long long is at least 64 bits
   static bool test_bits() { int bits32 = 8*sizeof(INT32); int bits64=8*sizeof(INT64); return (bits32==32) && (bits64>=(2*bits32)); }
public:
   LinearCongruential32(INT32 s=123456789) : a(1664525), c(1013904223) { seed(s); }
   void set_state(INT32 multiplier, INT32 increment) { a=multiplier; c=increment; }
   void seed(INT32 s) { x=s; }
   float operator()() { x = (a*x+c); return static_cast<float>(x)/static_cast<float>(std::numeric_limits<INT32>::max()); }
   float details()
   {
      //This gives output of intermediate results, specifically to look x_i as two 32-bit numbers (hi and lo)
      //The normal implementation just ignores the "carry" into hi
      //INT64 xi = a*x+c;
      //INT32 hi = hi = (xi>>32);  
      //INT32 lo = static_cast<INT32>(lo);
      //float  f = static_cast<float>(x)/static_cast<float>(0xFFFFFFFFUL);
      std::cout << "LC32: max=" << std::hex << std::numeric_limits<INT32>::max() << std::endl;
      std::cout << "  x_i= " << std::dec << x << std::endl;
      INT64 xi = a*x+c;
      std::cout << "x_i++= " << std::hex << xi << "\t" << std::dec << xi << std::endl;
      INT32 hi = (xi>>32);   // static_cast<INT32>(xi/0x100000000ULL);
      INT32 lo = static_cast<INT32>(xi);
      std::cout << "     = " << std::hex << hi << " " << lo << std::dec << std::endl;
      x = lo;
      std::cout << "x_i++= " << std::hex << x  << std::dec << std::endl;
      float f = static_cast<float>(x)/static_cast<float>(0xFFFFFFFFUL);
      std::cout << "    f= " << f << std::endl;
      return f; 
   }
private:
   INT32 x;
   INT64 a;  //  If a mod 8 is 1 or 5 and c is odd, the resulting base 232 congruential sequence will have period 232.
   INT32 c;  
};


// A lag-r multiply-with-carry random number generator, following en.wikipedia.org/wiki/Multiply-with-carry
// Includes as set of r+1 seed values x_0, x_1, x_2, ... x_r-1, and carry c_n-1. 
// The sequence is  x_n = (a*x_n-r + c_n-1) mod b
//                  c_n = floor[ ( a*x_n-r + c_n-1 ) / b ]
template<short r>
class MultiplyWithCarry
{
public:
   typedef unsigned long INT32;              // long is at least 32 bits
   typedef unsigned long long INT64;         // long long is at least 64 bits
   MultiplyWithCarry(INT32 s=1234566789) : a(18782) { this->seed(s); }
   void seed(INT32 s)
   {
      // initial carry
      x[0] = 362436 % a; 
      // LinearCongruential32
      INT32 xi = s;
      for(i=1; i<=r; i++) x[i] = (xi = 1664525*xi + 101390422);
      i=1;
   }
   void set_state(INT32 multiplier) { a=multiplier; x[0] = x[0] % a; }
   float operator()()
   {
      INT64 y = a*x[i] + x[0];
      x[0] = (y>>32);                    // The carry
#     if 0
      x[i] = static_cast<INT32>(y); 
#     else
      x[i] = ~static_cast<INT32>(y);      // The compliment makes it easier to calculate period
#     endif
      float f = static_cast<float>(x[i])/static_cast<float>(std::numeric_limits<INT32>::max());
      (i=(i%r))++;                       // 1 <= i <= r
      return f;
   }
   void print(std::ostream& out)
   {
     out << "MultiplyWithCarry<" << r << ">: a=" << a << " x=";
     for(int j=0; j<=r; j++) out << " " << x[j];
     out << std::endl;
   }
private:
   short i;          // index that treats x as a circular buffer
   INT64 a;          // the multiplier, needs to be 64-bit to force 64-bit precisionin y=a*x[i]+x[0]
   INT32 x[r+1];     // the seeds and the carry
};



void print_sizeof()
{
  std::cout << "sizeof(float)    =" << sizeof(float) << std::endl;
  std::cout << "sizeof(double)   =" << sizeof(double) << std::endl;
  std::cout << "sizeof(char)     =" << sizeof(char) << std::endl;
  std::cout << "sizeof(short)    =" << sizeof(short) << std::endl;
  std::cout << "sizeof(int)      =" << sizeof(int) << std::endl;
  std::cout << "sizeof(long)     =" << sizeof(long) << std::endl;
  std::cout << "sizeof(long long)=" << sizeof(long long) << std::endl;
  std::cout << "sizeof(u  )      =" << sizeof(unsigned int) << std::endl;
  std::cout << "sizeof(ul )      =" << sizeof(unsigned long) << std::endl;
  std::cout << "sizeof(ull)      =" << sizeof(unsigned long long) << std::endl;
}


// Runs simple moment tests of random number generators
template<class RNG>
bool RNGTestMoments(RNG& random, bool verbose, std::ostream& os) {

  bool failed = false;
  {
    const int nsamples   = 1000000;
    const int max_moment = 5;
    //
    // collect information about population
    //
    double momu[max_moment+1];
    for(int i=0; i<= max_moment; i++) { momu[i]=0.; }
    for(int i=0; i<nsamples; i++) {
      float urand = random();
      if((urand<0.)||(urand>=1.)) {
	failed = true;
	os << "RNGTestMoments failed, urand=" << urand << " outside [0,1)\n" << std::flush;
      }
      for(int j=0; j<= max_moment; j++) {
	momu[j]+=std::pow(urand,j);
      }
    }
    //
    // analyze the population
    //
    if(momu[0]<=2) {failed = true; return failed; }
    for(int j=0; j<=max_moment; j++) {
      double estimate = momu[j]/momu[0];
      double theory   = 1./static_cast<double>(j+1);
      double percentE = (theory>0)? std::abs(estimate-theory)/theory * 100. : 0.;
      if(verbose) {
	os << "   "
	   << "moment " << j 
	   << "\t estimate = " << estimate 
	   << "\t theory=" << theory
	   << "\t %error=" << percentE 
	   << std::endl << std::flush;
      }
      if(percentE>10.) {
	failed = true;
	os << "   ***RNGTestMoments failed for moment " << j << std::endl << std::flush;
      }
    }
  }
  return failed;
}


// Generate random seed value from system clock.
// This gets the time using an time_t = time(&time_t) and then does some 
// scrambling and puts the fastest varying bits first.
long int SeedFromClock() {
  // Get time information
  time_t rseed;
  rseed=time(&rseed);
  unsigned long int iraw=rseed;
  long int iseed = iraw;
  // Scramble that information
  long int high = iseed % 1000;
  long int low  = ( iseed / 1000 ) % 1000000;
  low  = ( ( low + high ) * 611953 ) % 1000000;
  iseed = high * 1000000 + low;
  //
  if(iseed<0) iseed *= -1;
  return iseed;
}


#ifdef USE_C11
#include<random>
#endif

// Give a unique, predictable seed to each processor
long ParallelSeed(long baseSeed)
{
   long seed_val = baseSeed;
#  ifdef USE_MPI
   // Get process id and global seed value
   int iproc = 0;
   int init = false;
   MPI_Initialized(&init);
   if( init )
   {
      MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
      MPI_Bcast(&seed_val,1,MPI_LONG,0,MPI_COMM_WORLD);
   }
   // set seed, get iproc-th result as my seed
   const int stride = 10;
#  ifdef USE_C11
   std::mt19937 mrand(seed_val);
   std::uniform_int_distribution<long> uint_dist;
   for(int icall=0; icall<stride*iproc; icall++)
      seed_val = uint_dist(mrand);
#  else
   LinearCongruential32 urng(seed_val);
   for(int icall=0; icall<stride*iproc; icall++)
      seed_val = urng()*std::numeric_limits<long>::max();
#  endif
#  endif
   return seed_val;
}


bool RNGTestParallelSeed(bool verbose=false)
{
   bool failed = false;
   long iseed = 20140806; // A random date
   iseed = ParallelSeed(iseed);
   int iproc = 0;
   int nproc = 1;
   std::vector<long> seeds;
#  ifdef USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&iproc);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   if(iproc==0) seeds.resize(nproc);
   MPI_Gather(&iseed,1,MPI_LONG,&(seeds[0]),1,MPI_LONG,0,MPI_COMM_WORLD);
#  else
   seeds.resize(1);
   seeds[0] = iseed;
#  endif
   if( verbose && iproc==0 )
   {
      std::cout << "RNGTestParallelSeed: iproc, seed" << std::endl;
      for(int i=0; i<seeds.size(); i++)
         std::cout << i << " " << seeds[i] << std::endl;
   }
   return failed;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Examples of constructing objects and calling tests
////////////////////////////////////////////////////////////////////////////////////////////////////
#if 0
int main(int argc, char* argv[])
{
   print_sizeof();
   std::cout << std::endl << std::endl;

   std::cout << "LinearCongruential32:" << std::endl;
   LinearCongruential32 urng_lc32;
   std::cout << "32bits = " << urng_lc32.test_bits() << std::endl;
   for(int i=0; i<10; i++) std::cout << urng_lc32.details() << std::endl;
   for(int i=0; i<10; i++) std::cout << urng_lc32() << std::endl;
   RNGTestMoments(urng_lc32,true,std::cout);
   std::cout << std::endl << std::endl;

   std::cout << "MultiplyWithCarry-1:" << std::endl;
   MultiplyWithCarry<1> urng_mwc1;
   for(int i=0; i<10; i++) std::cout << urng_mwc1() << std::endl;
   RNGTestMoments(urng_mwc1,true,std::cout);
   std::cout << std::endl << std::endl;

   std::cout << "MultiplyWithCarry-4:" << std::endl;
   MultiplyWithCarry<4> urng_mwc4;
   for(int i=0; i<10; i++) std::cout << urng_mwc4() << std::endl;
   RNGTestMoments(urng_mwc1,true,std::cout);

}
#endif  // Omit main



#endif
