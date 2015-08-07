#ifndef STATREC_HPP_
#define STATREC_HPP_

#include<iostream>
#include<vector>
#include<cmath>

// StatRec.hpp -- Class that turns accumulated moments into processed statistics
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)
// October 3, 2014
//
// http://en.wikipedia.org/wiki/Central_moment
// mean mu = mup = E[X]
// central moment mu_n = E[ (X - E[X])^n ]
// noncentral moment mup_n = E[ X ]
// mu_2 = mup_2 - mu^2
// mu_3 = mup_3 - 3*mu*mup_2 + 2*mu^3
// mu_4 = mup_4 - 4*mu*mup_3 + 6*mu^2*mup_2-3*mu^4
// http://en.wikipedia.org/wiki/Cumulant
// cumulant 
// kappa_1 = mup_1
// kappa_2 = mu_2
// kappa_3 = mu_3
// kappa_4 = mu_4 - 3*mu_2^2
// Formulas for cumulants in terms of mup
// kappa_2 = mup_2 - mup_1^2
// kappa_3 = mup_3 - 3*mup_2*mup_1 + 2*mup_1^3
// kappa_4 = mup_4 - 4*mup_3*mup_1 - 3*mup2_^2 + 12*mup_2*mup_1^2-6*mup_1^4
// Moments in terms of cumulants
// mup_4 = kappa_4 + 4*kappa_3*kappa_1 + 3*kappa_2^2 + 6*kappa_2*kappa_1^2 + kappa_1^4
// mu_4 = kappa_4 + 3*kappa_2^2
//
// http://en.wikipedia.org/wiki/Kurtosis
// Excess kurtosis is gamma_2 = kappa_4/kappa_2^2 = mu_4/mu_2^2 - 3
//        kurtosis is beta_2  = mu_4/mu_2^2
//        skewness is         = mu_3/(mu_2)^(3/2)
//        bound beta_2 >= skew^2 + 1
//
struct StatRec
{
public:
   double mu,mu_2,mu_3,mu_4;  // mean and central moments
   double mup_2,mup_3,mup_4;  // noncentral moments
   double kappa_4;            // statistics cumulant
   double binder_4;           // Binder cumulant 1 - mu_4/(3*mu_2^2)
   double skew,kurt,excess;   // Skewness, kurtosis (beta_2), and excess curtosis (gamma_2)
public:
   StatRec(const std::vector<double>& xmom) { this->operator()(xmom); }
   void operator()(const std::vector<double>& xmom)
   {
      mup_2=mup_3=mup_4=mu_2=mu_3=mu_4=kappa_4=binder_4=skew=kurt=0;
      switch( xmom.size() )
      {
         default: mup_4 = xmom[4]/xmom[0];
         case 4:  mup_3 = xmom[3]/xmom[0];
         case 3:  mup_2 = xmom[2]/xmom[0];
         case 2:  mu = xmom[1]/xmom[0];
                  break;
         case 1:
         case 0:  return;
      }
      double mu2th = mu*mu;
      double mu4th = mu2th*mu2th;
      switch( xmom.size() )
      {
         default: mu_4 = mup_4 - 4*mu*mup_3 + 6*mu2th*mup_2 - 3*mu4th;
         case 4:  mu_3 = mup_3 - 3*mu*mup_2 + 2*mu2th*mu;
         case 3:  mu_2 = mup_2 - mu2th;
         case 2:
         case 1:
         case 0:
                  break;
      }
      if( xmom.size()>4 )
      {
         kappa_4 = mu_4 - 3*mu_2*mu_2;   
         binder_4 = 1 - mu_4/(3*mu_2*mu_2);
         skew = mu_3/std::pow(mu_2,1.5);
         kurt = mu_4/mu_2/mu_2;
         excess = kappa_4/mu_2/mu_2;
         if( false )
         {
            bool ok = true;
            // kappa_2 = mup_2 - mup_1^2
            double pkappa_2 = mup_2 - mu*mu;
            double diff = std::abs(pkappa_2-mu_2)/std::abs(pkappa_2);
            if( diff>1.e-5 )
            {
               std::cerr << "StatRec error in kappa_2, diff=" << diff << std::endl;
               ok = false;
            }
            // kappa_3 = mup_3 - 3*mup_2*mup_1 + 2*mup_1^3
            double pkappa_3 = mup_3 - 3*mup_2*mu + 2*mu*mu*mu;
            diff = std::abs(pkappa_3-mu_3)/std::abs(pkappa_3);
            if( diff>1.e-5 )
            {
               std::cerr << "StatRec error in kappa_3, diff=" << diff << std::endl;
               ok = false;
            }
            // kappa_4 = mup_4 - 4*mup_3*mup_1 - 3*mup2_^2 + 12*mup_2*mup_1^2-6*mup_1^4
            double pkappa_4 = mup_4 - 4*mup_3*mu - 3*mup_2*mup_2 + 12*mup_2*mu*mu - 6*mu4th;
            diff = std::abs(pkappa_4-kappa_4)/std::abs(pkappa_4);
            if( diff>1.e-5 )
            {
               std::cerr << "StatRec error in kappa_4, diff=" << diff << std::endl;
               ok = false;
            }
            if( mu_2<-1e-6 )
            {
               std::cerr << "StatRec error in mu_2, val=" << mu_2 << std::endl;
               ok = false;
            }
            if( mu_4<-1e-6 )
            {
               std::cerr << "StatRec error in mu_4, val=" << mu_4 << std::endl;
               std::cerr << "terms = " << mup_4 << " " << -4*mu*mup_3 << " " << 6*mu*mu*mup_2 << " " << -3*mu*mu*mu*mu << std::endl;
               std::cerr << "nested mults = " << mup_4 + mu*(-4*mup_3+mu*(6*mup_2-3*mu*mu)) << std::endl;
               std::cerr << "mu_4 = kappa_4 + 3*kappa_2^2 = " << kappa_4 + 3*mu_2*mu_2 << std::endl;
               ok = false;
            }
            if( binder_4>1 )
            {
               std::cerr << "StatRec error in binder_4, val=" << binder_4 << std::endl;
               std::cerr << " = 1 - " << (mu_4) << "/(  3*(" << mu_2 << ")^2 )" << std::endl;;
               ok = false;
            }
            if( !ok )
            {
               std::cerr << "mu=" << mu << std::endl
                    << "mu_2=" << mu_2 << std::endl
                    << "mu_3=" << mu_3 << std::endl
                    << "mu_4=" << mu_4 << std::endl
                    << "mup_2=" << mup_2 << std::endl
                    << "mup_3=" << mup_3 << std::endl
                    << "mup_4=" << mup_4 << std::endl
                    << "kappa_4=" << kappa_4 << std::endl
                    << "binder_4=" << binder_4 << std::endl
                    << "skew=" << skew << std::endl
                    << "kurtosis=" << kurt << std::endl
                    << "pkappa_2=" << pkappa_2 << ", diff=" << pkappa_2 - mu_2 << std::endl
                    << "pkappa_3=" << pkappa_3 << ", diff=" << pkappa_3 - mu_3 << std::endl
                    << "pkappa_4=" << pkappa_4 << ", diff=" << pkappa_4 - kappa_4 << std::endl;
            }
         }
      }
   }

   // adds a sample to accumulated moments array, does not affect StatRec object
   static void add_sample(double x, std::vector<double>& xmom)
   {
      int kmax=xmom.size();
      //if( kmax==0 ) { xmom.resize(5); for(int k=0; k<5; k++) xmom[k]=0; }
      double xk = 1;
      for(int k=0; k<kmax; k++)
      {
         xmom[k] += xk;
         xk *= x;
      }
   }

   // adds a histogram point to accumulated moments array, does not affect StatRec object
   static void add_histpt(double x, double pdf, std::vector<double>& xmom)
   {
      int kmax=xmom.size();
      //if( kmax==0 ) { xmom.resize(5); for(int k=0; k<5; k++) xmom[k]=0; }
      double xk = 1;
      for(int k=0; k<kmax; k++)
      {
         xmom[k] += pdf*xk;
         xk *= x;
      }
   }
};


#endif
