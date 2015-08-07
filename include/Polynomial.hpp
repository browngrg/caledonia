#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_ 


#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<cmath>


// Taylor-expansion polynomials p(x) = sum_i a_i x^i
class Polynomial : public std::vector<double> {

public:

   Polynomial(int degree=0) : std::vector<double>(degree+1,0) { std::fill(this->begin(),this->end(),0); }

   Polynomial(const std::vector<double>& coeff) : std::vector<double>(coeff) {;}

   const std::vector<double>& operator=(const std::vector<double>& coeff) 
   { 
      this->resize(coeff.size());
      for(int i=0; i<coeff.size(); i++) (*this)[i] = coeff[i];
      return coeff;
   }

   // Evaluate the polynomial at x
   double operator()(double x) const { return value(*this,x); }

public:

   // Find the value of the polynomial at x
   static double value(const std::vector<double>& coeff, double x) {
      double val = 0;
      double xk = 1; for(int k=0; k<coeff.size(); k++) { val += coeff[k]*xk; xk *= x; }
      return val;
   }

   // Find the derivative of the polynomial
   static std::vector<double> derivative(const std::vector<double>& coeff)
   {
      if(coeff.size()<2) return std::vector<double>();
      std::vector<double> d(coeff.size()-1);
      for(int k=1; k<coeff.size(); k++)
         d[k-1] = k*coeff[k];
      return d;
   }

   // Find the antiderivative of the polynomial
   static std::vector<double> antiderivative(const std::vector<double>& coeff) 
   {
      std::vector<double> a(coeff.size()+1,0);
      for(int k=1; k<a.size(); k++)
         a[k] = coeff[k-1]/static_cast<double>(k);
      a[0] = 0;
      return a;
   }

   // Find the integral of polynomial over range
   static double integral(const std::vector<double>& coeff, double a, double b)
   {
      int n = degree(coeff)+1;
      double I = 0;
      double aip1 = 1;
      double bip1 = 1;
      for(int i=0; i<n; i++) {
         aip1 *= a;
         bip1 *= b;
         I += coeff[i]/static_cast<double>(i+1) * (bip1 - aip1);
      }
      return I;
   }

   // Find the degree of a polynomial
   static int degree(const std::vector<double>& coeff)
   {
      int degree = coeff.size()-1;
      while( degree>=0 && coeff[degree]==0) degree--;
      return degree; 
   }

   // Find the product of two polynomials
   static std::vector<double> multiply(const std::vector<double>& a, const std::vector<double>& b)
   {
      int adeg = degree(a);
      int bdeg = degree(b);
      int cdeg = adeg + bdeg;
      std::vector<double> c(cdeg+2,0);
      for(int ia=0; ia<=adeg; ia++)
         for(int ib=0; ib<=bdeg; ib++)
            c[ia+ib] += a[ia]*b[ib];
      return c;
   }

   // Find the root of the polynomial
   static double root(const std::vector<double>& coeff, double guess)
   {
      int deg = degree(coeff);
      if( deg<1 ) return guess;              // constant function
      if( deg==1 )return -coeff[0]/coeff[1]; // line
      Polynomial f(coeff);
      Polynomial df = derivative(f);
      double x = guess;
      double h = 0;
      if( df(x)!=0 ) h = -f(x)/df(x);
      int i = 0;
      while( std::fabs(h)>1.e-12 && i<10000 ) {
         x += h;
         if( df(x)!=0 ) h = -f(x)/df(x);
         i++;
      }
      return x;
   }

   // Find a value of independent variable where polynomial has derivative of zero
   static double extremum(const std::vector<double>& coeff, double xguess=0)
   {
      if( degree(coeff)<3 ) return xguess;
      Polynomial df = derivative(coeff);
      return root(df,xguess);
   }


   friend std::ostream& operator<<(std::ostream&,const Polynomial&);
   friend std::istream& operator>>(std::istream&,const Polynomial&);

};

std::ostream& operator<<(std::ostream& out,const Polynomial& p) 
{ out << p.size(); for(int k=0; k<p.size(); k++) out << " " << p[k]; return out; }

std::istream& operator<<(std::istream&  in, Polynomial& p) 
{ int ps; in >>ps; p.resize(ps); for(int k=0; k<p.size(); k++) in >> p[k]; return in;}

Polynomial derivative(const Polynomial& p) { return Polynomial::derivative(p); }

Polynomial antiderivative(const Polynomial& p) { return Polynomial::antiderivative(p); }

double root(const Polynomial& p, double guess=0) { return Polynomial::root(p,guess); }

double extremum(const Polynomial& p, double xguess=0) 
{ Polynomial dp = Polynomial::derivative(p); return Polynomial::root(dp,xguess); }




// when setting t[0] by hand, need to multiply by 2, since here T_0 = 1/2
// (as opposed to the usual T_0 = 1;
class Chebyshev {

private:

   double xm,xp,hdx;

public:

   // Specify the domain x \in [a,b] the polynomials describe
   void set_domain(double a, double b)
   {
      xm = a;
      xp = b;
      if( xm>xp ) { double tmp=xm; xm=xp; xp=tmp; }
      hdx = (xp-xm)/2.;
   }

   // Convert x \in [a,b] to c \in [-1,+1]
   double scale(double x) const { return (x-xm)/hdx - 1.; }

   // Convert c \in [-1,+1] to x \in [a,b]
   double unscale(double c) const { return (c+1.)*hdx + xm; }

public:

   // The coefficients of the expansion
   std::vector<double> t;

   // Default constructor: a[0] = 0
   Chebyshev(double a=-1, double b=+1) 
   { 
      set_domain(a,b); 
      this->t.resize(1,0); 
   }

   // Construct with given coefficients
   Chebyshev(int n, double* aptr, double a=-1, double b=+1)
   {
      set_domain(a,b);
      t.resize(n);
      for(int i=0; i<n; i++) t[i]=aptr[i];
   }

   // Copy constructor
   Chebyshev(const Chebyshev& orig) { copy(orig); }

   // Assignment operator
   const Chebyshev& operator=(const Chebyshev& orig) { copy(orig); return orig; }

   void copy(const Chebyshev& orig) 
   { 
      xm=orig.xm; xp=orig.xp; hdx=orig.hdx; t=orig.t; 
   }

   // Value of the polynomial at x \in [a,b]
   double value(double x) const
   {
      // This is the correct 5/17/2013
      double y = scale(x);
      // Clenshaw's algorithm
      int n = t.size() - 1;
      double bkp2 = 0;
      double bkp1 = 0;
      double bk   = 0;
      for(int k=n; k>=0; k--)
      {
         bkp2 = bkp1;
         bkp1 = bk;
         bk = 2*y*bkp1 - bkp2 + t[k];
      }
      return (bk-bkp2+t[0])/2.;
   }

   double operator()(double x) const { return value(x); }

   // Return polynomial describing product of two polynomials
   static Chebyshev multiply(const Chebyshev& r, const Chebyshev& q)
   {
      //  Fundamental property is 2*T_r(x)*T_q(x) = T_{r+q}(x) + T_{|r-q|}(x)
      //  Assume equivalent domain
      Chebyshev m(r);
      m.t.resize( r.t.size()+q.t.size()-1 );
      for(int k=0; k<m.t.size(); k++) m.t[k]=0;
      for(int i=0; i<r.t.size(); i++) {
         for(int j=0; j<q.t.size(); j++) {
            m.t[i+j] += r.t[i]*q.t[j];
            int aimj = i-j; if(aimj<0) aimj*=-1;
            m.t[aimj] += r.t[i]*q.t[j];
         } 
      }
      for(int i=0; i<m.t.size(); i++) m.t[i] /= 2.;
      return m;
   }


   // Returns Chebyshev polynomial for the derivative of the polynomial
   Chebyshev derivative() const
   {
      Chebyshev df(xm,xp);
      df.t.resize(t.size()-1);
      for(int k=(df.t.size()-1); k>=0; k--)
      {
         if( (k+1)<t.size() )
         {
            df.t[k] = static_cast<double>(2*(k+1))*t[k+1];
         }
         else
         {
            df.t[k] = 0;
         }
         if((k+2)<(df.t.size())) df.t[k] += df.t[k+2];
      }
      for(int k=0; k<df.t.size(); k++) df.t[k] /= hdx;
      df.t[0] /= 2.;
      return df;
   }

   // Returns Chebyshev polynomial for the antiderivative of the polynomial
   Chebyshev antiderivative() const
   {
      // integral Tk dx = (1/2) ( T_{k+1}/(k+1) - T_{k-1}/(k-1) )
      // F = sum b_j Tj = integral sum a_k Tk
      // F = sum b_j Tj = sum a_k integral Tk = sum a_k ( T_{k+1}/(k+1) - T_{k-1}/(k-1) )
      // gather terms with like Tj
      //         b_j Tj =  a_{j-1}/2 Tj/j - a_{j+1}/2 Tj/j
      //      2j b_j    =  a_{j-1} - a_{j+1}
      // This agrees with numerical recipes: C_i = ( c_{i-1} - c_{i+1} ) / (2*i)
      Chebyshev fi(xm,xp);
      fi.t.resize(t.size()+1);
      fi.t[1] = 2*t[0];                    // t[0] is always an issue
      for(int k=1; k<t.size(); k++)
      {
         fi.t[k-1] -= t[k];
         fi.t[k+1]  = t[k];
      }
      for(int k=1; k<fi.t.size(); k++) fi.t[k] *= hdx/static_cast<double>(2*k);
      fi.t[0] = 0;
      return fi;
   }

   // Build antiderivative brute force (used for debugging)
   Chebyshev antiderivative2() const
   {
      std::vector<double> x(201);
      for(int i=0; i<x.size(); i++)
         x[i] = unscale( -1. + 2.*static_cast<double>(i)/200. );
      std::vector<double> y(x.size());
      y[0] = this->value(x[0]);
      for(int i=1; i<x.size(); i++)
      {
         double xm = (x[i]+x[i-1])/2.;
         y[i] = y[i-1] + (x[i]-x[i-1])*this->value(xm);
      }
      Chebyshev fi;
      fi.t.resize(t.size()+1);
      fi.fit(x,y,xm,xp);
      if( true ) 
      {
         Chebyshev df = fi.derivative();
         Chebyshev ad = this->antiderivative();
         Chebyshev da = ad.derivative();
         std::cout << "# brute: " << fi << std::endl;
         std::cout << "# libry: " << ad << std::endl;
         std::cout << "# x, euler, brute, libry, f, de/dx, db/dx, dl/dx" << std::endl;
         for(int i=1; i<x.size(); i++)
         {
            double xi = x[i];
            double de = (y[i]-y[i-1])/(x[i]-x[i-1]);
            std::cout << xi << " " << y[i]-y[0] << " " << fi(xi)-fi(x[0]) << " " << ad(xi)-ad(x[0]) << " " << value(xi) << " " << de << " " << df(xi) << " " << da(xi) << std::endl;
         }
      }
      return fi;
   }

   // Returns the definite integral of the polynomial
   double integral(double a, double b) const
   {
      Chebyshev F = antiderivative();
      return F(b) - F(a);
   }

   // Calculate the Taylor-expansion coefficients for a Chebyshev polynomial
   std::vector<double> taylor(int N=7) const
   {
      std::vector<double> taylor(N,0);
      taylor[0] = t[0];
      taylor[1] = t[1];
      // based on code tested in ChebyshevTaylor.cc
      // static int chebunit[4] = { +1, 0, -1, 0 };
      std::vector<int> cheb(N,0);
      std::vector<int> chebm1(N,0);
      std::vector<int> chebm2(N,0);
      chebm2[0] = 1;   // T_0(x) = 1
      chebm1[1] = 1;   // T_1(x) = x
      for(int j=2; j<N; j++)
      {
         // T_j(x)
         cheb[0] = -chebm2[0];
         for(int i=1; i<=j; i++)    
            cheb[i] = 2*chebm1[i-1] - chebm2[i];  // T_j(x) = 2x*T_{j-1}(x) - T_{j-1}(x)
         chebm2 = chebm1;
         chebm1 = cheb;
         // Transfer results to new polynomial
         // At this point "cheb" is an array of the general coefficients for T_j(x) = sum_i cheb[i]*x^i
         // and t[j] are the particular expansion coefficients for this object f(x) = sum_j t[j] T_j(x)
         // taylor[0] += chebunit[j%4]*t[j];
         for(int i=0; i<=j; i++)
            taylor[i] += cheb[i]*t[j];
      } 
      return taylor;
   }

   // Find the root of the polynomial
   double root(double guess) const
   {
      int deg = Polynomial::degree(t);
      if( deg<1 ) return guess;              // constant function
      if( deg==1 )return -t[0]/t[1];         // line
      const Chebyshev& f(*this);
      Chebyshev df = derivative();
      double x = guess;
      double h = 0;
      if( df(x)!=0 ) h = -(x)/df(x);
      int i = 0;
      while( std::fabs(h)>1.e-16 && i<10000 ) {
         x += h;
         if( df(x)!=0 ) h = -f(x)/df(x);
         i++;
      }
      return x;
   }

   // Find a value of independent variable where polynomial has derivative of zero
   double extremum(double xguess=0) const
   {
      if( Polynomial::degree(t)<3 ) return xguess;
      Chebyshev df = derivative();
      return df.root(xguess);
   }

   // Find the smallest value of the polynomial
   double minimum()
   {
      Chebyshev df = derivative();
      // adjust the maximum iterations and step size if missing extrema
      double xguess = xm;
      double minval = value(xguess);
      double minloc = xm;
      double xlast  = xm;
      int iter = 0;
      const int niter = 1000;
      const bool accelguess = true;  
      while( xguess<xp && iter<niter ) 
      {
         if( !accelguess )
            xguess = xm + (xp-xm)*static_cast<float>(iter)/static_cast<float>(niter);
         double xnext = df.root(xguess);
         double xdiff2 = (xlast-xnext)*(xlast-xnext)/hdx/hdx;
         if( xdiff2>1.e-12 && xnext>=xm && xnext<=xp )
         {
            xlast = xnext;
            if( accelguess && xnext>xguess ) xguess=xnext;
            double xval = value(xnext);
            if( xval<minval ) 
            {
               minval = xval;
               minloc = xnext;
            }
         }
         xguess += 0.0001*hdx;
         iter++;
      } 
      return minval;
   }


   // Calculate the coefficients of Chebyshev Polynomial base on an interpolation table
   void fit(const std::vector<double>& xvec, const std::vector<double>& yvec, double a=-1, double b=+1)
   {
      int ncoeff = t.size();
      if( xvec.size()<ncoeff ) ncoeff = xvec.size();
      set_domain(a,b);
      t.resize(ncoeff);
      for(int j=0; j<ncoeff; j++)
      {
         double tj = 0;
         for(int k=0; k<ncoeff; k++)
         {
            double arg = M_PI*static_cast<double>(2*k+1)/static_cast<double>(2*ncoeff);
            double x = unscale(std::cos(arg));
            int ipt = xvec.size() - 2;
            while( ipt>0 && x<xvec[ipt] ) ipt--;  // k increasing == x decreasing
            double y = yvec[ipt] + (yvec[ipt+1]-yvec[ipt])*(x-xvec[ipt])/(xvec[ipt+1]-xvec[ipt]);
            double fx = y;
            tj += fx*std::cos(static_cast<double>(j)*arg);
         }
         t[j] = tj*2./static_cast<double>(ncoeff);
      }
      t[0] /= 2.;
   }

   friend std::ostream& operator<<(std::ostream&,const Chebyshev&);
   friend std::istream& operator>>(std::istream&,const Chebyshev&);

};

std::ostream& operator<<(std::ostream& out,const Chebyshev& p) 
{ out << p.t.size(); for(int k=0; k<p.t.size(); k++) out << " " << p.t[k]; return out; }

std::istream& operator<<(std::istream&  in, Chebyshev& p) 
{ int ps; in >>ps; p.t.resize(ps); for(int k=0; k<p.t.size(); k++) in >> p.t[k]; return in;}





// Calculate the coefficients of Chebyshev Polynomial based on given polynomial
Chebyshev chebyshev(const std::vector<double>& coeff, double a=-1, double b=+1)
{ 
   int ncoeff = coeff.size();
   Chebyshev cheb(a,b);
   cheb.t.resize(ncoeff);
   for(int j=0; j<ncoeff; j++)
   {
      double tj = 0;
      for(int k=0; k<ncoeff; k++)
      {
         double arg = M_PI*static_cast<double>(2*k+1)/static_cast<double>(2*ncoeff);
         double x = cheb.unscale(std::cos(arg));
         double fx = Polynomial::value(coeff,x); 
         tj += fx*std::cos(static_cast<double>(j)*arg);
      }
      cheb.t[j] = tj*2./static_cast<double>(ncoeff);
   }
   cheb.t[0] /= 2.;
   return cheb;
}


// Calculate the coefficients of Chebyshev Polynomial base on an interpolation table
Chebyshev chebyshev_interpolate(const std::vector<double>& xvec, const std::vector<double>& yvec, double a=-1, double b=+1)
{
   Chebyshev cheb(a,b);
   cheb.fit(xvec,yvec,a,b);
   return cheb;
}


// Needed to match interface of generic polynomial
Chebyshev chebyshev(const Chebyshev& orig) { return orig; }


// Convert Chebyshev to new range
Chebyshev chebyshev(const Chebyshev& orig, double a, double b) 
{ 
   int ncoeff = orig.t.size();
   Chebyshev cheb(a,b);
   cheb.t.resize(ncoeff);
   for(int j=0; j<ncoeff; j++)
   {
      double tj = 0;
      for(int k=0; k<ncoeff; k++)
      {
         double arg = M_PI*static_cast<double>(2*k+1)/static_cast<double>(2*ncoeff);
         double x = cheb.unscale(std::cos(arg));
         double fx = orig(x); 
         tj += fx*std::cos(static_cast<double>(j)*arg);
      }
      cheb.t[j] = tj*2./static_cast<double>(ncoeff);
   }
   cheb.t[0] /= 2.;
   return cheb;
}


// Add a count to a chebyshev histogram
// i.e, tabulate Chebyshev integration of moments
// This uses the identity h(x) = sum_j c_j T_j(x) = p(x) = sum_i=1,N delta(x_i)
// and the orthogonality integral int_-1^+1 T_k(x) T_j(x) / sqrt(1-x^2) = pi delta_k,j; k>0
// To derive c_k = ([2|1])/pi * sum_i=1,N T_k(x_i)/sqrt(1-x_i^2)
void add_count(Chebyshev& cheb, double xval)
{
   double xwin = cheb.scale(xval);
   if( xwin<-1 || xwin>+1 ) return;
   double wcheb = 1./std::sqrt(1.-xwin*xwin)/M_PI; // M_PI gets the scale of histogram right in counts
   double chebxm1 = 1.;
// cheb.t[0] += 2.*wcheb;                          // This gives histograms that don't go to zero
   cheb.t[0] += 0.5*wcheb;                         // This gives histograms that go to zero
   double chebxm  = xwin;
   cheb.t[1] += wcheb*xwin;
   for(int i=2; i<cheb.t.size(); i++)
   {
      double chebx = 2*xwin*chebxm - chebxm1;
      chebxm1 = chebxm;
      chebxm  = chebx;
      cheb.t[i] += wcheb*chebx;
   }
}

Chebyshev scale_hist(const Chebyshev& cheb, double counts=1)
{
   Chebyshev scaled = cheb;
   double s = std::fabs(scaled.t[1]);
   if( s!=0 )
      for(int k=0; k<scaled.t.size(); k++) scaled.t[k] /= s;
   Chebyshev F = scaled.antiderivative();
for(int k=0; k<F.t.size(); k++) std::cout << " " << F.t[k];
std::cout << std::endl;
   F.t[0] = 0;
/*
   s = std::fabs(F.t[1]);
   if( s!=0 )
      for(int k=0; k<F.t.size(); k++) F.t[k] /= s;
*/
   double lo = scaled.unscale(-1);
   double hi = scaled.unscale(+1);
   double norm = counts/( F(hi) - F(lo) );
std::cout << "hist: " << norm << " " << counts << " " << F(hi) << " " << F(lo) << " " << lo << " " << hi << std::endl;
   if( norm!=0 )
      for(int k=0; k<scaled.t.size(); k++) scaled.t[k] /= norm;
   return scaled;
}

// Find the best Chebyshev for the natural logarithm of given polynomial
Chebyshev chebyshev_logpoly(const std::vector<double>& coeff, double a=-1, double b=+1)
{
   int ncoeff = coeff.size();
   Chebyshev cheb(a,b);
   cheb.t.resize(ncoeff);
   for(int j=0; j<ncoeff; j++)
   {
      double tj = 0;
      for(int k=0; k<ncoeff; k++)
      {
         double arg = M_PI*static_cast<double>(2*k+1)/static_cast<double>(2*ncoeff);
         double x = cheb.unscale(std::cos(arg));
         double fx = Polynomial::value(coeff,x); 
         if( fx>0 ) {
            double logfx = std::log(fx);
            tj += logfx*std::cos(static_cast<double>(j)*arg);
         }
      }
      cheb.t[j] = tj*2./static_cast<double>(ncoeff);
   }
   cheb.t[0] /= 2.;
   return cheb;
}


Chebyshev chebyshev_logpoly(Chebyshev poly)
{
   Chebyshev cheb(poly);
   int ncoeff = cheb.t.size();
   cheb.t.resize(ncoeff);
   for(int j=0; j<ncoeff; j++)
   {
      double tj = 0;
      for(int k=0; k<ncoeff; k++)
      {
         double arg = M_PI*static_cast<double>(2*k+1)/static_cast<double>(2*ncoeff);
         double x = cheb.unscale(std::cos(arg));
         double fx = poly(x); 
         if( fx>0 ) {
            double logfx = std::log(fx); 
            tj += logfx*std::cos(static_cast<double>(j)*arg);
         }
      }
      cheb.t[j] = tj*2./static_cast<double>(ncoeff);
   }
   cheb.t[0] /= 2.;
   return cheb;
}


// some values very small, use a least squares fit to other values
// fit to least squares fit of exponential (lnh = mx+b)
// This is specific to Wang-Landau pileup
Chebyshev chebyshev_logpoly2(Chebyshev poly)
{
   Chebyshev cheb(poly);
   int ncoeff = cheb.t.size();
   cheb.t.resize(ncoeff);
   double xm = poly.unscale(-1);
   double xp = poly.unscale(+1);
   double dx = 1.e-4*(xp-xm);
   double fxm = poly(xm);
   double fxp = poly(xp);
   double x = xm;
   double fx = fxm;
   if( fxp>fxm ) { dx *= -1; x = xp; fx = fxp; }
   double thresh = fx/10.;  // fit over a decade
   double S,Sx,Sy,Sxx,Sxy,Syy; 
   S=Sx=Sy=Sxx=Sxy=Syy=0;
   while( fx>thresh )
   {
      double logfx = std::log(fx);
      S+=1; Sx += x; Sy += logfx; Sxx += x*x; Sxy += x*logfx; Syy += logfx*logfx;
      x += dx;
      fx = poly(x);
   }
   double Delta = S*Sxx-Sx*Sx;
   double b = (Sxx*Sy-Sxy*Sx)/Delta;
   double m = (Sxy*S-Sx*Sy)/Delta;
   // fit polynomial to this line
   for(int i=2; i<cheb.t.size(); i++) cheb.t[i] = 0;
   cheb.t[0] = b;
   cheb.t[1] = m;
   return cheb;
}



#endif  /* POLYHOMIAL_H_ */
