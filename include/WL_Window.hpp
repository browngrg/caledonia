#ifndef WL_WINDOW_HPP_
#define WL_WINDOW_HPP_

// Ditch when don't need for debugging
#include<iostream>

#include<vector>

#ifdef USE_MPI
#include<mpi.h>
#endif

// Class WL_Window -- Details about the energy window each walker operates in
class WL_Window
{

public:
    // -------------------- Constructors and such -------------------- //

   // The constructor
   WL_Window(float Emin=0, float Emax=1, int iwin=0, int numbin=200) { set(Emin,Emax,iwin,numbin); }

   // Set the values
   void set(float Emin, float Emax, int iwin, int numbin);
   void set_numbr(float Emin, float Emax, int iwin, int numbin);
   void set_delta(float Emin, float Emax, int iwin, float Ebin);

   // Copy constructor
   WL_Window(const WL_Window& orig) { this->copy(orig); }

   // Copy a window
   void copy(const WL_Window& orig);

   // Assignment operator
   const WL_Window& operator=(const WL_Window& orig) { copy(orig); return orig; }

public:
   //-------------------- Binning within Window --------------------//

   // convert E into bin index on [0,NBin-1]
   int bin(float E) const;

   // return value along scaled coordinate of bin
   // ubin=0 is the left, ubin=0.5 the center, ubin=1 the right of bin
   // Elo = unbin(0,0.) while Ehi = unbin(NBin,0.) = unbin(NBin-1,1.)
   //float unbin(int ibin, float ubin=0.5) const { return Elo+(ibin+ubin)*Ebin; }
   // 2/3/2015 fixed ubin = 0.5
   float unbin(int ibin) const { return Elo+(ibin+0.5)*Ebin; }

   // test if point is in range (excluding buffer region)
   bool in(float E) const { return E>=Elo && E<=Ehi; }

   // Return a container of bin centers
   void energies(std::vector<double>& energy) const { energy.resize(NBin); for(int i=0; i<NBin; i++) energy[i]=this->unbin(i); }

public: 
   //--------------- Implementation details ---------------//

   float Elo,Ehi;                      // Limits of the Window
   float Ebin;                         // Energy width of a Wang-Landau bin
   int   NBin;                         // Number of Wang-Landau Bins
 
   int iwindow;                        // Index for this window

};


void WL_Window::set(float Emin, float Emax, int iwin, int numbin)
{
   this->set_numbr(Emin,Emax,iwin,numbin); // Legacy call
}

void WL_Window::set_numbr(float Emin, float Emax, int iwin, int numbin)
{
   Ebin = (Emax-Emin)/static_cast<float>(numbin);
   Elo = Emin; // - 0.5*Ebin;
   Ehi = Emin + numbin*Ebin;
   NBin = numbin; // +1;
   iwindow = iwin;
}

void WL_Window::set_delta(float Emin, float Emax, int iwin, float _Ebin)
{
   const float Erange = Emax - Emin;
   Ebin = _Ebin;
   if( Ebin<=0 ) Ebin = Erange/100.;
   if( Ebin>(Erange/10.) ) Ebin = Erange/10.;
   NBin = Erange/Ebin;
   Elo = Emin;
   Ehi = Elo + Ebin*NBin;
   iwindow = iwin;
}


void WL_Window::copy(const WL_Window& orig)
{
   Elo = orig.Elo;
   Ehi = orig.Ehi;
   Ebin = orig.Ebin;
   NBin = orig.NBin;
   iwindow = orig.iwindow;
}


int WL_Window::bin(float E) const  
{
   int ibin = static_cast<int>((E-Elo)/Ebin); 
   if(ibin<0) ibin=0;
   if(ibin>=NBin) ibin=NBin-1; // {  std::cout << "ibin=" << ibin << "/" << NBin-1 << " E=" << E <<std::endl; ibin=NBin-1; }
   return ibin;
}

#endif  // WLWalker_HPP_
