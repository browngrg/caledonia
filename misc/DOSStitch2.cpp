// Convert as DOS file to Thermodynamic data
// This version also fits polynomials to the data

// Does not assume same bin size for different files

#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
#include<string>
#include<limits>
#include<algorithm>
#include<unistd.h>
#include<stdio.h>


// Reads even a simple x,y file (2010-12-17 GPB)
// Added selecting index of columns  (2014-06-23 GPB)
template<typename T>
void loadtxt_xy(std::string fname_in, int ix, int iy, std::vector<T>& x, std::vector<T>& y)
{
   const bool verbose = false;
   int imax = std::max(ix,iy);
   if( access(fname_in.c_str(),F_OK) == -1 )
   {
      std::cout << "\"" << fname_in << "\" does not exist" << std::endl;
      return;
   }
   std::ifstream input(fname_in.c_str());
   if(verbose) std::cout << "Processing " << fname_in << " " << static_cast<bool>(input) << " " << input.is_open() << " " << !input.eof() << std::endl;
   std::string buffer;
   const std::string separator="\t,; ";
   std::vector<T> number;
   while( input.is_open() && !input.eof() )
   {
      std::getline(input,buffer);
      // remove comments
      size_t hash = buffer.find_first_of("#");
      if( hash!=std::string::npos )
         buffer = buffer.substr(0,hash);
      // get numbers by column
      number.resize(0);
      size_t l = buffer.find_first_not_of(separator);
      while( l<buffer.size() )
      {
         size_t r = buffer.find_first_of(separator,l);
         if( r>l )
         {
            double v = 0;
            std::istringstream strm(buffer.substr(l,r-l));
            strm >> v;
            number.push_back(v);
         }
         l = buffer.find_first_not_of(separator,r);
      }
      // use the results
      if( number.size()>imax )
      {
         x.push_back( number[ix] );
         y.push_back( number[iy] );
      }
   }
   if(verbose)
   {
      std::cout << "     found " << x.size() << " data lines for columns " << ix << "," << iy << " imax=" << imax << std::endl;
      std::cout << "     last line = \"" << buffer << "\"" << std::endl;
   }
}



struct WinRec 
{ 
public: 
   std::vector<double> E; 
   std::vector<double> lng; 
   std::vector<double> slope;
public:
   WinRec(const std::vector<double>& Eval, const std::vector<double>& lngval) { E=Eval; lng=lngval; find_slope(); }
   void find_slope()
   {
      int npts = E.size();
      slope.resize(npts);
      for(int i=1; i<npts; i++)
         slope[i] = (lng[i]-lng[i-1])/(E[i]-E[i-1]);
      slope[0] = slope[1];
   }
   void smooth_slope(int npt=3)
   {
      std::vector<double> smooth(slope.size());
      for(int i=npt; i<(slope.size()-npt); i++)
      {
         smooth[i] = 0;
         for(int j=i-npt; j<i+npt; j++)
            smooth[i] += slope[j];
         smooth[i] /= static_cast<double>(2*npt+1);
      }
      for(int i=npt; i<(slope.size()-npt); i++)
         slope[i] = smooth[i]; 
   }
   int find_energy(double Eval, int istart=0)
   {
      while( istart<E.size() && E[istart]<Eval ) istart++; 
      return istart;
   }
};


int main(int argc, char* argv[])
{

   using namespace std;

   if(true)
   {
      cout << argv[0] << " DOSFile1 DOSFile2" << endl;
      cout << "Stitches to DOS Files together" << std::endl;
   }
   if(argc<3) return 0;

   // Get valid files
   vector<string> fname_in;
   for(int i=1; i<argc; i++)
   {
      if( access(argv[i],F_OK) == -1 )
      {
         cout << "\"" << argv[i] << "\" does not exist" << endl;
      }
      else
      {
         fname_in.push_back( string(argv[i]) );
      }
   }

   vector<WinRec> window;

   int ECol = 0;
   int LngCol = 6;
   vector<double> E,lng;
   int NFILE = fname_in.size();
   cout << "NFILE=" << NFILE << endl;
   for(int i=0; i<NFILE; i++)
   {
      E.resize(0);
      lng.resize(0);
      loadtxt_xy(fname_in[i],ECol,LngCol,E,lng);
      window.push_back( WinRec(E,lng) );
   }

   int NWin = 5;
   for(int i=0; i<NFILE; i++)
      window[i].smooth_slope(NWin);

   ofstream fout("DOSStitchWin.csv");
   fout << "# Read Input with ECol=" << ECol+1 << " LngCol=" << LngCol+1 << endl;
   fout << "# Column 1: Eleft" << endl;
   fout << "# Column 2: slope" << endl;
   fout << "# Column 3: lng" << endl;
   fout << "# Column 4: Erght" << endl;
   fout << "# Column 5: slope" << endl;
   fout << "# Column 6: lng" << endl;
   vector<double> Ejoin(NFILE);  Ejoin[0] = window[0].E.front();
   vector<double> shift(NFILE);  shift[0] = 0;
   int ileft_start = 0;
   int iwr = 0;
   int iwl = 0;
   vector<double> merge_E;
   vector<double> merge_lng; 
   for(int i=0; i<(NFILE-1); i++)
   {
      iwr = i+1;
      iwl = i;
      double meanhi = 0;
      for(int j=0; j<window[iwr].E.size(); j++) meanhi += window[iwr].E[j]; 
      meanhi /= static_cast<double>(window[iwr].E.size());
      double meanlo =  0;
      for(int j=0; j<window[iwl].E.size(); j++) meanlo += window[iwl].E[j]; 
      meanlo /= static_cast<double>(window[iwl].E.size());
      if( meanhi<meanlo ) std::swap(iwr,iwl);
      double Elo = window[iwr].E.front();  // lowest energy of right window
      double Ehi = window[iwl].E.back();   // highest energy of left window
      if( Elo>Ehi )
      {
         cout << "Windows do not overlap i=" << i << endl;
         return 1;
      }
      cout << "i=" << i << " Elo=" << Elo << " Ehi=" << Ehi << endl;
      int ileft   = window[iwl].find_energy(Elo);
      for(int ipt=ileft_start; ipt<ileft; ipt++)
      {
         merge_E.push_back( window[iwl].E[ipt] );
         merge_lng.push_back( window[iwl].lng[ipt] + shift[i] );
      }
      int irght   = window[iwr].find_energy(window[iwl].E[ileft]);
      double diff = window[iwr].slope[irght]- window[iwl].slope[ileft];  
      int foo = 2*(diff>0) - 1;                                        // -1 or 1
      while( (foo*diff)>0 && ileft<window[iwl].E.size() )
      {
         merge_E.push_back( window[iwl].E[ileft] );
         merge_lng.push_back( window[iwl].lng[ileft] + shift[i] );
         irght = window[iwr].find_energy( window[iwl].E[ileft] );
         diff = window[iwl].slope[ileft] - window[iwr].slope[irght];
         ileft++;
      }
      ileft_start = irght;
      Ejoin[i+1] = window[iwl].E[ileft];
      shift[i+1] =  ( window[iwl].lng[ileft] + shift[i] ) - window[iwr].lng[irght];
      cout << "i=" << i << " iright_start=" << ileft_start << " Ejoin=" << Ejoin[i+1] << " shift=" << shift[i+1] << endl;
      cout << "  shift[i]=" << shift[i] << " lng_left=" << window[iwl].lng[ileft] << " lng_right=" << window[iwr].lng[irght] << endl;
      cout << "  Eleft=" << window[iwl].E[ileft] << " Eright=" << window[iwr].E[irght] << endl;
      // This corrupts ileft and irght
      while( ileft<window[iwl].E.size() && window[iwl].E[ileft]<Ehi )
      {
         fout << window[iwl].E[ileft] << " " << window[iwl].slope[ileft] << " " << window[iwl].lng[ileft] << " "
              << window[iwr].E[irght] << " " << window[iwr].slope[irght] << " " << window[iwr].lng[irght] << endl;
         ileft++; irght++;
      }
   }
   for(int ipt=ileft_start; ipt<window[iwr].E.size(); ipt++)
   {
      merge_E.push_back( window[iwr].E[ipt] );
      merge_lng.push_back( window[iwr].lng[ipt] + shift[iwr] );
   }
   fout.close();

   fout.open("DOSStitchResult.csv");
   fout << "# Column 1: E" << endl;
   fout << "# Column 2: lng" << endl;
   fout << "# Column 3: iwin" << endl;
   fout << "# Column 4: shift" << endl;
   int iwin = 0;
   for(int ipt=0; ipt<merge_E.size(); ipt++)
   {
      if( merge_E[ipt] > Ejoin[iwin+1] ) iwin++;
      fout << merge_E[ipt] << " " << merge_lng[ipt] << " " << iwin << " " << shift[iwin] << endl;
   }



}
