#ifndef ITTM_HPP_
#define ITTM_HPP_

#include "MPI_Gang.hpp"

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm>

// A measurement object for the Infinite Temperature Transition Matrix method of calculating DOS
// Methods for Infinite-temperature Transition Matrix results
//
// Combined WL-Transition Matrix: RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
//
// C[I,J] is the number of proposed transition from I to J
// cbands  This is stored as CWIDTH diagonals, diagonal index = J-I-CHWIDTH
//         each diagonal is NBin long. To avoid knowing NBin, i = (J-I-CHWIDTH) + ibin*CWIDTH
//         so the diagonal index is the fastest increasing index. (3/3/2015 not sure "diagonal" is the right word)
// cob     The number of proposed transitions that could were out of the saved range
// 
// Combined WL-Transition Matrix: RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
// wleta comes from Eq. (24) S_WL[I] <- (1-wleta)*S_WL[I] + wleta*S_TM[I]
class ITTM
{
public:

   double Elo, Ehi, EBin;                   // The energy range and bin width
   int NBinE;

   enum { CHWIDTH = 5 };                    // Number of (super/sub)diagonals to keep in C matrix
   enum { CWIDTH = 2*CHWIDTH+1 };           // Number of diagonals to keep in C matrix (this is calculated)
   std::vector<long long> C;                // Number of proposed transitions from I to J
                                            // This is stored as CWIDTH diagonals, diagonal index = J-I-CHWIDTH
                                            // each diagonal is NBin long. To avoid knowing NBin, i = (J-I-CHWIDTH) + ibin*CWIDTH
                                            // so the diagonal index is the fastest increasing index.
   std::vector<long long> Cob;              // Number of proposed transitions out of saved diagonals
   std::vector<double>  Sittm;              // A local version of the Infinite-temperature transition matrix result
   std::vector<double>  Vittm;              // A local version of the Vittm error matrix

   bool WLTMMC_FULL;
   std::vector<long long> Cfull;            // Full matrix corresponding to data structure C (for debugging/analysis)
   long long CJhi, CJlo;                    // J outside window

   double    EStepSum;                      // Statistics of how big a MC step is
   long long EStepCnt;

   
   bool PD;                                 // Projective Dynamics data
   std::vector<long long> GrowCts;          // Binned growing and shrinking
   std::vector<long long> ShrkCts;

   MPI_Gang mp;

public:

   /// Constructor
   ITTM() { WLTMMC_FULL = false; PD=true; EBin = 0; }

   void init(double Emin, double Emax, double EBin);

   void write() const;

   void read();

   template<typename MCWalker>
   void add_sample(const MCWalker& walker, bool at_equilibrium=true);

   void add_sample(int J, int I);

   void get_global(ITTM& global) const;

   void calc_SITTM() { this->calc_SITTM(Sittm); }

   // Calculate S=lng from proposed moves (infinite temperature transition matrix)
   template<typename F>
   void calc_SITTM(std::vector<F>& SITTM) const;

   // Calculate S=lng from proposed moves (infinite temperature transition matrix)
   // Eq. (19) of RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
   template<typename L, typename F>
   static void calc_SITTM(const std::vector<L>& C, const std::vector<L>& Cob, int CWIDTH, std::vector<F>& SITTM);

   template<typename F>
   void calc_SITTM_Rayleigh(std::vector<F>& SITTM) const;

   void calc_VITTM() { this->calc_VITTM(Vittm); }

   template<typename F>
   void calc_VITTM(std::vector<F>& vttt) const;

   template<typename L>
   void calc_H(std::vector<L>& H) const;

   template<typename L>
   void calc_H(std::vector<L>& H, std::vector<L>& Hfull) const;

   template<typename L>
   static void calc_H(const std::vector<L>& C, const std::vector<L>& Cob, int CWIDTH, std::vector<L>& H);

};


void ITTM::init(double Emin, double Emax, double _EBin)
{
   EBin  = _EBin;
   NBinE = 0;
   if( EBin>0 ) 
   {
      NBinE = static_cast<int>((Emax-Emin)/EBin+0.5); 
      if( NBinE<1 ) NBinE = 1;
   }
   Elo = Emin;
   Ehi = Elo + NBinE*EBin;
   C.resize(CWIDTH*NBinE);
   for(int i=0; i<CWIDTH*NBinE; i++) C[i]=0;
   Cob.resize(NBinE);
   for(int i=0; i<NBinE; i++) Cob[i]=0;
   if( WLTMMC_FULL )
   {
      Cfull.resize(NBinE*NBinE);
      for(long i=0; i<NBinE*NBinE; i++) Cfull[i]=0;
   }
   CJlo = 0;
   CJhi = 0;
   EStepSum = 0;
   EStepCnt =0 ;
   if( PD )
   {
      GrowCts.resize(NBinE);
      for(int i=0 ;i<NBinE; i++) GrowCts[i]=0;
      ShrkCts.resize(NBinE);
      for(int i=0 ;i<NBinE; i++) ShrkCts[i]=0;
   }
} 


void ITTM::read()
{
   if( mp.iproc_pool==0 || C.size()==0 ) return;
   std::ifstream fin("ITTM_CBand.csv");
   if( fin && fin.is_open() )
   {
      for(int i=0; i<NBinE; i++)
      {
         long long ti;
         fin >> ti;
         for(int icol=0; icol<CWIDTH; icol++)
            fin >> C[ icol+i*CWIDTH ];
         fin >> Cob[i];
      }
      std::cout << __FILE__ << ":" << __LINE__ << " read C matrix" << std::endl;
      fin.close();
   }
   if( Cfull.size()==0  ) return;
   fin.open("ITTM_CFull.csv");
   if( fin && fin.is_open() )
   {
      long long ti;
      fin >> ti;
      for(int i=0; i<NBinE*NBinE; i++) fin >> Cfull[i];
      fin.close();
   }
}

void ITTM::get_global(ITTM& global) const
{
   global.Elo = Elo;
   global.Ehi = Ehi;
   global.EBin = EBin;
   global.NBinE = NBinE;
#  ifndef USE_MPI
   global.C = C;
   global.Cob = Cob;
   global.Cfull = Cfull;
   global.CJhi = CJhi;
   global.CJlo = CJlo;
   global.GrowCts = GrowCts;
   global.ShrkCts = ShrkCts;
   return;
#  else
   int NBuff = C.size();
   global.C.resize(NBuff);
   for(int i=0; i<NBuff; i++) global.C[i] = 0;
   MPI_Allreduce(&(C[0]),&(global.C[0]),NBuff,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
   NBuff = Cob.size();
   global.Cob.resize(NBuff);
   for(int i=0; i<NBuff; i++) global.Cob[i] = 0;
   MPI_Allreduce(&(Cob[0]),&(global.Cob[0]),NBuff,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
   global.EStepSum = 0;
   MPI_Allreduce(&(EStepSum),&(global.EStepSum),1,MPI_DOUBLE,MPI_SUM,mp.comm_pool);
   global.EStepCnt = 0;
   MPI_Allreduce(&(EStepCnt),&(global.EStepCnt),1,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
   if(PD)
   {
      // MPICH on linux fails without the print statements
      NBuff = GrowCts.size();
      global.GrowCts.resize(NBuff); for(int i=0; i<NBuff; i++) global.GrowCts[i]=0;
      global.ShrkCts.resize(NBuff); for(int i=0; i<NBuff; i++) global.ShrkCts[i]=0;
      //std::cout << __FILE__ << ":" << __LINE__ << "(" << mp.iproc_pool << ") NBuff=" << NBuff << std::endl;
      MPI_Allreduce(&(GrowCts[0]),&(global.GrowCts[0]),NBuff,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
      //std::cout << __FILE__ << ":" << __LINE__ << "(" << mp.iproc_pool << ") NBuff=" << NBuff << std::endl;
      MPI_Allreduce(&(ShrkCts[0]),&(global.ShrkCts[0]),NBuff,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
      //std::cout << __FILE__ << ":" << __LINE__ << "(" << mp.iproc_pool << ") NBuff=" << NBuff << std::endl;
   }
   if( Cfull.size()>0 )
   {
      NBuff = Cfull.size();
      global.Cfull.resize(NBuff);
      for(int i=0; i<NBuff; i++) global.Cfull[i] = 0;
      MPI_Allreduce(&(Cfull[0]),&(global.Cfull),NBuff,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
      std::vector<long long> blocal(2),bglobal(2);
      blocal[0] = CJlo;
      blocal[1] = CJhi;
      MPI_Allreduce(&(blocal[0]),&(bglobal[0]),2,MPI_LONG_LONG,MPI_SUM,mp.comm_pool);
      global.CJlo = bglobal[0];
      global.CJhi = bglobal[1];
   }
#  endif
}

void ITTM::write() const
{
   ITTM global;
   this->get_global(global);
   if( mp.iproc_pool == 0 )
   {
      std::ofstream fout("ITTM_CBand.csv");
      for(int i=0; i<NBinE; i++)
      {
         fout << i;
         for(int icol=0; icol<CWIDTH; icol++)
            fout << " " << global.C[ icol + i*CWIDTH ];
         fout << " " << global.Cob[ i ] << std::endl;
      }
   }
   if( mp.iproc_pool == 0 )
   {
      global.calc_SITTM();
      global.calc_VITTM();
      double smax = *(std::max_element(global.Sittm.begin(),global.Sittm.end()));
      for(int i=0; i<global.Sittm.size(); i++) global.Sittm[i] -= smax;
      double EStepMean = EStepSum/static_cast<double>(EStepCnt);
      double EStepFrac = EStepMean/EBin;
      bool EBinOK = (EStepFrac>1) && (EStepFrac<10);
      std::vector<double> SPD(NBinE);
      std::vector<double> htot(NBinE);
      std::vector<double> grow(NBinE);
      std::vector<double> shrk(NBinE);
      if( PD )
      {
         for(int ibin=0; ibin<NBinE; ibin++)
         {
            htot[ibin] = static_cast<double>(GrowCts[ibin]+ShrkCts[ibin]);
            grow[ibin] = static_cast<double>(GrowCts[ibin])/htot[ibin];
            shrk[ibin] = static_cast<double>(ShrkCts[ibin])/htot[ibin];
         }
         SPD[0] = 0;
         for(int ibin=1; ibin<NBinE; ibin++)
         {
            double delta_lng = 0;
            if( shrk[ibin]>0 )
               delta_lng = std::log(grow[ibin]/shrk[ibin]);
            SPD[ibin] = SPD[ibin-1] + EBin*delta_lng;
         }
         smax = *(std::max_element(SPD.begin(),SPD.end()));
         for(int i=0; i<SPD.size(); i++) SPD[i] -= smax;
      } 
      std::ofstream fout("ITTM_DOS.csv");
      fout << "# Infinite temperature transition matrix results" << std::endl;
      fout << "# Average Energy step = " << EStepMean << " = " << EStepFrac << " of EBin" << std::endl;
      fout << "# Column 1: Energy" << std::endl;
      fout << "# Column 2: Sittm, estimated density of states (ok=" << EBinOK << ")" << std::endl;
      fout << "# Column 3: Vittm, violation of the TTT property" << std::endl;
      if( PD ) 
      {
      fout << "# Column 4: Spd, estimated from infinite temperature Projective Dynamics" << std::endl;
      fout << "# Column 5: Grow, probability of E_proposed>E_old" << std::endl;
      fout << "# Column 6: Shrink, probability of E_proposed<=E_old" << std::endl;
      fout << "# Column 7: Cts, number of counts used to calculate Columns 4-7" << std::endl;
      }
      for(int ibin=0; ibin<NBinE; ibin++)
      {
         fout << Elo+(ibin+0.5)*EBin << " " << global.Sittm[ibin] << " " << global.Vittm[ibin];
         if( PD )
            fout << " " << SPD[ibin] << " " << grow[ibin] << " " << shrk[ibin] << " " << htot[ibin];
         fout << std::endl;
      }
   }
   if( Cfull.size()>0 && mp.iproc_pool == 0 )
   {
      // This branch not debugged. GPB 7/31/2015
      std::ofstream fout("ITTM_CFull.csv");
      fout << NBinE << std::endl;
      for(int i=0; i<NBinE*NBinE; i++) fout << " " << global.Cfull[i] << std::endl;
   }
}



template<typename MCWalker>
void ITTM::add_sample(const MCWalker& walker, bool at_equilibrium)
{
   double Enow = walker.now.E;
   double Eold = walker.old.E;
   EStepSum += std::fabs(Enow-Eold);
   EStepCnt++;
   if( Eold<Elo || Eold>Ehi ) return;
   int Jbin = static_cast<int>((Enow-Elo)/EBin);
   int Ibin = static_cast<int>((Eold-Elo)/EBin);
   add_sample(Jbin,Ibin);
   if( PD )
   {
      if( Enow>Eold )
      {
         GrowCts[Ibin]++;
      }
      else
      {
         ShrkCts[Ibin]++;
      }
   }
}


void ITTM::add_sample(int J, int I)
{
   if( I<0 || I>=NBinE ) return;
   if(J<0) { CJlo++; J=0; }
   if(J>=NBinE) { CJhi++; J=NBinE-1; }
   if(true && C.size()<(NBinE*CWIDTH) ) std::cout << __FILE__ << ":" << __LINE__ << " C.size()=" << C.size() << std::endl;
   int IJ_diff = J-I;
   int IJ_indx = IJ_diff + CHWIDTH;                        // Shift into storage scheme indexing
   if( IJ_indx>=0 && IJ_indx<CWIDTH) 
      C[ IJ_indx + I*CWIDTH ] += 1;          
   else
      Cob[ I ] += 1;
   if( false && Cfull.size()>0 )
   {
      int idx = I + J*NBinE;
      if( idx<0 ) std::cout << __FILE__ << ":" << __LINE__ << " index=" << idx << " I,J=" << I << "," << J << std::endl;
      if( idx>Cfull.size() ) std::cout << __FILE__ << ":" << __LINE__ << " index=" << idx << "/" << Cfull.size() << " I,J=" << I << "," << J << std::endl;
   }
   if( Cfull.size()>0 )
      Cfull[ I + J*NBinE ] += 1;
}



template<typename L>
void ITTM::calc_H(std::vector<L>& H) const
{
   // Find total transition for each macrostate
   H.resize(NBinE);
   for(int i=0; i<NBinE; i++) H[i] = Cob[i];    
   for(int i=0; i<NBinE; i++)
      for(int j=0; j<CWIDTH; j++) H[i] += C[j+i*CWIDTH];
   bool anyzero = false;
   for(int i=0; !anyzero && i<NBinE; i++) anyzero = (H[i]==0);
   if( anyzero ) 
   {
      bool allzero = true;
      for(int i=0; allzero && i<NBinE; i++) allzero = (H[i]==0);
      if(allzero)
      {
         std::cout << __FILE__ << ":" << __LINE__ << " all counting arrays are zero" << std::endl;
      }
      else
      {
         if(false) std::cout << __FILE__ << ":" << __LINE__ << " empty bin will cause division by zero "
                             << "NBinE=" << NBinE << " H[0]=" << H[0] << " H[NBinE-1]=" << H[NBinE-1] << std::endl;
      }
   }
}


template<typename L>
void ITTM::calc_H(std::vector<L>& H, std::vector<L>& Hfull) const
{
   // Find total transition for each macrostate
   H.resize(NBinE);
   for(int i=0; i<NBinE; i++) H[i] = Cob[i];    
   for(int i=0; i<NBinE; i++)
      for(int j=0; j<CWIDTH; j++) H[i] += C[j+i*CWIDTH];
   bool anyzero = false;
   for(int i=0; !anyzero && i<NBinE; i++) anyzero = (H[i]==0);
   if( anyzero ) 
   {
      bool allzero = true;
      for(int i=0; allzero && i<NBinE; i++) allzero = (H[i]==0);
      if(allzero)
      {
         std::cout << __FILE__ << ":" << __LINE__ << " all counting arrays are zero" << std::endl;
      }
      else
      {
         if(false) std::cout << __FILE__ << ":" << __LINE__ << " empty bin will cause division by zero "
                             << "NBinE=" << NBinE << " H[0]=" << H[0] << " H[NBinE-1]=" << H[NBinE-1] << std::endl;
      }
   }
   if( Cfull.size()>0 )
   {
      Hfull.resize(NBinE);
      for(int irow=0; irow<NBinE; irow++)
      {
         for(int jcol=0; jcol<NBinE; jcol++)
            Hfull[irow] += Cfull[irow+NBinE*jcol];
         if( Hfull[irow] != H[irow] )
            std::cout << __FILE__ << ":" << __LINE__ << " irow=" << irow << " Hfull, H = " << Hfull[irow] << ", " << H[irow] << std::endl;
      }
   }
}


// Calculate S=lng from proposed moves (infinite temperature transition matrix)
template<typename F>
void ITTM::calc_SITTM(std::vector<F>& SITTM) const
{
   // TODO: See if this short circuit works
   if(false)
   {
      calc_SITTM(this->C,this->Cob,CWIDTH,SITTM); 
      return;
   }
   // Eq. (19) of RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
   SITTM.resize(NBinE);
   if(false) std::cout << __FILE__ << ":" << __LINE__ << " NBinE=" << NBinE << std::endl;
   if(NBinE<3) return; 
   std::vector<long long> H;
   calc_H(H);
   // calculate entropy from infinite temperature transition matrix
   // M(I,J) -> M[ (J-I+CHWIDTH) + I*CWIDTH ]
   // shift index I+1 -> i : I -> i-1
   SITTM[0] = 1;
   for(int i=1; i<NBinE; i++)
   {
      //if( H[i]>0 && H[(i-1)]>0 && C[(CHWIDTH-1)+i*CWIDTH]>0 )
      if( C[(CHWIDTH+1)+(i-1)*CWIDTH]>0 && C[(CHWIDTH-1)+i*CWIDTH]>0 )
      {
         double Timi = static_cast<double>(C[(CHWIDTH+1)+(i-1)*CWIDTH])/static_cast<double>(H[(i-1)]);  // Superdiagonal
         double Tiim = static_cast<double>(C[(CHWIDTH-1)+i*CWIDTH])/static_cast<double>(H[i]);          // Subdiagonal
         SITTM[i] = SITTM[i-1] + std::log(Timi/Tiim);    
      }
      else
      {
         //std::cout << __FILE__ << ":" << __LINE__ << " C==0 ibin=" << i << std::endl;
         SITTM[i] = 0;
      }
   }
}


template<typename F>
void ITTM::calc_VITTM(std::vector<F>& vttt) const
{
   // Eq. (16) of RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28, 715-726 (2007)
   // vttt = Violation of the total transition matrix (TTT)
   // This is called detailed balance violation
   // Eq. (11) of J.-S. Wang, L.W. Lee, Computer Physics Communications 127, 131-136 (2000)
   // which references J.-S. Wang, Eur. Phys. J. B 8, 287 (1999)
   // This doesn't work if bins are large and C matrix is tridiagonal, see Fig. 3 of Ghulghazaryan
   // The I is diagonal, I+1 the adjacent diagonals, and I+2 is the next diagonal.
   // Can still calculate S(I) when C is tridiagonal, but won't have error estimate.
   if(false) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   if(vttt.size()<NBinE) vttt.resize(NBinE);
   if( NBinE<3 ) return;
   std::vector<long long> H;
   calc_H(H);
   for(int i=0; i<(NBinE-3); i++)
   {
#        if 0
      double T01 = static_cast<double>(C[(CHWIDTH+1)+(i+0)*CWIDTH])/static_cast<double>(H[i]);       // Tinf(I,I+1)
      double T12 = static_cast<double>(C[(CHWIDTH+2)+(i+1)*CWIDTH])/static_cast<double>(H[i+1]);     // Tinf(I+1,I+2)
      double T20 = static_cast<double>(C[(CHWIDTH+0)+(i+2)*CWIDTH])/static_cast<double>(H[i+2]);     // Tinf(I+2,I)
      double T10 = static_cast<double>(C[(CHWIDTH+0)+(i+1)*CWIDTH])/static_cast<double>(H[i+1]);     // Tinf(I,I+1)
      double T21 = static_cast<double>(C[(CHWIDTH+1)+(i+2)*CWIDTH])/static_cast<double>(H[i+2]);     // Tinf(I+1,I+2)
      double T02 = static_cast<double>(C[(CHWIDTH+2)+(i+0)*CWIDTH])/static_cast<double>(H[i]);       // Tinf(I+2,I)
      vttt[i] = 1 - (T01*T12*T20)/(T02*T21*T10);
#        else
      // This version avoids division by zero
      double T01 = static_cast<double>(C[(CHWIDTH+1)+(i+0)*CWIDTH]); // /static_cast<double>(H[i]);       // Tinf(I,I+1)
      double T02 = static_cast<double>(C[(CHWIDTH+2)+(i+0)*CWIDTH]); // /static_cast<double>(H[i]);       // Tinf(I+2,I)
      double T12 = static_cast<double>(C[(CHWIDTH+2)+(i+1)*CWIDTH]); // /static_cast<double>(H[i+1]);     // Tinf(I+1,I+2)
      double T20 = static_cast<double>(C[(CHWIDTH+0)+(i+2)*CWIDTH]); // /static_cast<double>(H[i+2]);     // Tinf(I+2,I)
      double T10 = static_cast<double>(C[(CHWIDTH+0)+(i+1)*CWIDTH]); // /static_cast<double>(H[i+1]);     // Tinf(I,I+1)
      double T21 = static_cast<double>(C[(CHWIDTH+1)+(i+2)*CWIDTH]); // /static_cast<double>(H[i+2]);     // Tinf(I+1,I+2)
      double numer = T01*T12*T20;
      double denom = T02*T21*T10;
      if( denom==0 )
         vttt[i] = 0;
      else
         vttt[i] = 1 - numer/denom;
#        endif
   }
   for(int i=(NBinE-3); i<NBinE; i++) vttt[i] = 0;
   if(false) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
}


template<typename L>
void ITTM::calc_H(const std::vector<L>& C, const std::vector<L>& Cob, int CWIDTH, std::vector<L>& H)
{
   // Find total transition for each macrostate
   int nbin = Cob.size();
   H.resize(nbin);
   for(int i=0; i<nbin; i++) H[i] = Cob[i];    
   for(int i=0; i<nbin; i++)
      for(int j=0; j<CWIDTH; j++) H[i] += C[j+i*CWIDTH];
   bool anyzero = false;
   for(int i=0; !anyzero && i<nbin; i++) anyzero = (H[i]==0);
   if( anyzero ) 
   {
      bool allzero = true;
      for(int i=0; allzero && i<nbin; i++) allzero = (H[i]==0);
      if(allzero)
      {
         std::cout << __FILE__ << ":" << __LINE__ << " all counting arrays are zero" << std::endl;
      }
      else
      {
         if(false) std::cout << __FILE__ << ":" << __LINE__ << " empty bin will cause division by zero "
                             << "nbin=" << nbin << " H[0]=" << H[0] << " H[nbin-1]=" << H[nbin-1] << std::endl;
      }
   }
}


// Calculate S=lng from proposed moves (infinite temperature transition matrix)
// Eq. (19) of RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007)
template<typename L, typename F>
void ITTM::calc_SITTM(const std::vector<L>& C, const std::vector<L>& Cob, int CWIDTH, std::vector<F>& SITTM)
{
   int CHWIDTH = (CWIDTH-1)/2;
   int nbin = Cob.size();
   SITTM.resize(nbin);
   if(false) std::cout << __FILE__ << ":" << __LINE__ << " nbin=" << nbin << std::endl;
   if(nbin<3) return; 
   std::vector<long long> H;
   calc_H(C,Cob,CWIDTH,H);
   // calculate entropy from infinite temperature transition matrix
   // M(I,J) -> M[ (J-I+CHWIDTH) + I*CWIDTH ]
   // shift index I+1 -> i : I -> i-1
   SITTM[0] = 1;
   for(int i=1; i<nbin; i++)
   {
      //if( H[i]>0 && H[(i-1)]>0 && C[(CHWIDTH-1)+i*CWIDTH]>0 )
      if( C[(CHWIDTH+1)+(i-1)*CWIDTH]>0 && C[(CHWIDTH-1)+i*CWIDTH]>0 )
      {
         double Timi = static_cast<double>(C[(CHWIDTH+1)+(i-1)*CWIDTH])/static_cast<double>(H[(i-1)]);  // Superdiagonal
         double Tiim = static_cast<double>(C[(CHWIDTH-1)+i*CWIDTH])/static_cast<double>(H[i]);          // Subdiagonal
         SITTM[i] = SITTM[i-1] + std::log(Timi/Tiim);    
      }
      else
      {
         SITTM[i] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Legacy Codes
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef WL_WALKER_RAYLEIGH_HPP
#define WL_WALKER_RAYLEIGH_HPP

#include<cmath>

// DGEBAL - balance a general real matrix A
#ifdef USE_BLAS
extern "C" void dgeev_(char* jobvl, char* jobvr, int* N, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int *lwork, int* info);
extern "C" void dgemm_(char* transa, char* transb ,int* M, int* N, int* K, const double* alpha, const double* A, int* lda, const double* B, int* ldb, const double* beta, double* C, int* ldc);
extern "C" void dgemv_(char* trans, int* M, int* N, const double* alpha, const double* A, int* lda, const double* X, int* incx, const double* beta, double* y, int* incy);
extern "C" void dgetrf_(int* M, int* N, double* A, int* lda, int* ipivot, int* info);
extern "C" void dgetri_(int* N, double* A, int* lda, int* ipivot, double* work, int* lwork, int* info);
#endif

// Returns the the Rayleigh quotient of x with respect to A, r(x) = xT*A*x/(xT*x)
double Rayleigh(const std::vector<double>& A, const std::vector<double>& x)
{
#  ifndef USE_BLAS
   std::cout << __FILE__ << ":" << __LINE__ << " member function Rayleigh(A,x) requires -DUSE_BLAS" << std::endl;
   return 1;
#  endif
   // This implementation assumes xT*x = 1
   int N = x.size();
   char trans = 'N';
   int incr = 1;
   double alpha = 1;
   double beta = 0;
   std::vector<double> Ax(N,0);
   dgemv_(&trans,&N,&N,&alpha,&(A[0]),&N,&(x[0]),&incr,&beta,&(Ax[0]),&incr);
   double xTAx = 0;
   double xTx = 0;
   for(int i=0; i<N; i++) 
   {
      xTAx += x[i]*Ax[i];
      xTx  += x[i]*x[i];
   }
   return xTAx/xTx;
}

double RayleighQuotientIteration(const std::vector<double>& mat, std::vector<double>& evec, int niter=20)
{
#  ifndef USE_BLAS
   std::cout << __FILE__ << ":" << __LINE__ << " member function RayleighQuotientIteration requires -DUSE_BLAS" << std::endl;
   return 0;
#  endif
   // Converge to the eigenvector using the Rayleigh quotient method
   // Returns the associated eigenvalue
   // This is known to have problems for non-symmetric matrice
   // Maysum Panju, "Iterativec Methods for Computing Eigenvalues and Eigenvectors,"
   // The Waterloo Mathematics Review 1, 9-18 (2011).
   int N = evec.size();
   double vmag = 0; for(int i=0; i<N; i++) vmag += evec[i]*evec[i];      // Normalize the eigenvector
   vmag = sqrt(vmag);
   for(int i=0; i<N; i++) evec[i] /= vmag;
   double lambda = Rayleigh(mat,evec);
   //cout << "# converge init lambda=" << lambda << std::endl;
   for(int k=0; k<niter; k++)
   {
      // Solve (A-lambda*I)*w = v
      std::vector<double> matinv = mat;
      for(int i=0; i<N; i++) matinv[i+i*N] -= lambda;                                    // mat = A - lambda*I
      std::vector<int> ipivot(N);
      int info;
      dgetrf_(&N,&N,&(matinv[0]),&N,&(ipivot[0]),&info);                                 // LU factorization of matrix A=P*L*U
      std::vector<double> work(10*N);
      int lwork = work.size();
      dgetri_(&N,&(matinv[0]),&N,&(ipivot[0]),&(work[0]),&lwork,&info);                  // Computes inverse matrix of LU factorized matrix
      char trans = 'N';
      int incr = 1;
      double alpha = 1;
      double beta = 0;
      std::vector<double> w(N);
      dgemv_(&trans,&N,&N,&alpha,&(matinv[0]),&N,&(evec[0]),&incr,&beta,&(w[0]),&incr);  // w = mat^{-1}*v
      double wmag = 0; for(int i=0; i<N; i++) wmag += w[i]*w[i];
      wmag = sqrt(wmag);
      for(int i=0; i<N; i++) evec[i] = w[i]/wmag;
      lambda = Rayleigh(mat,evec);
   }
   return lambda;
}

#endif
template<typename F>
void ITTM::calc_SITTM_Rayleigh(std::vector<F>& _SITTM) const
{
   // This algorithm is copied after AnalyzeFullC::main() 
   if( Cfull.size()==0 ) return;
   int nbin = NBinE;
   // look for empty rows and columns
   int minbin = 0;
   long long cts = 0;
   while( minbin<nbin && cts==0 )
   {
      for(int jcol=minbin; jcol<nbin; jcol++)
         cts += Cfull[minbin+nbin*jcol];
      minbin++;
   }
   if(cts>0) minbin--;
   cts = 0;
   while( minbin<nbin && cts==0 )
   {
      for(int irow=minbin; irow<nbin; irow++)
         cts += Cfull[irow+nbin*minbin];
      minbin++;
   }
   if(cts>0) minbin--;
   // Construct infinite-temperature transition matrix (possibly smaller than Cfull)
   // Transpose so can work with right eigenvector
   // Normalize rows
   int nold = nbin;
   nbin = nbin - minbin;
   std::vector<double> Tinf(nbin*nbin);
   std::vector<long long> Hsum(nbin);
   for(int irow=0; irow<nbin; irow++)
   {
      Hsum[irow] = 0;
      for(int jcol=0; jcol<nbin; jcol++)
         Hsum[irow] += Cfull[(jcol+minbin)+nold*(irow+minbin)];
      for(int jcol=0; jcol<nbin; jcol++)
         Tinf[irow+nbin*jcol] = static_cast<double>(Cfull[(jcol+minbin)+nold*(irow+minbin)])/static_cast<double>(Hsum[irow]);
   }
   // Calculate beta from Eq. (3) of PMC de Oliveira, TJP Penna, HJ Herrmann Eur Phys J B 1, 205 (1998)
   std::vector<double> beta(nbin);   // actually deltaE*beta(E)
   std::vector<double> lng3(nbin);    // from all up/down
   lng3[0] = 0;
   for(int irow=0; irow<(nbin-1); irow++)
   {
      beta[irow] =  -log(Tinf[irow+nbin*(irow+1)]/Tinf[(irow+1)+nbin*irow]); // minus sign cuz transpose of Tinf in paper
      if( beta[irow]!=beta[irow] ) beta[irow]=0;
      lng3[irow+1] = lng3[irow] + beta[irow];
   }
   beta[nbin-1] = beta[nbin-2];
   if( lng3[0]!=0 ) for(int irow=nbin-1; irow>=0; irow--) lng3[irow] -= lng3[0];   // fix lng so it starts at zero and is positive
   if( lng3[nbin-1]<0 ) for(int irow=0; irow<nbin; irow++) lng3[irow] *= -1;
   // It gives different results from Eq. (19) of RG Ghulghazaryan, S. Hayryan, CK Hu, J Comp Chem 28,  715-726 (2007) (calc_SITTM)
   //    double Timi = static_cast<double>(C[(CHWIDTH+1)+(i-1)*CWIDTH])/static_cast<double>(H[(i-1)]);  // Superdiagonal
   //    double Tiim = static_cast<double>(C[(CHWIDTH-1)+i*CWIDTH])/static_cast<double>(H[i]);          // Subdiagonal
   //    SITTM[i] = SITTM[i-1] + std::log(Timi/Tiim);    
   std::vector<double> lng19(nbin);
   calc_SITTM(lng19);
   if( lng19[0]!=0 ) for(int irow=nbin-1; irow>=0; irow--) lng19[irow] -= lng19[0];
   if( lng19[nbin-1]<0 ) for(int irow=0; irow<nbin; irow++) lng19[irow] *= -1;
   // Improve estimate by finding eigenvector of transition matrix
   std::vector<double> lng = lng3;
   double gmax = std::max(lng[0],lng[nbin-1]);
   for(int irow=0; irow<nbin; irow++) lng[irow] = exp(lng[irow]-gmax);
   double lambda = RayleighQuotientIteration(Tinf,lng);
   for(int irow=0; irow<nbin; irow++) lng[irow] = log(fabs(lng[irow]));
   if( lng[0]!=0 ) for(int irow=nbin-1; irow>=0; irow--) lng[irow] -= lng[0];
   if( lng[nbin-1]<0 ) for(int irow=0; irow<nbin; irow++) lng[irow] *= -1;
   // TODO: MPI_Reduce(SUM), only print if iproc 0;
   if( true )
   {
      char buff[500];
      sprintf(buff,"WalkerDOS-%02d.csv",mp.iproc_pool);
      std::ofstream fout(buff);
      fout << "# Data from Analysis of Full C Matrix" << std::endl;;
      fout << "# orig nbin = " << nold <<" analyze nbin=" << nbin << std::endl;;
      fout << "# lambda = " << lambda << std::endl;;
      fout << "# Column 1: ibin" << std::endl;;
      fout << "# Column 2: iterative eigenvector estimate" << std::endl;;
      fout << "# Column 3: Broad-histogram estimate of lng(E) Tinf Eq. 3" << std::endl;;
      fout << "# Column 4: ITTM estiamte of lng(E) Tinf Eq. (19)" << std::endl;;
      fout << "# Column 5: deltaE*beta" << std::endl;;
      fout << "# Column 6: H total counts in bin" << std::endl;;
      fout << std::setprecision(16);
      for(int irow=0; irow<nbin; irow++)
      {
         fout << irow << " " << lng[irow] << " " << lng3[irow] << " " << lng19[irow] << " " << beta[irow] << " " << Hsum[irow] << std::endl;
      }
   }
   _SITTM = lng;
}

#endif  //  ITTM_HPP_
