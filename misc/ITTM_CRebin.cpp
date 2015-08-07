
#include"ITTM.hpp"

#include<iostream>
#include<vector>
#include<sstream>
#include<string>


int main(int argc, char* argv[])
{
   std::vector<long long> cbands;
   std::vector<long long> cob;
   int ncol;   
   read_C(std::cin,cbands,cob,ncol);
   int nrow = cob.size();


#  if 1
   // Only keep every nth bin
   int rFac = 2;
   int odd_row = 0;
   if(argc>1) 
   {
      sscanf(argv[1],"%d",&rFac);
   }
   else
   {
      int first_row = 0;
      int idiag = (ncol-1)/2;
      while( first_row<nrow && cbands[idiag+first_row*ncol]==0 ) first_row++;
      int next_row = odd_row+1;
      while( next_row<nrow && cbands[idiag+next_row*ncol]==0 ) next_row++;
      rFac = next_row - first_row;
      odd_row = ( next_row - first_row ) % rFac;
   }
   int nrow2 = nrow/rFac;
   int nhcol = (ncol-1)/2;
   int odd_col = nhcol % rFac;
   int nhcol2 = nhcol/rFac;
   int ncol2 = 2*nhcol2+1;
#  else
   // Drop the even bins
   int odd_row = 0;
   int nrow2 = nrow/2;
   int nhcol = (ncol-1)/2;
   int odd_col = nhcol%2;
   int nhcol2 = nhcol/2;
   int ncol2 = 2*nhcol2+1;
#  endif
   if(false) std::cout << "rFac=" << rFac << " orow=" << odd_row << "ncol = " << ncol << " ncol2=" << ncol2 << std::endl;
   std::vector<long long> cbands2( nrow2*ncol2 );
   std::vector<long long> cob2( nrow2 );
   for(int irow=0; irow<nrow2; irow++)
   {
      cob2[irow] = cob[rFac*irow+odd_row];
      for(int icol=0; icol<ncol2; icol++)
         cbands2[icol+irow*ncol2] = cbands[(rFac*icol+odd_col)+ncol*(rFac*irow+odd_row)];
   }

   write_C(std::cout,cbands2,cob2,ncol2);

   return 0;
}
