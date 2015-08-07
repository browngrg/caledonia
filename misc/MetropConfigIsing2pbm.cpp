// Converts a MC_Metropolis configuration line to a pbm picture


#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>


int main(int argc, char* argv[])
{
   int lineno=1;
   if(argc>1) sscanf(argv[1],"%d",&lineno);

   std::string buffer;
   int iline = 0;
   while( iline<lineno && !std::cin.eof() )
   {
       getline(std::cin,buffer);
       iline++;
   }
   if( iline==lineno && buffer.size()>0 )
   {
      std::istringstream fin(buffer);
      float kT,E,imcs;
      fin >> kT >> E >> imcs;
      std::vector<int> bmap;
      while( !fin.eof() ) 
      {
         int bit;
         fin >> bit;
         bmap.push_back(bit);
      }
      int L2 = bmap.size();
      int LY = static_cast<int>(std::sqrt(L2)+.01);
      int LX = L2/LY;
      std::cout << "P1 " << LX << " " << LY << std::endl;
      std::cout << "# kT=" << kT << std::endl;
      for(int i=0; i<L2; i++) bmap[i] = (bmap[i]+1)/2;
      int ii=0;
      for(int j=0; j<LY; j++)
      {
         std::cout << bmap[ii++];
         for(int i=1; i<LX; i++)
            std::cout << " " << bmap[ii++];
         std::cout << std::endl;
      }
   }

}
