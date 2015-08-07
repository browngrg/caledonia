// Converts EMMap.csv to a greyscale picture


#include<iostream>
#include<vector>
#include<string>


int main(int argc, char* argv[])
{


   std::vector<float> array;
   int NX = 0;
   int NY = 0;

   float lastx;
   while( !std::cin.eof() )
   {
      std::string buffer;
      getline(std::cin,buffer);
      if( buffer.size()==0 ) continue;
      if( buffer[0]=='#' ) continue;
      float E,M,norm;
      long raw;
      sscanf(buffer.c_str(),"%f %f %f %ld",&E,&M,&norm,&raw);
      if(array.size()==0) lastx = E;
      if( E!=lastx ) NY++;
      array.push_back(norm);
      lastx = E;
   }
   NY++;

   NX = array.size()/NY;
   std::cout << "P2 " << NX << " " << NY << " 1000" << std::endl;
   int iel = 0;
   for(int iy=0; iy<NY; iy++)
   {
      std::cout << static_cast<int>(1000*array[iel++]);
      for(int ix=1; ix<NX; ix++)
         std::cout << " " << static_cast<int>(1000*array[iel++]);
      std::cout << std::endl;
   }

}
