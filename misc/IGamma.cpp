

#include<iostream>
#include<cmath>

// Return integer index of wlgamma assuming wlgamma_i+1 = wlgamma_i/2
int igamma(double wlgamma)
{
   int p2 = static_cast<int>(-std::log(wlgamma)/std::log(2.)+0.01);
   return p2;
}


int main(int argc, char* argv[])
{

   std::cout << "# Test of igamma(wlgamma) function" << std::endl;
   std::cout << "# Column 1: iteration" << std::endl;
   std::cout << "# Column 2: wlgamma" << std::endl;
   std::cout << "# Column 3: function return value" << std::endl;
   double wlgamma = 1;
   for(int i=0; i<20; i++)
   {
      int ig= igamma(wlgamma);
      std::cout << i << " " << wlgamma << " " << ig << std::endl;
      wlgamma /= 2.;
   }
}

