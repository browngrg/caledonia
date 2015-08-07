
#include<iostream>
#include<sstream>
#include<string>


// Get X,Y output from the input file
void grabcol(int xcol, int ycol, std::istream& fin, std::ostream& fout)
{
   xcol--;
   ycol--;
   std::string buff;
   double vals[50];
   while( !fin.eof() )
   {
      std::getline(fin,buff);
      if( buff.size()==0 ) continue;
      if( buff[0]=='#' ) continue;
      std::istringstream sin(buff);
      int icol=0;
      while( !sin.eof() ) sin >> vals[icol++];
      fout << vals[xcol] << " " << vals[ycol] << std::endl;
   }
}


int main(int argc, char* argv[])
{
   grabcol(1,7,std::cin,std::cout);
}
