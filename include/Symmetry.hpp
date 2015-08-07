#ifndef SYMMETRY_HPP_
#define SYMMETRY_HPP_

#include"Vec.hpp"

#include<string>
#include<vector>
#include<algorithm>
#include<iostream>
#include<sstream>

// References:
// Hans Wondratschek, "Matrices, Mappings and Crystallographic Symmetry"
// Bilboa Crystallographic Server:  http://www.cryst.ehu.es
// NRL pages: http://cst-www.nrl.navy.mil/lattice

namespace Sym
{

   // ==========================================================================
   // ===== General information about Bravais lattices and crystal systems
   // ==========================================================================

   /// Adjectives describing  objects have context-sensitive meanings
   enum Adjective
   {
      INVALID        = 9,
      TRICLINIC      = 0,
      MONOCLINIC     = 1,
      ORTHORHOMBIC   = 2,
      TETRAGONAL     = 3,
      TRIGONAL       = 4,
      HEXAGONAL      = 5,
      CUBIC          = 6,
      RHOMBOHEDRAL   = 7,
      ORTHOHEXAGONAL = 8
   };

   /// Total number of crystal systems
   enum { NCrystalSystem = 7 };

   /// Compile-time name resolution
   template<int i> std::string CrystalSystemString()            { return "invalid";     }
   template<> std::string CrystalSystemString<TRICLINIC>()      { return "triclinic";   }
   template<> std::string CrystalSystemString<MONOCLINIC>()     { return "monoclinic";  }
   template<> std::string CrystalSystemString<ORTHORHOMBIC>()   { return "orthorhombic";}
   template<> std::string CrystalSystemString<TETRAGONAL>()     { return "tetragonal";  }
   template<> std::string CrystalSystemString<TRIGONAL>()       { return "trigonal";    }
   template<> std::string CrystalSystemString<HEXAGONAL>()      { return "hexagonal";   }
   template<> std::string CrystalSystemString<ORTHOHEXAGONAL>() { return "hexagonal";   }
   template<> std::string CrystalSystemString<CUBIC>()          { return "cubic";       }

   /// Run-time name resolution
   std::string CrystalSystemString(Adjective i)
   {
      switch(i)
      {
         case TRICLINIC:      return "triclinic";    break;
         case MONOCLINIC:     return "monoclinic";   break;
         case ORTHORHOMBIC:   return "orthorhombic"; break;
         case TETRAGONAL:     return "tetragonal";   break;
         case TRIGONAL:       return "trigonal";     break;
         case HEXAGONAL:
         case ORTHOHEXAGONAL: return "hexagonal";    break;
         case CUBIC:          return "cubic";        break;
         default:             return "invalid";
      }
   }

   // Convert string to value
   Adjective ParseCrystalSytem(std::string adj)
   {
      Adjective ival = INVALID;
      for(int i=0; i<9 && ival==INVALID; i++)
         if( adj.compare(CrystalSystemString(static_cast<Adjective>(i)))==0 ) ival=static_cast<Adjective>(i);
      return ival;
   }

   template<int i> std::string LatticeSystemString()            { return "invalid";      }
   template<> std::string LatticeSystemString<TRICLINIC>()      { return "triclinic";    }
   template<> std::string LatticeSystemString<MONOCLINIC>()     { return "monoclinic";   }
   template<> std::string LatticeSystemString<ORTHORHOMBIC>()   { return "orthorhombic"; }
   template<> std::string LatticeSystemString<TETRAGONAL>()     { return "tetragonal";   }
   template<> std::string LatticeSystemString<RHOMBOHEDRAL>()   { return "rhombohedral"; }
   template<> std::string LatticeSystemString<HEXAGONAL>()      { return "hexagonal";    }
   template<> std::string LatticeSystemString<ORTHOHEXAGONAL>() { return "hexagonal";    }
   template<> std::string LatticeSystemString<CUBIC>()          { return "cubic";        }
  
   std::string LatticeSystemString(Adjective i)
   {
      switch(i)
      {
         case TRICLINIC:      return "triclinic";    break;   // alpha != beta != gamma      a, b, c
         case MONOCLINIC:     return "monoclinic";   break;   // alpha==gamma==90 != beta    a, b, c
         case ORTHORHOMBIC:   return "orthorhombic"; break;   // alpha==beta==gamma==90      a, b, c
         case TETRAGONAL:     return "tetragonal";   break;   // alpha==beta==gamma==90      a, a, c
         case RHOMBOHEDRAL:   return "rhombohedral"; break;   // alpha==beta==90  != gamma   a, a, a
         case HEXAGONAL:                                      // alpha==beta==90 gamma==120  a, a, c
         case ORTHOHEXAGONAL: return "hexagonal";    break;   // alpha==beta==gamma==90
         case CUBIC:          return "cubic";        break;   // alpha==beta==gamma==90      a, a, a
         default:             return "invalid";
      }
   }

   // Convert string to value
   Adjective ParseLatticeSystem(std::string adj)
   {
      Adjective ival = INVALID;
      for(int i=0; i<9 && ival==INVALID; i++)
         if( adj.compare(LatticeSystemString(static_cast<Adjective>(i)))==0 ) ival=static_cast<Adjective>(i);
      return ival;
   }

   // Lattice angles in degrees
   Adjective DeduceLatticeSystem(float a[3], float alpha[3])
   {
      using std::fabs;
      int sides = fabs((a[0]-a[1])<1.e-4) + fabs((a[0]-a[2])<1.e-4);
      int right = fabs((alpha[0]-90.)<1.e-4) + fabs((alpha[1]-90.)<1.e-4) + fabs((alpha[2]-90.)<1.e04);
      switch(right)
      {
      case 3:
         if(sides==2) return CUBIC;
         if(sides==1) return TETRAGONAL;
                      return ORTHORHOMBIC;
      case 2:
         if(sides==2) return RHOMBOHEDRAL;
         if(fabs((alpha[1]-90.)>1.e-4)) return MONOCLINIC;
         if(sides==1 && fabs((alpha[2]-120.)<1.e-4)) return HEXAGONAL;
      case 1:
      case 0:
         if(true)
         std::cout << "DeduceLatticeSystem(a,alpha) could not figure out lattice system for" << std::endl
                   << "lattice constants " << a[0] << " " << a[1] << " " << a[2] << std::endl
                   << "lattice angles " << alpha[0] << " " << alpha[1] << " " << alpha[2] << std::endl
                   << "equal sides = " << sides << std::endl
                   << "right angles = " << right << std::endl;
         return TRICLINIC;
      }
      // gamma = 60 is Trigonal?
      // What is ORTHOHEXAGONAL? 
     return TRICLINIC;
   }

   /// Possible lattice centerings
   enum CenteringType
   {
      PRIMITIVE    = 0,  // b0=(0,0,0)
      BODYCENTERED = 1,  // b0=(0,0,0), b1=(1/2,1/2,1/2)
      FACECENTERED = 2,  // b0=(0,0,0), b1=(0,1/2,1/2), b2=(1/2,0,1/2), b3=(1/2,1/2,0)
      BASECENTERED = 3   // b0=(0,0,0), b1=(1/2,1/2,0)
   };

   /// Compile-time name resolution
   template<int i> std::string CenteringString()                { return "invalid";      }
   template<> std::string CenteringString<PRIMITIVE>()          { return "primitive";    }
   template<> std::string CenteringString<BODYCENTERED>()       { return "bodycentered"; }
   template<> std::string CenteringString<FACECENTERED>()       { return "facecentered"; }
   template<> std::string CenteringString<BASECENTERED>()       { return "basecentered"; }
   
   /// Run-time name resolution
   std::string CenteringString(CenteringType c)
   {
      switch(c)
      {
         case PRIMITIVE:    return "primitive";    break;
         case BODYCENTERED: return "bodycentered"; break;
         case FACECENTERED: return "facecentered"; break;
         case BASECENTERED: return "basecentered"; break;
         defalut:           return "invalid";      break;
      }
   }

   /// Return centering type
   CenteringType ParseCenteringType(std::string name)
   {
      CenteringType ival = PRIMITIVE;
      for(int i=0; i<4; i++)
         if( name.compare(CenteringString(static_cast<CenteringType>(i)))==0 ) ival=static_cast<CenteringType>(i);
      return ival;
   }

   int DeduceCentering(Vec<float,3> p0, Vec<float,3> p1, Vec<float,3> p2)
   {
      using std::fabs;
      float mag0 = p0*p0;
      float mag1 = p1*p1;
      float mag2 = p2*p2;
      int axial = (fabs(mag0-1.00)<1.e-2) + (fabs(mag1-1.00)<1.e-2) + (fabs(mag2-1.00)<1.e-2);
      int face  = (fabs(mag0-0.50)<1.e-2) + (fabs(mag1-0.50)<1.e-2) + (fabs(mag2-0.50)<1.e-2);
      int body  = (fabs(mag0-0.75)<1.e-2) + (fabs(mag1-0.75)<1.e-2) + (fabs(mag2-0.75)<1.e-2);
      int check = axial + face + body;
      if( check!=3 )
      {
      std::cout << "DeduceCentering(p0,p1,p2) only classified " << check << " vectors" << std::endl
                   << "p0 = " << p0[0] << " " << p0[1] << " " << p0[2] << " |p0|^2=" << mag0 << std::endl
                   << "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << " |p1|^2=" << mag1 << std::endl
                   << "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << " |p2|^2=" << mag2 << std::endl
                   << "axial = " << axial << " face = " << face << " body = " << body << std::endl;
      }
      if(axial==3) return Sym::PRIMITIVE;
      if(face ==3) return Sym::FACECENTERED;
      if(body ==3) return Sym::BODYCENTERED;
      if(face ==2) return Sym::FACECENTERED;
      std::cout << "DeduceCentering(p0,p1,p2) could not classify" << std::endl
                << "p0 = " << p0[0] << " " << p0[1] << " " << p0[2] << " |p0|^2=" << mag0 << std::endl
                << "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << " |p1|^2=" << mag1 << std::endl
                << "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << " |p2|^2=" << mag2 << std::endl
                << "axial = " << axial << " face = " << face << " body = " << body << std::endl;
   }


   // ==========================================================================
   // ===== Point Groups
   // ==========================================================================

   /// The number of point groups
   enum { NPointGroup = 32 };

   /// Returns the local number of point group of a space group (ordered by International Table)
   int PointGroup(int itseqno) 
   {
           if( itseqno<  1 )  return -1;
      else if( itseqno<  2 )  return  0;  // Triclinic
      else if( itseqno<  3 )  return  1;
      else if( itseqno<  6 )  return  2;  // Monoclinic
      else if( itseqno< 10 )  return  3; 
      else if( itseqno< 16 )  return  4;
      else if( itseqno< 25 )  return  5;  // Orthorhombic 
      else if( itseqno< 47 )  return  6; 
      else if( itseqno< 75 )  return  7;
      else if( itseqno< 81 )  return  8;  // Tetragonal 
      else if( itseqno< 83 )  return  9; 
      else if( itseqno< 89 )  return 10; 
      else if( itseqno< 99 )  return 11; 
      else if( itseqno<111 )  return 12; 
      else if( itseqno<123 )  return 13; 
      else if( itseqno<143 )  return 14;
      else if( itseqno<147 )  return 15;  // Trigonal 
      else if( itseqno<149 )  return 16; 
      else if( itseqno<156 )  return 17; 
      else if( itseqno<162 )  return 18; 
      else if( itseqno<168 )  return 19;
      else if( itseqno<174 )  return 20;  // Hexagonal 
      else if( itseqno<175 )  return 21; 
      else if( itseqno<177 )  return 22; 
      else if( itseqno<183 )  return 23; 
      else if( itseqno<187 )  return 24; 
      else if( itseqno<191 )  return 25; 
      else if( itseqno<195 )  return 26;
      else if( itseqno<200 )  return 27;  // Cubic 
      else if( itseqno<207 )  return 28; 
      else if( itseqno<215 )  return 29; 
      else if( itseqno<221 )  return 30; 
      else if( itseqno<231 )  return 31;
      else                    return -1;
   }

   /// Compile-time resolution of the symbol associated with a point-group index
   template<int i> std::string SchoenfleisSymbol() { return "invalid"; }
   template<> std::string SchoenfleisSymbol<0>()   { return "C1";      } //  0 Triclinic
   template<> std::string SchoenfleisSymbol<1>()   { return "Ci";      } //  1
   template<> std::string SchoenfleisSymbol<2>()   { return "C2";      } //  2 Monoclinic
   template<> std::string SchoenfleisSymbol<3>()   { return "Cs";      } //  3
   template<> std::string SchoenfleisSymbol<4>()   { return "C2h";     } //  4
   template<> std::string SchoenfleisSymbol<5>()   { return "D2";      } //  5 Orthorhombic
   template<> std::string SchoenfleisSymbol<6>()   { return "C2v";     } //  6
   template<> std::string SchoenfleisSymbol<7>()   { return "D2h";     } //  7
   template<> std::string SchoenfleisSymbol<8>()   { return "C4";      } //  8 Tetragonal
   template<> std::string SchoenfleisSymbol<9>()   { return "S4";      } //  9
   template<> std::string SchoenfleisSymbol<10>()  { return "C4h";     } // 10
   template<> std::string SchoenfleisSymbol<11>()  { return "D4";      } // 11
   template<> std::string SchoenfleisSymbol<12>()  { return "C4v";     } // 12
   template<> std::string SchoenfleisSymbol<13>()  { return "D2d";     } // 13
   template<> std::string SchoenfleisSymbol<14>()  { return "D4h";     } // 14
   template<> std::string SchoenfleisSymbol<15>()  { return "C3";      } // 15 Trigonal
   template<> std::string SchoenfleisSymbol<16>()  { return "S6";      } // 16
   template<> std::string SchoenfleisSymbol<17>()  { return "D3";      } // 17
   template<> std::string SchoenfleisSymbol<18>()  { return "C3v";     } // 18
   template<> std::string SchoenfleisSymbol<19>()  { return "D3d";     } // 19
   template<> std::string SchoenfleisSymbol<20>()  { return "C6";      } // 20 Hexagonal
   template<> std::string SchoenfleisSymbol<21>()  { return "C3h";     } // 21
   template<> std::string SchoenfleisSymbol<22>()  { return "C6h";     } // 22
   template<> std::string SchoenfleisSymbol<23>()  { return "D6";      } // 23
   template<> std::string SchoenfleisSymbol<24>()  { return "C6v";     } // 24
   template<> std::string SchoenfleisSymbol<25>()  { return "D3h";     } // 25
   template<> std::string SchoenfleisSymbol<26>()  { return "D6h";     } // 26
   template<> std::string SchoenfleisSymbol<27>()  { return "T";       } // 27 Cubic
   template<> std::string SchoenfleisSymbol<28>()  { return "Th";      } // 28
   template<> std::string SchoenfleisSymbol<29>()  { return "O";       } // 29
   template<> std::string SchoenfleisSymbol<30>()  { return "Td";      } // 30
   template<> std::string SchoenfleisSymbol<31>()  { return "Oh";      } // 31

   /// Returns the symbol associated with a point-group index
   std::string SchoenfleisSymbol(int ipoint)
   {
      std::string const value[] =
      {
         "C1",   //  0 Triclinic
         "Ci",   //  1
         "C2",   //  2 Monoclinic
         "Cs",   //  3
         "C2h",  //  4
         "D2",   //  5 Orthorhombic
         "C2v",  //  6
         "D2h",  //  7
         "C4",   //  8 Tetragonal
         "S4",   //  9
         "C4h",  // 10
         "D4",   // 11
         "C4v",  // 12
         "D2d",  // 13
         "D4h",  // 14
         "C3",   // 15 Trigonal
         "S6",   // 16
         "D3",   // 17
         "C3v",  // 18
         "D3d",  // 19
         "C6",   // 20 Hexagonal
         "C3h",  // 21
         "C6h",  // 22
         "D6",   // 23
         "C6v",  // 24
         "D3h",  // 25
         "D6h",  // 26
         "T",    // 27 Cubic
         "Th",   // 28
         "O",    // 29
         "Td",   // 30
         "Oh"    // 31
      };
      if( ipoint<0 || ipoint>31 )
         return "invalid";
      else
         return value[ipoint];
   }

   int PointGroup(std::string schoenfleis)
   {
      std::string const value[] =
      {
         "C1",   //  0 Triclinic
         "Ci",   //  1
         "C2",   //  2 Monoclinic
         "Cs",   //  3
         "C2h",  //  4
         "D2",   //  5 Orthorhombic
         "C2v",  //  6
         "D2h",  //  7
         "C4",   //  8 Tetragonal
         "S4",   //  9
         "C4h",  // 10
         "D4",   // 11
         "C4v",  // 12
         "D2d",  // 13
         "D4h",  // 14
         "C3",   // 15 Trigonal
         "S6",   // 16
         "D3",   // 17
         "C3v",  // 18
         "D3d",  // 19
         "C6",   // 20 Hexagonal
         "C3h",  // 21
         "C6h",  // 22
         "D6",   // 23
         "C6v",  // 24
         "D3h",  // 25
         "D6h",  // 26
         "T",    // 27 Cubic
         "Th",   // 28
         "O",    // 29
         "Td",   // 30
         "Oh"    // 31
      };
      int i=0;
      while( i<32 && value[i].compare(schoenfleis)!=0) i++;
      return i;
   }

   void _pgscopy(std::vector<std::string>& result, std::string* array, int nel)
   {
      result.resize(nel);
      for(int i=0; i<nel; i++) result[i] = array[i];
   }

   // Return the short hand notation for symmetry operations of point group
   std::vector<std::string> PointGroupShorthand(int iptgroup)
   {
      std::string C1[1]   = { "x,y,z" };
      std::string Ci[2]   = { "x,y,z", "-x,-y,-z" };
      std::string C2[2]   = { "x,y,z", "-x,y,-z" };
      std::string Cs[2]   = { "x,y,z", "x,-y,z" };
      std::string C2h[4]  = { "x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z" };
      std::string D2[4]   = { "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z" };
      std::string C2v[4]  = { "x,y,z", "-x,-y,z", "x,-y,z", "-x,y,z" };
      std::string D2h[8]  = { "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "-x,-y,-z", "x,y,-z", "x,-y,z", "-x,y,z" };
      std::string C4[4]   = { "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z" };
      std::string S4[4]   = { "x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z" };
      std::string C4h[8]  = { "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z", "-x,-y,-z", "x,y,-z", "y,-x,-z", "-y,x,-z" };
      std::string D4[8]   =  { "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z", "-x,y,-z", "x,-y,-z", "y,x,-z", "-y,-x,-z" };
      std::string C4v[8]  = { "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z", "x,-y,z", "-x,y,z", "-y,-x,z", "y,x,z" };
      std::string D2d[8]  = { "x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z", "-x,y,-z", "x,-y,-z", "-y,-x,z", "y,x,z" };
      std::string D4h[16] = { "x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z", "-x,y,-z", "x,-y,-z", "y,x,-z", "-y,-x,-z", "-x,-y,-z", "x,y,-z", "y,-x,-z", "-y,x,-z", "x,-y,z", "-x,y,z", "-y,-x,z", "y,x,z" };
      std::string C3[3]   = { "x,y,z", "-y,x-y,z", "-x+y,-x,z" };
      std::string C3i[6]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z" };
      std::string D3[6]   = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z" };
      std::string C3v[6]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,z", "-x+y,y,z", "x,x-y,z" };
      std::string D3d[12] = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z", "y,x,z", "x-y,-y,z", "-x,-x+y,z" };
      std::string C6[6]   = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z" };
      std::string C3h[6]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z" };
      std::string C6h[12] = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z" };
      std::string D6[12]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z" };
      std::string C6v[12]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z", "-y,-x,z", "-x+y,y,z", "x,x-y,z", "y,x,z", "x-y,-y,z", "-x,-x+y,z" };
      std::string D3h[12]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z", "-y,-x,z", "-x+y,y,z", "x,x-y,z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z"};
      std::string D6h[24]  = { "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z", "-y,-x,z", "-x+y,y,z", "x,x-y,z", "y,x,z", "x-y,-y,z", "-x,-x+y,z" };
      std::string T[12]  = { "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "z,x,y", "z,-x,-y", "-z,-x,y", "-z,x,-y", "y,z,x", "-y,z,-x", "y,-z,-x", "-y,-z,x" };
  std::string Th[24]  = { "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "z,x,y", "z,-x,-y", "-z,-x,y", "-z,x,-y", "y,z,x", "-y,z,-x", "y,-z,-x", "-y,-z,x", "-x,-y,-z", "x,y,-z", "x,-y,z", "-x,y,z", "-z,-x,-y", "-z,x,y", "z,x,-y", "z,-x,y", "-y,-z,-x", "y,-z,x", "-y,z,x", "y,z,-x", };
      std::string  O[24] = { "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "z,-x,-y", "-z,x,-y", "-z,-x,y", "y,z,x", "-y,z,-x", "-y,-z,x", "y,-z,-x", "y,x,-z", "-y,-x,-z", "-y,x,z", "y,-x,z", "x,z,-y", "-x,z,y", "x,-z,y", "-x,-z,-y", "z,y,-x", "z,-y,x", "-z,-y,-x", "-z,y,x" };
      std::string Td[24] = { "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "z,-x,-y", "-z,x,-y", "-z,-x,y", "y,z,x", "-y,z,-x", "-y,-z,x", "y,-z,-x", "y,x,z", "-y,-x,z", "-y,x,-z", "y,-x,-z", "x,z,y", "-x,z,-y", "x,-z,-y", "-x,-z,y", "z,y,x", "z,-y,-x", "-z,-y,x", "-z,y,-x" };
      std::string Oh[48] = { "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "z,-x,-y", "-z,x,-y", "-z,-x,y", "y,z,x", "-y,z,-x", "-y,-z,x", "y,-z,-x", "y,x,-z", "-y,-x,-z", "-y,x,z", "y,-x,z", "x,z,-y", "-x,z,y", "x,-z,y", "-x,-z,-y", "z,y,-x", "z,-y,x", "-z,-y,-x", "-z,y,x", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "-z,-x,-y", "-z,x,y", "z,-x,y", "z,x,-y", "-y,-z,-x", "y,-z,x", "y,z,-x", "-y,z,x", "-y,-x,z", "y,x,z", "y,-x,-z", "-y,x,-z", "-x,-z,y", "x,-z,-y", "-x,z,-y", "x,z,y", "-z,-y,x", "-z,y,-x", "z,y,x", "z,-y,-x" };
      // create the vector
      std::vector<std::string> result;
      switch(iptgroup)
      {
      case  0: _pgscopy(result, C1,  1); break;
      case  1: _pgscopy(result, Ci,  2); break;
      case  2: _pgscopy(result, C2,  2); break;
      case  3: _pgscopy(result, Cs,  2); break;
      case  4: _pgscopy(result, C2h, 4); break;
      case  5: _pgscopy(result, D2,  4); break;
      case  6: _pgscopy(result, C2v, 4); break;
      case  7: _pgscopy(result, D2h, 8); break;
      case  8: _pgscopy(result, C4,  4); break;
      case  9: _pgscopy(result, S4,  4); break;
      case 10: _pgscopy(result, C4h, 8); break;
      case 11: _pgscopy(result, D4,  8); break;
      case 12: _pgscopy(result, C4v, 8); break;
      case 13: _pgscopy(result, D2d, 8); break;
      case 14: _pgscopy(result, D4h,16); break;
      case 15: _pgscopy(result, C3,  3); break;
      case 16: _pgscopy(result, C3i, 6); break;
      case 17: _pgscopy(result, D3,  6); break;
      case 18: _pgscopy(result, C3v, 6); break;
      case 19: _pgscopy(result, D3d,12); break;
      case 20: _pgscopy(result, C6,  6); break;
      case 21: _pgscopy(result, C3h, 6); break;
      case 22: _pgscopy(result, C6h,12); break;
      case 23: _pgscopy(result, D6, 12); break;
      case 24: _pgscopy(result, C6v,12); break;
      case 25: _pgscopy(result, D3h,12); break;
      case 26: _pgscopy(result, D6h,24); break;
      case 27: _pgscopy(result, T,  12); break;
      case 28: _pgscopy(result, Th, 24); break;
      case 29: _pgscopy(result, O,  24); break;
      case 30: _pgscopy(result, Td, 24); break;
      case 31: _pgscopy(result, Oh, 48); break;
      }
      return result;
   }


   // ==========================================================================
   // ===== Space Groups
   // ==========================================================================

   /// The total number of possible space groups
   enum { NSpaceGroup = 230 };

   /// Returns symbol associated with a particular ITA sequence number
   std::string HermannMauginSymbol(int itseqno);

   /// Returns the centering associated with each space group
   //    P — Primitive
   //    I — Body centered (from the German "Innenzentriert")
   //    F — Face centered (from the German "Flächenzentriert")
   //    A — Base centered on A faces only
   //    B — Base centered on B faces only
   //    C — Base centered on C faces only
   //    R — Rhombohedral
   int Centering(int itseqno)
   {
      int A[4]  = { 38, 39, 40, 41 };
      int C[16] = { 5, 8, 9, 12, 15, 20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68 };
      int F[16] = { 22, 42, 43, 69, 70, 196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228 };
      int I[38] = { 23, 24, 44, 45, 46, 71, 72, 73, 74, 79, 80, 82, 87, 88, 97, 98, 107, 108,
                    109, 110, 115, 120, 121, 122, 139, 140, 141, 142, 197, 199, 204, 206, 211,
                    214, 217, 220, 229, 230 };
      int R[7]  = { 146, 148, 155, 160, 161, 166, 167 };
      int P[149]= { 1, 2, 3, 4, 6, 7, 10, 11, 13, 14, 16, 17, 18, 19, 25, 26, 27, 28, 29, 30, 31,
                    32, 33, 34,   47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62,
                    75, 76, 77, 78, 81, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100,
                    101, 102, 103, 104, 105, 106, 111, 112, 113, 114, 116, 117, 118, 119, 123, 124,
                    125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 143, 144,
                    145, 147, 149, 150, 151, 152, 153, 154, 156, 157, 158, 159, 162, 163, 164, 165,
                    168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                    184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 198, 200, 201, 205,
                    207, 208, 212, 213, 215, 218, 221, 222, 223, 224 };
      // Test 
      bool found = false;
      for(int i=0; i<4  && !found; i++) found = A[i]==itseqno;
      for(int i=0; i<16 && !found; i++) found = C[i]==itseqno;
      if( found ) return BASECENTERED;
      for(int i=0; i<16 && !found; i++) found = F[i]==itseqno;
      if( found ) return FACECENTERED;
      for(int i=0; i<38 && !found; i++) found = I[i]==itseqno;
      if( found ) return BODYCENTERED;
      return PRIMITIVE;
   }

   /// Returns the crystal system of each space group
   int CrystalSystem(int itseqno) 
   {
           if(itseqno<  1)  return 10;
      else if(itseqno<  3)  return TRICLINIC;
      else if(itseqno< 16)  return MONOCLINIC;
      else if(itseqno< 75)  return ORTHORHOMBIC;
      else if(itseqno<143)  return TETRAGONAL;
      else if(itseqno<168)  return TRIGONAL;
      else if(itseqno<195)  return HEXAGONAL;
      else if(itseqno<231)  return CUBIC;
      else                  return 10;
   }

   /// Returns the lattice system of each space group
   int LatticeSystem(int itseqno)
   {
          if(itseqno==146 || 
             itseqno==148 || 
             itseqno==155 || 
             itseqno==160 || 
             itseqno==166 || 
             itseqno==167)  return RHOMBOHEDRAL; 
           if(itseqno<  1)  return 10;
      else if(itseqno<  3)  return TRICLINIC;
      else if(itseqno< 16)  return MONOCLINIC;
      else if(itseqno< 75)  return ORTHORHOMBIC;
      else if(itseqno<143)  return TETRAGONAL;
      else if(itseqno<195)  return HEXAGONAL;
      else if(itseqno<231)  return CUBIC;
      else                  return 10;
   }


// The Hermann-Maugin symbols indexed by International Table sequence number
// Adapted from http://www.ruppweb.org/ray/101index.html

// Reference:
// Hans Wondratschek, "Matrices, Mappings and Crystallographic Symmetry"
//
// *** The Fundamental Matrix of the Coordinate Basis G ***
// The matrix Gik is defined Gik = a[i]*a[k]*cos(alpha[j])
// The distance between two sites is r^2 = r_i G_ik r_k, with r_i = y_i - x_i
// The angle between two vectors is cos(Phi) = ( r_i G_ik t_k ) / sqrt( r^2 * t^2 )
// So G is the operator that finds the scalar product for sites specified by IVecs
// For an orthonormal baiss G = I
// The volume of a unit cell is det(G) = V^2
//
// A matrix is orthogonal if A^{-1} = A^T
//
// *** Definitions of Space group and Point group ***
// The set R of all symmetry operations of a crystal pattern is called the SPACE GROUP
// of the crystal pattern. The set of all element of R which leave a given point P fixed
// is the called the POINT GROUP (SITE SYMMETRY) S of P with regard to the space group R.
//
// HM = Hermann-Maugin Symbol
// The isometries of a POINT GROUP can be   (W^k = I defines the order k)
// 1. Identity, order 1
// 2. Inversion, order 2
// 3. Rotation,  order N=2,3,4, or 6      
//         
//      HM symbol N^j         Phi = j*360^o/N; 
//      -------------         ---------------
//           1^1    =    1   =     0^o
//           2^1    =    2   =   180^o
//           3^1    =    3   =   120^o
//           3^2    =        =   240^o
//           4^1    =    4   =    90^o
//           4^3    =        =   270^o
//           6^1    =    6   =    60^o
//           6^5    =        =   300^o,
// 4. Rotoinversions
//      \bar{N}^j              Phi = j*360^o/N; 
//      -------------         ---------------
//      \bar{1}^1 = inversion  =   0^o
//      \bar{2}^1 = reflection = 180^o
//      \bar{3}^1 = \bar{3}    = 120^o
//      \bar{3}^3 = inversion
//      \bar{3}^6 = I
//      \bar{4}^1 = \bar{4}    =  90^o
//      \bar{4}^2 = 2
//      \bar{4}^3 = \bar{4}^3  = 270^o
//      \bar{6}^1 = \bar{6}^1  =  60^o
//      \bar{6}^2 = 3
//      \bar{6}^3 = reflection 
// 5. Reflections, order 2
//      m = \bar{2}  
//

template<int i> struct SGTraits {};

// TRICLINIC
// POINT = 0
template<> struct SGTraits<1>
{
   enum { LATTICE = TRICLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 0 };
   static std::string HermannMaugin() { return "P 1"; }            //   1
};

template<> struct SGTraits<2>
{
   enum { LATTICE = TRICLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 0 };
   static std::string HermannMaugin() { return "P -1"; }          //   2
};

// MONOCLINIC
// POINT = 1
template<> struct SGTraits<3>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 1 };
   static std::string HermannMaugin() { return "P 2"; }           //   3
};

template<> struct SGTraits<4>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 1 };
   static std::string HermannMaugin() { return "P 21"; }          //   4
};

template<> struct SGTraits<5>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 1 };
   static std::string HermannMaugin() { return "C 2"; }           //   5
};

// POINT = 2
template<> struct SGTraits<6>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 2 };
   static std::string HermannMaugin() { return "P M"; }           //   6
};

template<> struct SGTraits<7>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 2 };
   static std::string HermannMaugin() { return "P C"; }           //   7
};

template<> struct SGTraits<8>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 2 };
   static std::string HermannMaugin() { return "C M"; }           //   8
};

template<> struct SGTraits<9>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 2 };
   static std::string HermannMaugin() { return "C C"; }           //   9
};

   // POINT = 3
template<> struct SGTraits<10>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "P 2/M"; }         //  10
};

template<> struct SGTraits<11>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "P 21/M"; }        //  11
};

template<> struct SGTraits<12>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "C 2/M"; }         //  12
};

template<> struct SGTraits<13>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "P 2/C"; }         //  13
};

template<> struct SGTraits<14>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "P 21/C"; }        //  14
};

template<> struct SGTraits<15>
{
   enum { LATTICE = MONOCLINIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 3 };
   static std::string HermannMaugin() { return "C 2/C"; }         //  15
};

// ORTHORHOMBIC
// POINT = 4
template<> struct SGTraits<16>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "P 2 2 2"; }       //  16
};

template<> struct SGTraits<17>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "P 2 2 21"; }      //  17
};

template<> struct SGTraits<18>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "P 21 21 2"; }     //  18
};

template<> struct SGTraits<19>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "P 21 21 21"; }    //  19
};

template<> struct SGTraits<20>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "C 2 2 21"; }      //  20
};

template<> struct SGTraits<21>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "C 2 2 2"; }       //  21
};

template<> struct SGTraits<22>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "F 2 2 2"; }       //  22
};

template<> struct SGTraits<23>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "I 2 2 2"; }       //  23
};

template<> struct SGTraits<24>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 4 };
   static std::string HermannMaugin() { return "I 21 21 21"; }    //  24
};

// POINT = 5
template<> struct SGTraits<25>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P M M 2"; }       //  25
};

template<> struct SGTraits<26>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P M C 21"; }      //  26
};

template<> struct SGTraits<27>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P C C 2"; }       //  27
};

template<> struct SGTraits<28>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P M A 2"; }       //  28
};

template<> struct SGTraits<29>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P C A 21"; }      //  29
};

template<> struct SGTraits<30>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P N C 2"; }       //  30
};

template<> struct SGTraits<31>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P M N 21"; }      //  31
};

template<> struct SGTraits<32>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P B A 2"; }       //  32
};

template<> struct SGTraits<33>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P N A 21"; }      //  33
};

template<> struct SGTraits<34>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "P N N 2"; }       //  34    
};

template<> struct SGTraits<35>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "C M M 2"; }       //  35
};

template<> struct SGTraits<36>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "C M C 21"; }      //  36
};

template<> struct SGTraits<37>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "C C C 2"; }       //  37
};

template<> struct SGTraits<38>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "A M M 2"; }       //  38
};

template<> struct SGTraits<39>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "A B M 2"; }       //  39
};

template<> struct SGTraits<40>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "A M A 2"; }       //  40
};

template<> struct SGTraits<41>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "A B A 2"; }       //  41
};

template<> struct SGTraits<42>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "F M M 2"; }       //  42
};

template<> struct SGTraits<43>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "F D D 2"; }       //  43
};

template<> struct SGTraits<44>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "I M M 2"; }       //  44
};

template<> struct SGTraits<45>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "I B A 2"; }       //  45
};

template<> struct SGTraits<46>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 5 };
   static std::string HermannMaugin() { return "I M A 2"; }       //  46
};

// POINT = 6
template<> struct SGTraits<47>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P M M M"; }       //  47
};

template<> struct SGTraits<48>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P N N N"; }       //  48
};

template<> struct SGTraits<49>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P C C M"; }       //  49
};

template<> struct SGTraits<50>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P B A N"; }       //  50
};

template<> struct SGTraits<51>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P M M A"; }       //  51
};

template<> struct SGTraits<52>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P N N A"; }       //  52
};

template<> struct SGTraits<53>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P M N A"; }       //  53
};

template<> struct SGTraits<54>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P C C A"; }       //  54
};

template<> struct SGTraits<55>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P B A M"; }       //  55
};

template<> struct SGTraits<56>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P C C N"; }       //  56
};

template<> struct SGTraits<57>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P B C M"; }       //  57
};

template<> struct SGTraits<58>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P N N M"; }       //  58
};

template<> struct SGTraits<59>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P M M N"; }       //  59
};

template<> struct SGTraits<60>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P B C N"; }       //  60
};

template<> struct SGTraits<61>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P B C A"; }       //  61
};

template<> struct SGTraits<62>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "P N M A"; }       //  62
};

template<> struct SGTraits<63>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C M C M"; }       //  63
};

template<> struct SGTraits<64>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C M C A"; }       //  64
};

template<> struct SGTraits<65>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C M M M"; }       //  65
};

template<> struct SGTraits<66>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C C C M"; }       //  66
};

template<> struct SGTraits<67>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C M M A"; }       //  67
};

template<> struct SGTraits<68>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BASECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "C C C A"; }       //  68
};

template<> struct SGTraits<69>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "F M M M"; }       //  69
};

template<> struct SGTraits<70>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "F D D D"; }       //  70
};

template<> struct SGTraits<71>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "I M M M"; }       //  71
};

template<> struct SGTraits<72>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "I B A M"; }       //  72
};

template<> struct SGTraits<73>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "I B C A"; }       //  73
};

template<> struct SGTraits<74>
{
   enum { LATTICE = ORTHORHOMBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 6 };
   static std::string HermannMaugin() { return "I M M A"; }       //  74
};

// TETRAGONAL
// POINT = 7
template<> struct SGTraits<75>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "P 4"; }           //  75
};

template<> struct SGTraits<76>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "P 41"; }          //  76
};

template<> struct SGTraits<77>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "P 42"; }          //  77
};

template<> struct SGTraits<78>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "P 43"; }          //  78
};

template<> struct SGTraits<79>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "I 4"; }           //  79
};

template<> struct SGTraits<80>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 7 };
   static std::string HermannMaugin() { return "I 41"; }          //  80
};

// POINT = 8
template<> struct SGTraits<81>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 8 };
   static std::string HermannMaugin() { return "P -4"; }          //  81
};

template<> struct SGTraits<82>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 8 };
   static std::string HermannMaugin() { return "I -4"; }          //  82
};

// POINT = 9
template<> struct SGTraits<83>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "P 4/M"; }         //  83
};

template<> struct SGTraits<84>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "P 42/M"; }        //  84
};

template<> struct SGTraits<85>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "P 4/N"; }         //  85
};

template<> struct SGTraits<86>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "P 42/N"; }        //  86
};

template<> struct SGTraits<87>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "I 4/M"; }         //  87
};

template<> struct SGTraits<88>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 9 };
   static std::string HermannMaugin() { return "I 41/A"; }        //  88
};

// POINT = 10
template<> struct SGTraits<89>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 4 2 2"; }       //  89
};

template<> struct SGTraits<90>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 4 21 2"; }      //  90
};

template<> struct SGTraits<91>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 41 2 2"; }      //  91
};

template<> struct SGTraits<92>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 41 21 2"; }     //  92
};

template<> struct SGTraits<93>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 42 2 2"; }      //  93
};

template<> struct SGTraits<94>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 42 21 2"; }     //  94
};

template<> struct SGTraits<95>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 43 2 2"; }      //  95
};

template<> struct SGTraits<96>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "P 43 21 2"; }     //  96
};

template<> struct SGTraits<97>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "I 4 2 2"; }       //  97
};

template<> struct SGTraits<98>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 10 };
   static std::string HermannMaugin() { return "I 41 2 2"; }      //  98
};

// POINT = 11
template<> struct SGTraits<99>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 4 M M"; }       //  99
};

template<> struct SGTraits<100>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 4 B M"; }       // 100
};

template<> struct SGTraits<101>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 42 C M"; }      // 101
};

template<> struct SGTraits<102>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 42 N M"; }      // 102
};

template<> struct SGTraits<103>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 4 C C"; }       // 103
};

template<> struct SGTraits<104>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 4 N C"; }       // 104
};

template<> struct SGTraits<105>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 42 M C"; }      // 105
};

template<> struct SGTraits<106>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "P 42 B C"; }      // 106
};

template<> struct SGTraits<107>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "I 4 M M"; }       // 107
};

template<> struct SGTraits<108>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "I 4 C M"; }       // 108
};

template<> struct SGTraits<109>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "I 41 M D"; }      // 109
};

template<> struct SGTraits<110>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 11 };
   static std::string HermannMaugin() { return "I 41 C D"; }      // 110
};

// POINT = 12
template<> struct SGTraits<111>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 2 M"; }      // 111
};

template<> struct SGTraits<112>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 2 C"; }      // 112
};

template<> struct SGTraits<113>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 21 M"; }     // 113
};

template<> struct SGTraits<114>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 21 C"; }     // 114
};

template<> struct SGTraits<115>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "I -4 M 2"; }      // 115
};

template<> struct SGTraits<116>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 C 2"; }      // 116
};

template<> struct SGTraits<117>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 B 2"; }      // 117
};

template<> struct SGTraits<118>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 N 2"; }      // 118
};

template<> struct SGTraits<119>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "P -4 M 2"; }      // 119
};

template<> struct SGTraits<120>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "I -4 C 2"; }      // 120
};

template<> struct SGTraits<121>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "I -4 2 M"; }      // 121
};

template<> struct SGTraits<122>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 12 };
   static std::string HermannMaugin() { return "I -4 2 D"; }      // 122
};

// POINT = 13
template<> struct SGTraits<123>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/M M M"; }     // 123
};

template<> struct SGTraits<124>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/M C C"; }     // 124
};

template<> struct SGTraits<125>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/N B M"; }     // 125
};

template<> struct SGTraits<126>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/N N C"; }     // 126
};

template<> struct SGTraits<127>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/M B M"; }     // 127
};

template<> struct SGTraits<128>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/M N C"; }     // 128
};

template<> struct SGTraits<129>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/N M M"; }     // 129
};

template<> struct SGTraits<130>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 4/N C C"; }     // 130
};

template<> struct SGTraits<131>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/M M C"; }    // 131
};

template<> struct SGTraits<132>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/M C M"; }    // 132
};

template<> struct SGTraits<133>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/N B C"; }    // 133
};

template<> struct SGTraits<134>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/N N M"; }    // 134
};

template<> struct SGTraits<135>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/M B C"; }    // 135
};

template<> struct SGTraits<136>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/M N M"; }    // 136
};

template<> struct SGTraits<137>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/N M C"; }    // 137
};

template<> struct SGTraits<138>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "P 42/N C M"; }    // 138
};

template<> struct SGTraits<139>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "I 4/M M M"; }     // 139
};

template<> struct SGTraits<140>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "I 4/M C M"; }     // 140
};

template<> struct SGTraits<141>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "I 41/A M D"; }    // 141
};

template<> struct SGTraits<142>
{
   enum { LATTICE = TETRAGONAL };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 13 };
   static std::string HermannMaugin() { return "I 41/A C D"; }    // 142
};

// TRIGONAL
// POINT = 14
template<> struct SGTraits<143>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 14 };
   static std::string HermannMaugin() { return "P 3"; }           // 143
};

template<> struct SGTraits<144>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 14 };
   static std::string HermannMaugin() { return "P 31"; }          // 144
};

template<> struct SGTraits<145>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 14 };
   static std::string HermannMaugin() { return "P 32"; }          // 145
};

template<> struct SGTraits<146>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 14 };
   static std::string HermannMaugin() { return "R 3"; }           // 146
};

// POINT = 15
template<> struct SGTraits<147>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 15 };
   static std::string HermannMaugin() { return "P -3"; }          // 147
};

template<> struct SGTraits<148>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 15 };
   static std::string HermannMaugin() { return "R -3"; }          // 148
};

template<> struct SGTraits<149>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 15 };
   static std::string HermannMaugin() { return "P 3 1 2"; }       // 149
};

// POINT = 16
template<> struct SGTraits<150>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "P 3 2 1"; }       // 150
};

template<> struct SGTraits<151>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "P 31 1 2"; }      // 151
};

template<> struct SGTraits<152>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "P 31 2 1"; }      // 152
};

template<> struct SGTraits<153>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "P 32 1 2"; }      // 153
};

template<> struct SGTraits<154>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "P 32 2 1"; }      // 154
};

template<> struct SGTraits<155>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 16 };
   static std::string HermannMaugin() { return "R 3 2"; }         // 155
};

// POINT = 17
template<> struct SGTraits<156>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "P 3 M 1"; }       // 156
};

template<> struct SGTraits<157>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "P 3 1 M"; }       // 157
};

template<> struct SGTraits<158>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "P 3 C 1"; }       // 158
};

template<> struct SGTraits<159>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "P 3 1 C"; }       // 159
};

template<> struct SGTraits<160>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "R 3 M"; }         // 160
};

template<> struct SGTraits<161>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 17 };
   static std::string HermannMaugin() { return "R 3 C"; }         // 161
};

// POINT = 18
template<> struct SGTraits<162>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "P -3 1 M"; }      // 162
};

template<> struct SGTraits<163>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "P -3 1 C"; }      // 163
};

template<> struct SGTraits<164>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "P -3 M 1"; }      // 164
};

template<> struct SGTraits<165>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "P -3 C 1"; }      // 165
};

template<> struct SGTraits<166>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "R -3 M"; }        // 166
};

template<> struct SGTraits<167>
{
   enum { LATTICE = TRIGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 18 };
   static std::string HermannMaugin() { return "R -3 C"; }        // 167
};

// HEXAGONAL
// POINT = 19
template<> struct SGTraits<168>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 6"; }           // 168
};

template<> struct SGTraits<169>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 61"; }          // 169
};

template<> struct SGTraits<170>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 65"; }          // 170
};

template<> struct SGTraits<171>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 62"; }          // 171
};

template<> struct SGTraits<172>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 64"; }          // 172
};

template<> struct SGTraits<173>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 19 };
   static std::string HermannMaugin() { return "P 63"; }          // 173
};

// POINT = 20
template<> struct SGTraits<174>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 20 };
   static std::string HermannMaugin() { return "P -6"; }          // 174
};

// POINT = 21
template<> struct SGTraits<175>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 21 };
   static std::string HermannMaugin() { return "P 6/M"; }         // 175
};

template<> struct SGTraits<176>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 21 };
   static std::string HermannMaugin() { return "P 63/M"; }        // 176
};

// POINT = 22
template<> struct SGTraits<177>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 6 2 2"; }       // 177
};

template<> struct SGTraits<178>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 61 2 2"; }      // 178
};

template<> struct SGTraits<179>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 65 2 2"; }      // 179
};

template<> struct SGTraits<180>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 62 2 2"; }      // 180
};

template<> struct SGTraits<181>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 64 2 2"; }      // 181
};

template<> struct SGTraits<182>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 22 };
   static std::string HermannMaugin() { return "P 63 2 2"; }      // 182
};

// POINT = 23
template<> struct SGTraits<183>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 23 };
   static std::string HermannMaugin() { return "P 6 M M"; }       // 183
};

template<> struct SGTraits<184>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 23 };
   static std::string HermannMaugin() { return "P 6 C C"; }       // 184
};

template<> struct SGTraits<185>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 23 };
   static std::string HermannMaugin() { return "P 63 C M"; }      // 185
};

template<> struct SGTraits<186>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 23 };
   static std::string HermannMaugin() { return "P 63 M C"; }      // 186
};

// POINT = 24
template<> struct SGTraits<187>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 24 };
   static std::string HermannMaugin() { return "P -6 M 2"; }      // 187
};

template<> struct SGTraits<188>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 24 };
   static std::string HermannMaugin() { return "P -6 C 2"; }      // 188
};

template<> struct SGTraits<189>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 24 };
   static std::string HermannMaugin() { return "P -6 2 M"; }      // 189
};

template<> struct SGTraits<190>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 24 };
   static std::string HermannMaugin() { return "P -6 2 C"; }      // 190
};

// POINT = 25
template<> struct SGTraits<191>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 25 };
   static std::string HermannMaugin() { return "P 6/M M M"; }     // 191
};

template<> struct SGTraits<192>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 25 };
   static std::string HermannMaugin() { return "P 6/M C C"; }     // 192
};

template<> struct SGTraits<193>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 25 };
   static std::string HermannMaugin() { return "P 63/M C M"; }    // 193
};

template<> struct SGTraits<194>
{
   enum { LATTICE = HEXAGONAL };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 25 };
   static std::string HermannMaugin() { return "P 63/M M C"; }    // 194
};

// CUBIC (minus sign in front of triade optional)
// POINT 26
template<> struct SGTraits<195>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 26 };
   static std::string HermannMaugin() { return "P 2 3"; }         // 195
};

template<> struct SGTraits<196>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 26 };
   static std::string HermannMaugin() { return "F 2 3"; }         // 196
};

template<> struct SGTraits<197>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 26 };
   static std::string HermannMaugin() { return "I 2 3"; }         // 197
};

template<> struct SGTraits<198>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 26 };
   static std::string HermannMaugin() { return "P 21 3"; }        // 198
};

template<> struct SGTraits<199>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 26 };
   static std::string HermannMaugin() { return "I 21 3"; }        // 199
};

// POINT = 27
template<> struct SGTraits<200>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "P M 3"; }         // 200
};

template<> struct SGTraits<201>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "P N 3"; }         // 201
};

template<> struct SGTraits<202>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "F M 3"; }         // 202
};

template<> struct SGTraits<203>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "F D 3"; }         // 203
};

template<> struct SGTraits<204>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "I M 3"; }         // 204
};

template<> struct SGTraits<205>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "P A 3"; }         // 205
};

template<> struct SGTraits<206>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 27 };
   static std::string HermannMaugin() { return "I A 3"; }         // 206
};

// POINT = 28
template<> struct SGTraits<207>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "P 4 3 2"; }       // 207
};

template<> struct SGTraits<208>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "P 42 3 2"; }      // 208
};

template<> struct SGTraits<209>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "F 4 3 2"; }       // 209
};

template<> struct SGTraits<210>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "F 41 3 2"; }      // 210
};

template<> struct SGTraits<211>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "I 4 3 2"; }       // 211
};

template<> struct SGTraits<212>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "P 43 3 2"; }      // 212
};

template<> struct SGTraits<213>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "P 41 3 2"; }      // 213
};

template<> struct SGTraits<214>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 28 };
   static std::string HermannMaugin() { return "I 41 3 2"; }      // 214
};

// POINT = 29
template<> struct SGTraits<215>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 29 };
   static std::string HermannMaugin() { return "P -4 3 M"; }      // 215
};

template<> struct SGTraits<216>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 29 };
   static std::string HermannMaugin() { return "F -4 3 M"; }      // 216
};

template<> struct SGTraits<217>
{
   enum { LATTICE = CUBIC };
   enum { CENTER  = BODYCENTERED };
   enum { POINT   = 29 };
   enum { BRAVAIS = 11 };
   static std::string HermannMaugin() { return "I -4 3 M"; }      // 217
};

template<> struct SGTraits<218>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 29 };
   static std::string HermannMaugin() { return "P -4 3 N"; }      // 218
};

template<> struct SGTraits<219>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 29 };
   static std::string HermannMaugin() { return "F -4 3 C"; }      // 219
};

template<> struct SGTraits<220>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 29 };
   static std::string HermannMaugin() { return "I -4 3 D"; }      // 220
};

// POINT = 30
template<> struct SGTraits<221>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "P M 3 M"; }       // 221
};

template<> struct SGTraits<222>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "P N 3 N"; }       // 222
};

template<> struct SGTraits<223>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "P M 3 N"; }       // 223
};

template<> struct SGTraits<224>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = PRIMITIVE };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "P N 3 M"; }       // 224
};

template<> struct SGTraits<225>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "F M 3 M"; }       // 225
};

template<> struct SGTraits<226>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "F M 3 C"; }       // 226
};

template<> struct SGTraits<227>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "F D 3 M"; }       // 227
};

template<> struct SGTraits<228>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = FACECENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "F D 3 C"; }       // 228
};

template<> struct SGTraits<229>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "I M 3 M"; }       // 229
};

template<> struct SGTraits<230>
{
   enum { LATTICE = CUBIC };
   enum { CENTER = BODYCENTERED };
   enum { POINT = 30 };
   static std::string HermannMaugin() { return "I A 3 D"; }       // 230
};


// References:
// Hans Wondratschek, "Matrices, Mappings and Crystallographic Symmetry"
// Bilboa Crystallographic Server:  http://www.cryst.ehu.es
// NRL pages: http://cst-www.nrl.navy.mil/lattice

// In general, symmetry operations of crystallography are of the form
// \tilde{x} = W*x + w, with \tilde{x}, x, and w vectors; W a 3x3 matrix 
// or SEITZ notation  \tilde{x} = (W,w)*x
// Law of Composition:  (W,w) = (V,v)*(U,u) = (V*U,V*u+v)
// Reversion: (W,w)^{-1} = (W^{-1},-W^{-1}*w)
//
// Vectors are not affected by the translation w.
// If only use primitive lattice basis then only need W, and all elements are integers.
// 

template<typename TA, typename TB, typename TC>
void sym_MatrixMultiply(TA const A[], TB const B[], TC C[])
{
     C[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6]; C[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7]; C[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
     C[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6]; C[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7]; C[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
     C[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6]; C[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7]; C[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}



class Operator
{
public:

   /// Used as the denominator to represent column elements as a fraction
   //  Maximum denominator is 6, 60 = 2^2 * 3 * 5
   enum { DENOM = 60 };

   ~Operator() { ; }

   // Default constructor is identity operator
   Operator();

   // Construction through explicit specification of elements
   Operator(int  W11,  int  W12,  int  W13,
            int  W21,  int  W22,  int  W23,
            int  W31,  int  W32,  int  W33, 
            int  w1=0, int  w2=0, int  w3=0,
            int  determinant=1);

   // Construction through international table shorthand
   Operator(std::string shorthand) { this->assign(shorthand); }
   Operator(const char* shorthand) { this->assign(std::string(shorthand)); }

   // Construction through composition
   Operator(const Operator& U, const Operator& V);
   
   // Copy Constructor
   Operator(const Operator& S);
   
   // Assign the values of the elements  
   void assign(int  W11,  int  W12,  int  W13,
               int  W21,  int  W22,  int  W23,
               int  W31,  int  W32,  int  W33,
               int  w1=0, int  w2=0, int  w3=0,
               int  determinant=1);

   // Assign by shorthand notation
   void assign(std::string shorthand);

   // Assign only the matrix part
   void assign(const Operator& M);

   // Assign only the vector part
   void assign(const int column[3]);

   // Assign matrix and vector part
   void assign(const Operator& M, const int column[3]);

   // Apply the operator to a vector expressed in terms of primitive basis vectors
   void operator()(const int x[3], int r[3], bool posvec=false) const;
   void operator()(const Vec<int,3>& x, Vec<int,3>& r, bool posvec=false) const { this->operator()(&(x[0]),&(r[0]),posvec); }

   /// Determinant of matrix representation
   int det() const;

   /// Trace of the matrix representation
   int trace() const;

   /// Type of (improper) rotation associated with operation
   int type() const;

   /// Order of the symmetry operation, W^k = I
   int order() const;

   /// Inverse or reversion of symmetry operation
   Operator inv() const;

   /// Transpose of a symmetry matrix
   Operator T() const;
   
   /// check for equality
   bool operator==(const Operator& M) const;

   /// check for nonequality
   bool operator!=(const Operator& M) const;
 
   // These are legacy from somewhere 
   static int  frac_to_sixtieths(std::string frac);
   static void parse_shorthand(std::string sh_general, std::string& sh_matrix, int sixtieths[3]);

private:
   
   int e[12];   // The Elements of the transformation
   
};

Operator::Operator()
{
   this->assign(1,0,0,0,1,0,0,0,1);
}

Operator::Operator(int  W11,  int  W12,  int  W13,
                   int  W21,  int  W22,  int  W23,
                   int  W31,  int  W32,  int  W33, 
                   int  w1  , int  w2  , int  w3  ,
                   int  determinant)
{ this->assign(W11,W12,W13,W21,W22,W23,W33,w1,w2,w3, determinant); }

Operator::Operator(const Operator& A, const Operator& B) 
{ 
   // (A,a)(B,b) = ( AB, Ab+a )
   sym_MatrixMultiply(A.e,B.e,e); 
   e[ 9] = ( A.e[0]*B.e[ 9] + A.e[1]*B.e[10] + A.e[2]*B.e[11] + A.e[ 9] ) % DENOM;
   e[10] = ( A.e[3]*B.e[ 9] + A.e[4]*B.e[10] + A.e[5]*B.e[11] + A.e[10] ) % DENOM;
   e[11] = ( A.e[6]*B.e[ 9] + A.e[7]*B.e[10] + A.e[8]*B.e[11] + A.e[11] ) % DENOM;
   if(e[ 9]<0) e[ 9]+=DENOM;  
   if(e[10]<0) e[10]+=DENOM;  
   if(e[11]<0) e[11]+=DENOM;
}

Operator::Operator(const Operator& M) 
{ 
   assign(M);
}

void Operator::assign(int  W11,  int  W12,  int  W13,
                      int  W21,  int  W22,  int  W23,
                      int  W31,  int  W32,  int  W33, 
                      int  w1  , int  w2  , int  w3  ,
                      int  d)
{
      e[0]=d*W11; e[1]=W12/d; e[2]=W13/d;  e[3]=W21/d; e[4]=W22/d; e[5]=W23/d;  e[6]=W31/d; e[7]=W32/d; e[8]=W33/d; 
      e[ 9]= w1 % DENOM;
      e[10]= w2 % DENOM;
      e[11]= w3 % DENOM;
}


void Operator::assign(const Operator& M)
{
   for(int i=0; i<12; i++) 
      e[i]=M.e[i]; 
}

void Operator::assign(const Operator& M, const int* column) 
{ 
   for(int i=0; i<9; i++) 
      e[i]=M.e[i]; 
   assign(column);
}

void Operator::assign(const int* column )
{
   for(int i=0; i<3; i++)
      e[9+i] = column[i] % DENOM;
}

void Operator::operator()(const int* v, int* result, bool pos_vector) const 
{ 
   result[0] = e[0]*v[0]+e[1]*v[1]+e[2]*v[2];
   result[1] = e[3]*v[0]+e[4]*v[1]+e[5]*v[2];
   result[2] = e[6]*v[0]+e[7]*v[1]+e[8]*v[2];
   if( pos_vector ) 
   {
      result[0]  += e[9]; 
      result[1] += e[10]; 
      result[2] += e[11];
   }
}

void Operator::assign(std::string shorthand)
{
   // x,y,z
   if( shorthand.compare("x,y,z"   ) == 0 ) this->assign( 1, 0, 0,  0, 1, 0,  0, 0, 1 );
   if( shorthand.compare("-x,y,z"  ) == 0 ) this->assign(-1, 0, 0,  0, 1, 0,  0, 0, 1 );
   if( shorthand.compare("x,-y,z"  ) == 0 ) this->assign( 1, 0, 0,  0,-1, 0,  0, 0, 1 );
   if( shorthand.compare("-x,-y,z" ) == 0 ) this->assign(-1, 0, 0,  0,-1, 0,  0, 0, 1 );
   if( shorthand.compare("x,y,-z"  ) == 0 ) this->assign( 1, 0, 0,  0, 1, 0,  0, 0,-1 );
   if( shorthand.compare("-x,y,-z" ) == 0 ) this->assign(-1, 0, 0,  0, 1, 0,  0, 0,-1 );
   if( shorthand.compare("x,-y,-z" ) == 0 ) this->assign( 1, 0, 0,  0,-1, 0,  0, 0,-1 );
   if( shorthand.compare("-x,-y,-z") == 0 ) this->assign(-1, 0, 0,  0,-1, 0,  0, 0,-1 );
   // z, x, y
   if( shorthand.compare("z,x,y"   ) == 0 ) this->assign( 0, 0, 1,  1, 0, 0,  0, 1, 0 );
   if( shorthand.compare("-z,x,y"  ) == 0 ) this->assign( 0, 0,-1,  1, 0, 0,  0, 1, 0 );
   if( shorthand.compare("z,-x,y"  ) == 0 ) this->assign( 0, 0, 1, -1, 0, 0,  0, 1, 0 );
   if( shorthand.compare("-z,-x,y" ) == 0 ) this->assign( 0, 0,-1, -1, 0, 0,  0, 1, 0 );
   if( shorthand.compare("z,x,-y"  ) == 0 ) this->assign( 0, 0, 1,  1, 0, 0,  0,-1, 0 );
   if( shorthand.compare("-z,x,-y" ) == 0 ) this->assign( 0, 0,-1,  1, 0, 0,  0,-1, 0 );
   if( shorthand.compare("z,-x,-y" ) == 0 ) this->assign( 0, 0, 1, -1, 0, 0,  0,-1, 0 );
   if( shorthand.compare("-z,-x,-y") == 0 ) this->assign( 0, 0,-1, -1, 0, 0,  0,-1, 0 );
   // y, z, x 
   if( shorthand.compare("y,z,x"   ) == 0 ) this->assign( 0, 1, 0,  0, 0, 1,  1, 0, 0 );
   if( shorthand.compare("-y,z,x"  ) == 0 ) this->assign( 0,-1, 0,  0, 0, 1,  1, 0, 0 );
   if( shorthand.compare("y,-z,x"  ) == 0 ) this->assign( 0, 1, 0,  0, 0,-1,  1, 0, 0 );
   if( shorthand.compare("-y,-z,x" ) == 0 ) this->assign( 0,-1, 0,  0, 0,-1,  1, 0, 0 );
   if( shorthand.compare("y,z,-x"  ) == 0 ) this->assign( 0, 1, 0,  0, 0, 1, -1, 0, 0 );
   if( shorthand.compare("-y,z,-x" ) == 0 ) this->assign( 0,-1, 0,  0, 0, 1, -1, 0, 0 );
   if( shorthand.compare("y,-z,-x" ) == 0 ) this->assign( 0, 1, 0,  0, 0,-1, -1, 0, 0 );
   if( shorthand.compare("-y,-z,-x") == 0 ) this->assign( 0,-1, 0,  0, 0,-1, -1, 0, 0 );
   /// y, x, z
   if( shorthand.compare("y,x,z"   ) == 0 ) this->assign( 0, 1, 0,  1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("-y,x,z"  ) == 0 ) this->assign( 0,-1, 0,  1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("y,-x,z"  ) == 0 ) this->assign( 0, 1, 0, -1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("-y,-x,z" ) == 0 ) this->assign( 0,-1, 0, -1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("y,x,-z"  ) == 0 ) this->assign( 0, 1, 0,  1, 0, 0,  0, 0,-1 );
   if( shorthand.compare("-y,x,-z" ) == 0 ) this->assign( 0,-1, 0,  1, 0, 0,  0, 0,-1 );
   if( shorthand.compare("y,-x,-z" ) == 0 ) this->assign( 0, 1, 0, -1, 0, 0,  0, 0,-1 );
   if( shorthand.compare("-y,-x,-z") == 0 ) this->assign( 0,-1, 0, -1, 0, 0,  0, 0,-1 );
   /// z, y, x
   if( shorthand.compare("z,y,x"   ) == 0 ) this->assign( 0, 0, 1,  0, 1, 0,  1, 0, 0 );
   if( shorthand.compare("-z,y,x"  ) == 0 ) this->assign( 0, 0,-1,  0, 1, 0,  1, 0, 0 );
   if( shorthand.compare("z,-y,x"  ) == 0 ) this->assign( 0, 0, 1,  0,-1, 0,  1, 0, 0 );
   if( shorthand.compare("-z,-y,x" ) == 0 ) this->assign( 0, 0,-1,  0,-1, 0,  1, 0, 0 );
   if( shorthand.compare("z,y,-x"  ) == 0 ) this->assign( 0, 0, 1,  0, 1, 0, -1, 0, 0 );
   if( shorthand.compare("-z,y,-x" ) == 0 ) this->assign( 0, 0,-1,  0, 1, 0, -1, 0, 0 );
   if( shorthand.compare("z,-y,-x" ) == 0 ) this->assign( 0, 0, 1,  0,-1, 0, -1, 0, 0 );
   if( shorthand.compare("-z,-y,-x") == 0 ) this->assign( 0, 0,-1,  0,-1, 0, -1, 0, 0 );
   /// x, z, y
   if( shorthand.compare("x,z,y"   ) == 0 ) this->assign( 1, 0, 0,  0, 0, 1,  0, 1, 0 );
   if( shorthand.compare("-x,z,y"  ) == 0 ) this->assign(-1, 0, 0,  0, 0, 1,  0, 1, 0 );
   if( shorthand.compare("x,-z,y"  ) == 0 ) this->assign( 1, 0, 0,  0, 0,-1,  0, 1, 0 );
   if( shorthand.compare("-x,-z,y" ) == 0 ) this->assign(-1, 0, 0,  0, 0,-1,  0, 1, 0 );
   if( shorthand.compare("x,z,-y"  ) == 0 ) this->assign( 1, 0, 0,  0, 0, 1,  0,-1, 0 );
   if( shorthand.compare("-x,z,-y" ) == 0 ) this->assign(-1, 0, 0,  0, 0, 1,  0,-1, 0 );
   if( shorthand.compare("x,-z,-y" ) == 0 ) this->assign( 1, 0, 0,  0, 0,-1,  0,-1, 0 );
   if( shorthand.compare("-x,-z,-y") == 0 ) this->assign(-1, 0, 0,  0, 0,-1,  0,-1, 0 );
   /// x-y, x, z
   if( shorthand.compare("x-y,x,z" ) ==   0 ) this->assign( 1,-1, 0,  1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("x-y,x,-z") ==   0 ) this->assign( 1,-1, 0,  1, 0, 0,  0, 0,-1 );
   if( shorthand.compare("x-y,-y,z" ) ==  0 ) this->assign( 1,-1, 0,  0,-1, 0,  0, 0, 1 );
   if( shorthand.compare("x-y,-y,-z") ==  0 ) this->assign( 1,-1, 0,  0,-1, 0,  0, 0,-1 );
   if( shorthand.compare("-x+y,-x,z") ==  0 ) this->assign(-1, 1, 0, -1, 0, 0,  0, 0, 1 );
   if( shorthand.compare("-x+y,-x,-z") == 0 ) this->assign(-1, 1, 0, -1, 0, 0,  0, 0,-1 );
   if( shorthand.compare("-x+y,y,z") ==   0 ) this->assign(-1, 1, 0,  0, 1, 0,  0, 0, 1 );
   if( shorthand.compare("-x+y,y,-z") ==  0 ) this->assign(-1, 1, 0,  0, 1, 0,  0, 0,-1 );
   /// x, x-y, ,z
   if( shorthand.compare("x,x-y,z" ) ==   0 ) this->assign( 1, 0, 0,  1,-1, 0,  0, 0, 1 );
   if( shorthand.compare("x,x-y,-z") ==   0 ) this->assign( 1, 0, 0,  1,-1, 0,  0, 0,-1 );
   if( shorthand.compare("-y,x-y,z" ) ==  0 ) this->assign( 0,-1, 0,  1,-1, 0,  0, 0, 1 );
   if( shorthand.compare("-y,x-y,-z") ==  0 ) this->assign( 0,-1, 0,  1,-1, 0,  0, 0,-1 );
   if( shorthand.compare("-x,-x+y,z") ==  0 ) this->assign(-1, 0, 0, -1, 1, 0,  0, 0, 1 );
   if( shorthand.compare("-x,-x+y,-z") == 0 ) this->assign(-1, 0, 0, -1, 1, 0,  0, 0,-1 );
   if( shorthand.compare("y,-x+y,z") ==   0 ) this->assign( 0, 1, 0, -1, 1, 0,  0, 0, 1 );
   if( shorthand.compare("y,-x+y,-z") ==  0 ) this->assign( 0, 1, 0, -1, 1, 0,  0, 0,-1 );
}

int Operator::frac_to_sixtieths(std::string frac)
{
   int islash = frac.find_first_of("/");
   int numer = 0;
   if( islash>1 && islash!=std::string::npos )
   {
      std::istringstream nstr( frac.substr(0,islash) );
      nstr >> numer;
   }
   int denom = 1;
   islash++;
   if( islash<frac.length() && islash!=std::string::npos )
   {
      std::istringstream dstr( frac.substr(islash,frac.length()-islash) );
      dstr >> denom;
   }
   //while( numer<0 ) numer += DENOM;
   return ((DENOM*numer)/denom) % DENOM;
}

void Operator::parse_shorthand(std::string sh_general, std::string& sh_matrix, int sixtieths[3])
{
   std::ostringstream matrix_str;
   size_t iright = 0;
   size_t ileft = 0;
   size_t iplus = 0;
   // x
   ileft  = iright;
   iright = sh_general.find_first_of(",",ileft);
   iplus  = sh_general.find_first_of("+",ileft);
   if( iplus>iright ) iplus = iright;
   matrix_str << sh_general.substr(ileft,iplus-ileft) << ",";
   sixtieths[0] = frac_to_sixtieths( sh_general.substr(iplus,iright-iplus) );
   // y
   ileft  = iright+1;
   iright = sh_general.find_first_of(",",ileft);
   iplus  = sh_general.find_first_of("+",ileft);
   if( iplus>iright ) iplus = iright;
   matrix_str << sh_general.substr(ileft,iplus-ileft) << ",";
   sixtieths[1] = frac_to_sixtieths( sh_general.substr(iplus,iright-iplus) );
   // z
   ileft  = iright+1;
   iright = sh_general.length();
   iplus  = sh_general.find_first_of("+",ileft);
   if( iplus>iright ) iplus = iright;
   matrix_str << sh_general.substr(ileft,iplus-ileft);
   sixtieths[2] = frac_to_sixtieths( sh_general.substr(iplus,iright-iplus) );
   // result
   sh_matrix = matrix_str.str();
   // debugging
   if( false )
   {
      std::cout << sh_general << " -> " << sh_matrix << " " << sixtieths[0] << " " << sixtieths[1] << " " << sixtieths[2] << std::endl;
   }
}


bool Operator::operator==(const Operator& M) const { bool same=true; for(int i=0; same && i<12; i++) same = (e[i]==M.e[i]); return same; }

bool Operator::operator!=(const Operator& M) const { return !( this->operator==(M) ); }

std::vector<Operator> generate_group(const std::vector<std::string>& generators)
{
   // seed group
   std::vector<Operator> ops;
   ops.push_back( Operator("x,y,z") );
   for(int igen=0; igen<generators.size(); igen++)
      if( generators[igen].length()>4 && generators[igen]!="x,y,z" ) 
         ops.push_back( Operator(generators[igen]) );
   int oldnops = 0;
   int iteration = 0;
   while( oldnops<ops.size() && iteration<10 )
   {
      oldnops = ops.size();
      if( false ) std::cout << "stage " << iteration << " numops=" << oldnops << std::endl;
      for(int i=0; i<oldnops; i++)
      {
         for(int j=0; j<oldnops; j++)
         {
            Operator newop(ops[i],ops[j]);
            std::vector<Operator>::const_iterator iops = std::find(ops.begin(),ops.end(),newop);
            if( iops==ops.end() ) ops.push_back(newop);
         }
      }
      iteration++;

   }
   return ops;
}

std::vector<Operator> get_group(int ipoint)
{
   std::vector<std::string> sh = PointGroupShorthand(ipoint);
   std::vector<Operator> result(sh.size());
   for(int i=0; i<sh.size(); i++) { result[i].assign(sh[i]); }
   return result;
}

std::vector<Operator> get_group(std::string schoenfleis) { return get_group(PointGroup(schoenfleis)); }


}          // namespace Sym

#endif     // SYMMETRY_HPP_ 
