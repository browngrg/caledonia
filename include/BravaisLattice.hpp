#ifndef SYMMETRY_BRAVAIS_H_
#define SYMMETRY_BRAVAIS_H_

#include "Vec.hpp"

namespace Sym {


   // Information about the basis in Bravais Lattices, based on centering
   // Basis vectors are given in halves of the primitive vectors
   template<int C> class Basis {};

   template<> class Basis<PRIMITIVE>
   {
   public:
      enum { NBasis = 1 };
      static void basis(int i, int b[3]) { b[0]=0; b[1]=0; b[2]=0; }
      static std::vector< Vec<int,3> > basis_vectors() 
      { 
         std::vector< Vec<int,3> > result(NBasis); 
         for(int i=0; i<result.size(); i++) 
            basis(i,&(result[i][0])); return result; 
      }
      static int  ibasis(int b[3]) { return 0; }
   };

   template<> class Basis<BODYCENTERED>
   {
   public:
      enum { NBasis = 2 };
      static void basis(int i, int b[3]) 
      { 
         switch(i)
         {
         case  1: b[0]=1; b[1]=1; b[2]=1; break;
         default: b[0]=0; b[1]=0; b[2]=0; 
         }
      }
      static std::vector< Vec<int,3> > basis_vectors() 
      { 
         std::vector< Vec<int,3> > result(NBasis); 
         for(int i=0; i<result.size(); i++) 
            basis(i,&(result[i][0])); return result; 
      }
      static int ibasis(int b[3]) { if( b[0]==0 && b[1]==0 && b[2]==0 ) return 0; return 1; }
   };

   template<> class Basis<FACECENTERED>
   {
   public:
      enum { NBasis = 4 };
      static void basis(int i, int b[3]) 
      { 
         switch(i)
         {
         case  1: b[0]=0; b[1]=1; b[2]=1; break;
         case  2: b[0]=1; b[1]=0; b[2]=1; break;
         case  3: b[0]=1; b[1]=1; b[2]=0; break;
         default: b[0]=0; b[1]=0; b[2]=0; 
         }
      }
      static std::vector< Vec<int,3> > basis_vectors() 
      { 
         std::vector< Vec<int,3> > result(NBasis); 
         for(int i=0; i<result.size(); i++) 
            basis(i,&(result[i][0])); return result; 
      }
      static int ibasis(int b[3]) 
      { 
         if( b[0]!=0 || b[1]!=0 || b[2]!=0 )
         {
            if( b[0]==0 ) return 1;
            if( b[1]==0 ) return 2;
            if( b[2]==0 ) return 3;
         }
      }
   };

   template<> class Basis<BASECENTERED>
   {
   public:
      enum { NBasis = 2 };
      static void basis(int i, int b[3]) 
      { 
         switch(i)
         {
         case  1: b[0]=1; b[1]=1; b[2]=0; break;
         default: b[0]=0; b[1]=0; b[2]=0; 
         }
      }
      static std::vector< Vec<int,3> > basis_vectors() 
      { 
         std::vector< Vec<int,3> > result(NBasis); 
         for(int i=0; i<result.size(); i++) 
            basis(i,&(result[i][0])); return result; 
      }
      static int ibasis(int b[3]) { if( b[0]==0 && b[1]==0 && b[2]==0 ) return 0; return 1; }
   };



   std::vector< Vec<int,3> > BasisVectors(int center=PRIMITIVE)
   {
      switch(center)
      {
      case BODYCENTERED: return Basis<BODYCENTERED>::basis_vectors(); break;
      case FACECENTERED: return Basis<FACECENTERED>::basis_vectors(); break;
      case BASECENTERED: return Basis<BASECENTERED>::basis_vectors(); break;
      case PRIMITIVE:
      default:           return Basis<PRIMITIVE>::basis_vectors();
      }
   }


 
   // The number of possible Bravais Lattices 
   enum { NBraviasLattice = 14 };

   // Returns the most symmetric space group for each Bravais lattice.
   // Have not determined this for all Bravais lattices, so return a
   // very unsymmetric space group (itseqno=1) if unsure.
   int BravaisToITSeqNo(int ibravais)
   {
      int itseqno[15] = {
             -1,   //  Invalid
             -1,   //  <TRICLINIC,PRIMITIVE>            // 1
             -1,   //  <MONOCLINIC,PRIMITIVE>           // 2
             -1,   //  <MONOCLINIC,BASECENTERED>        // 3
             -1,   //  <ORTHORHOMBIC,PRIMITIVE>         // 4
             -1,   //  <ORTHORHOMBIC,BASECENTERED>      // 5
             -1,   //  <ORTHORHOMBIC,BODYCENTERED>      // 6
             -1,   //  <ORTHORHOMBIC,FACECENTERED>      // 7
             -1,   //  <TETRAGONAL,PRIMITIVE>           // 8
             -1,   //  <TETRAGONAL,BODYCENTERED>        // 9
            221,   //  <CUBIC,PRIMITIVE>                // 10
            229,   //  <CUBIC,BODYCENTERED>             // 11
            225,   //  <CUBIC,FACECENTERED>             // 12
            191,   //  <HEXAGONAL,PRIMITIVE>            // 13
             -1    //  <TRIGONAL,PRIMITIVE>             // 14
      };
      int val = -1;
      if( ibravais>0 && ibravais<15 ) val=itseqno[ibravais];
      if( val==-1 )
      {
         std::cout << "Sym::BravaisToITSeqNo(int ibravais): not implemented for ibravais=" << ibravais << std::endl;
         val = 1;
      }
      return val;
   }

   // Run-time conversion of Bravais lattice adjectives to index
   int BravaisIndex(int lattice, int center)
   {
      int ib = 0;
      if( center==PRIMITIVE )
      {
         switch(lattice)
         {
         case TRICLINIC:      ib=1;  break;
         case MONOCLINIC:     ib=2;  break;
         case ORTHORHOMBIC:   ib=4;  break;
         case TETRAGONAL:     ib=8;  break;
         case CUBIC:          ib=10; break;
         case HEXAGONAL:      ib=13; break;
         case TRIGONAL:       ib=14; break;
         }
      }
      if( center==BODYCENTERED )
      {
         switch(lattice)
         {
         case ORTHORHOMBIC:   ib=6;  break;
         case TETRAGONAL:     ib=9;  break;
         case CUBIC:          ib=11; break;
         }
      }
      if( center==FACECENTERED )
      {
         switch(lattice)
         {
         case ORTHORHOMBIC:   ib=7;  break;
         case CUBIC:          ib=12; break;
         }
      }
      if( center==BASECENTERED )
      {
         switch(lattice)
         {
         case MONOCLINIC:     ib=3;  break;
         case ORTHORHOMBIC:   ib=5;  break;
         }
      }
      return ib;
   }

 
   template<int lat, int cent> class BravaisLattice {};


   template<> class BravaisLattice<TRICLINIC,PRIMITIVE>    
   { 
   public:
      enum { INDEX =  1 }; 
      enum { LATTICE = TRICLINIC    }; 
      enum { CENTER = PRIMITIVE    };    
   }; // (1) triclinic P, 



   template<> class BravaisLattice<MONOCLINIC,PRIMITIVE>    
   { 
   public:
      enum { INDEX =  2 }; 
      enum { LATTICE = MONOCLINIC   }; 
      enum { CENTER = PRIMITIVE    };    
   }; // (2) monoclinic P, 



   template<> class BravaisLattice<MONOCLINIC,BASECENTERED> 
   { 
   public:
      enum { INDEX =  3 }; 
      enum { LATTICE = MONOCLINIC   }; 
      enum { CENTER = BASECENTERED };    
   }; // (3) monoclinic C, 



   template<> class BravaisLattice<ORTHORHOMBIC,PRIMITIVE>    
   { 
   public:
      enum { INDEX =  4 }; 
      enum { LATTICE = ORTHORHOMBIC }; 
      enum { CENTER = PRIMITIVE    };    }; 
   // (4) orthorhombic P, 


   template<> class BravaisLattice<ORTHORHOMBIC,BASECENTERED> 
   { 
   public:
      enum { INDEX =  5 }; 
      enum { LATTICE = ORTHORHOMBIC }; 
      enum { CENTER = BASECENTERED };    
   }; // (5) orthorhombic C, 


   
   template<> class BravaisLattice<ORTHORHOMBIC,BODYCENTERED> 
   { 
   public:
      enum { INDEX =  6 }; 
      enum { LATTICE = ORTHORHOMBIC }; 
      enum { CENTER = BODYCENTERED };    
   }; // (6) orthorhombic I,


   
   template<> class BravaisLattice<ORTHORHOMBIC,FACECENTERED> 
   { 
   public:
      enum { INDEX =  7 }; 
      enum { LATTICE = ORTHORHOMBIC }; 
      enum { CENTER = FACECENTERED };    
   }; // (7) orthorhombic F, 


   
   template<> class BravaisLattice<TETRAGONAL,PRIMITIVE>    
   { 
   public:
      enum { INDEX =  8 }; 
      enum { LATTICE = TETRAGONAL   }; 
      enum { CENTER = PRIMITIVE    };    
   }; // (8) tetragonal P, 


   
   template<> class BravaisLattice<TETRAGONAL,BODYCENTERED> 
   { 
   public:
      enum { INDEX =  9 }; 
      enum { LATTICE = TETRAGONAL   }; 
      enum { CENTER = BODYCENTERED };    
   }; // (9) tetragonal I, 


   
   template<> class BravaisLattice<HEXAGONAL,PRIMITIVE>    
   { 
   public:
      enum { INDEX = 13 }; 
      enum { LATTICE = HEXAGONAL    }; 
      enum { CENTER = PRIMITIVE    };    
   }; // (13) hexagonal P, 


   
   template<> class BravaisLattice<TRIGONAL,PRIMITIVE>    
   { 
   public:
      enum { INDEX = 14 }; 
      enum { LATTICE = HEXAGONAL    }; 
      enum { CENTER = PRIMITIVE    };    
   }; // (14) trigonal R,


   
   template<> class BravaisLattice<CUBIC,PRIMITIVE>    
   { 
   public:
      enum { INDEX = 10 }; 
      enum { LATTICE = CUBIC        };
      enum { CENTER = PRIMITIVE    };    
   }; // (10) cubic P, 


   
   template<> class BravaisLattice<CUBIC,BODYCENTERED> 
   { 
   public:
      enum { INDEX   = 11 }; 
      enum { LATTICE = CUBIC }; 
      enum { CENTER  = BODYCENTERED };    
      enum { ITSEQNO = 229 };
      enum { POINT   = SGTraits<ITSEQNO>::POINT };
   }; // (11) cubic I, 


   
   template<> class BravaisLattice<CUBIC,FACECENTERED> 
   { 
      enum { INDEX = 12 }; 
      enum { LATTICE = CUBIC        }; 
      enum { CENTER = FACECENTERED };    
   }; // (12) cubic F,


}     // namespace Sym

#endif  // SYMMETRY_BRAVAIS_H_
