#ifndef GRID_HPP_
#define GRID_HPP_

#include"Symmetry.hpp"
#include"BravaisLattice.hpp"
#include"Vec.hpp"

#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<iostream>


// A grid is a tessellation of a rectangular prism of space
// This class provides information about the topology of the tessellation
//
// The neighbors lists std::vector<int>& nindex produced by the methods named
// neighbors(jsite,displacements,nindex) have the following structure:
//    At the lowest level they are tuples of NBody-1 integers that give the index
//    of the other sites in a cluster. A value of -1 indicates that the cluster extends
//    outside the geometric boundaries of the grid. 
//    A segment of the list consists of the clusters with the jsite (int index) in the grid as the 
//    i-th element of the cluster (following the notation of van de Walle's cluster expansions)
//    There are NBody segments per neighbor list (2 for 2-body, three for 3-body, 4 for four-body)
//
class Grid
{
public:                        // index is the one-dimensional coordinate of the node used for indexing in data arrays
   typedef Vec<int,4> GIndex;  // the 4-dimensional coordinates of the node that specify Bravais cell and ibasis
   typedef Vec<int,3> NrmVec;  // the 3-dimensional coordinates of the node, in "normal coordinates"
                               //         here normal coordinates means specified in terms of the 3 Bravais vectors 
                               //         the "h" stands for half-units, so center of a cell is at (1,1,1) and the
                               //         corners are at (2,0,0), (0,2,0), (0,2,0),...
                               //         These coordinates cover all the possibilities for the Bravais Lattices in 3d
                               //         For the most general space groups would need to do in 60th's. 
                               //         hnorml is a natural system to express displacements and symmetry operations in
public:
   int     L[3];       // The extent of the region in unit cells
   bool    pbc[3];     // Periodic boundary conditions
   int     zbegin;     // The first index of objects belonging to this region
   int     zend;       // The end index of objects belonging to this region
public:
   Grid(int lx=10, int ly=10, int lz=10, int lattice=Sym::CUBIC, int center=Sym::PRIMITIVE);
   void set_bravais(int lattice, int center=Sym::PRIMITIVE);   // Specify the type of Bravais lattice for the grid
   void set_symmetry(int itseqno=221);                         // Specify the symmetry by space group
   int get_lattice()    const { return lsystem; }
   int get_center()     const { return cntring; }
   int get_spacegroup() const { return itseqno; }
   std::string to_string();
   void from_string(std::string str);
public:
   int  size() const { return nbasis*L[0]*L[1]*L[2]; }                      // Number of nodes in the Grid
   int  index(const GIndex& grid_index) const;                              // Map grid coordinate to 1d storage index
   int  index(const NrmVec& nrml_coords) const;                             // Map normal coordinate to storage index
   bool apply_pbc(Grid::NrmVec& nj, int twoL[3]) const;                     // Apply pbc to vector
   void coords(int index, GIndex& grid_index) const;                        // Map storage index to grid coordinates
   void hnorml(int index, NrmVec& nrml) const;                              // Map storage index to normal coordinates
   void hnorml(const GIndex& grid_index, NrmVec& nrml) const;               // Map grid coordinates to normal coordinates
   int  neigbhor(int index, const GIndex& disp) const; 
   int  neighbor(int index, const NrmVec& disp) const;
   void equiv_disp(const NrmVec& disp, std::vector<NrmVec>& alldisp) const;
   void neighbors(int index, const std::vector<NrmVec>& alldisp, std::vector<int>& nindex) const;
   void equiv_three(const Vec<NrmVec,2>& disp, std::vector< Vec<NrmVec,2> >& alldisp) const;                 // basic interface
   void equiv_three(std::vector< Vec<NrmVec,2> >& alldisp) const;                                            // recursive call
   void neighbors(int index, const std::vector< Vec<NrmVec,2> >& alldisp, std::vector<int>& nindex) const;   // list of three-body neighbors
   void equiv_four(const Vec<NrmVec,3>& disp, std::vector< Vec<NrmVec,3> >& alldisp) const;
   void equiv_four(std::vector< Vec<NrmVec,3> >& alldisp) const;                           
   void neighbors(int index, const std::vector< Vec<NrmVec,3> >& alldisp, std::vector<int>& nindex) const;   // list of four-body neighbors
private:
   int                        lsystem;    // Lattice system of tessellation
   int                        cntring;    // Centering of Lattice
   int                        itseqno;    // Description of symmetry by space group (International Table Sequence Number)
   int                        nbasis;     // Number of basis sites
   int                        ibasis[8];  // maps basis vector to ibasis
   std::vector< Vec<int,3> >  basis_vec;  // Vector to basis atoms in unit cell (in half units)
   std::vector<Sym::Operator> sym_op;     // Possible symmetry operations
};


Grid::Grid(int lx, int ly, int lz, int lattice, int center)
{
   L[0] = lx;
   L[1] = ly;
   L[2] = lz;
   pbc[0] = pbc[1]  = pbc[2] = true;
   this->set_bravais(lattice,center);
}

void Grid::set_bravais(int lattice, int center)
{

   lsystem=lattice;
   cntring=center;
   if( lsystem<0 || lsystem>8 ) lsystem=Sym::CUBIC;
   if( cntring<0 || cntring>3 ) cntring=Sym::PRIMITIVE;
   int ibravais = Sym::BravaisIndex(lsystem,cntring);
   this->set_symmetry( Sym::BravaisToITSeqNo(ibravais) );
   if( lsystem!=lattice || cntring!=center )
      std::cout << "Grid::set_bravais(lattice,center): Grid::set_symmetry(itseqno) changed Bravais lattice details" << std::endl;
}

void Grid::set_symmetry(int itseqno_val)
{
   itseqno = itseqno_val;
   lsystem = Sym::LatticeSystem(itseqno);
   cntring = Sym::Centering(itseqno);
   // get the basis vectors for this lattice
   basis_vec = Sym::BasisVectors(cntring);
   nbasis = basis_vec.size();
   sym_op = Sym::get_group( Sym::PointGroup(itseqno) );
   // create map of basis vector to ibasis
   for(int i=0; i<8; i++) ibasis[i] = 0;
   for(int i=0; i<nbasis; i++) 
   {
      Vec<int,3> b = basis_vec[i];
      int itable = b[0]+2*b[1]+4*b[2];
      ibasis[itable] = i;
   }
}

std::string Grid::to_string()
{
   std::ostringstream strm;
   strm << Sym::LatticeSystemString(static_cast<Sym::Adjective>(lsystem)) << " " 
        << Sym::CenteringString(static_cast<Sym::CenteringType>(cntring)) << " " 
        << itseqno;
   std::string str = strm.str();   
   return str;
}

void Grid::from_string(std::string str)
{
   std::string label,system,center;
   std::istringstream sin(str);
   sin >> system >> center;
   lsystem = Sym::ParseLatticeSystem(system);
   cntring = Sym::ParseCenteringType(center);
   set_bravais(lsystem,cntring);
   int spacegroup = -1;
   if( !sin.eof() ) 
   {
      sin >> spacegroup;
   }
   if( spacegroup>=0 ) set_symmetry(spacegroup);
   if( false )
   { 
      std::cout << "\"" << label << "\" \"" << system << "\" \"" << center << "\" \"" << spacegroup << "\"" << std::endl;
      std::cout << system << " " << lsystem << std::endl;
      std::cout << center << " " << cntring << std::endl;
      std::cout << spacegroup << " " << itseqno << std::endl;
   }
}

int Grid::index(const Grid::GIndex& g) const
{
   return g[3]+nbasis*(g[0]+L[0]*(g[1]+L[1]*g[2]));
}

int Grid::index(const Grid::NrmVec& h) const
{
   Vec<int,3> c = h/2;                // cell index
   Vec<int,3> b = h-2*c;              // vector to basis location
   int btable = b[0]+2*(b[1]+2*b[2]); // range = [0,7]
   // reimplement the idex calculation from Grid::index(GIndex)
   // to avoid cost of marshalling data
   return ibasis[btable] + nbasis*(c[0]+L[0]*(c[1]+L[1]*c[2]));
}

void Grid::coords(int index, Grid::GIndex& g) const
{
   int lda0 = nbasis*L[0];
   int lda1 = L[1]*lda0;
   g[2] = index/lda1;
   index -= g[2]*lda1;
   g[1] = index/lda0;
   index -= g[1]*lda0;
   g[0] = index/nbasis;
   g[3] = index - g[0]*nbasis;
}

void Grid::hnorml(int index, Grid::NrmVec& n) const
{
   // get cell coords, shadows(reimplements) Grid::coords(index,GIndex)
   int lda0 = nbasis*L[0];
   int lda1 = L[1]*lda0;
   n[2] = index/lda1;
   index -= n[2]*lda1;
   n[1] = index/lda0;
   index -= n[1]*lda0;
   n[0] = index/nbasis;
   // get basis
   int b = index - n[0]*nbasis;
   // convert to normalized coordinates in units of half cells
   n *= 2;
   n += basis_vec[b];
}

void Grid::hnorml(const Grid::GIndex& g, Grid::NrmVec& n) const
{
   // n = 2*cell + basis
   n[0] = 2*g[0];
   n[1] = 2*g[1];
   n[2] = 2*g[2];
   n += basis_vec[g[3]];
}

int Grid::neigbhor(int index, const Grid::GIndex& disp) const
{
   NrmVec b = basis_vec[disp[3]];
   NrmVec ndisp( disp[0]+b[0], disp[1]+b[1], disp[2]+b[2] );
   return neighbor(index,ndisp);
}

void _grid_pbc(int& i, int L) 
{
   i -= L*((i+(i<0))/L-(i<0));
   // i = i - L*((i+(i<0))/L-(i<0));
   // i = (i<0)? i+(-(i+1)/L+1)*L : i-(i/L)*L;
}

// apply periodic boundary conditions, return false if out of bounds
bool Grid::apply_pbc(Grid::NrmVec& nj, int twoL[3]) const
{
   bool valid = true;
   for(int i=0; i<3; i++) 
   {
      if(pbc[i]) 
         _grid_pbc( nj[i], twoL[i] ); 
      else
         valid = valid && (0<=nj[i]) && (nj[i]<twoL[i]);
   }
   return valid;
}

int Grid::neighbor(int index, const Grid::NrmVec& disp) const
{
   NrmVec ni;
   hnorml(index,ni);
   ni += disp;
   // apply periodic boundary conditions
   bool valid = true;
   for(int i=0; i<3; i++) 
   {
      if(pbc[i]) 
         _grid_pbc( ni[i], 2*L[i] ); 
      else
         valid = valid && (0<=ni[i]) && (ni[i]<=(2*L[i]));
   }
   return (valid)? this->index(ni) : -1; 
}

void Grid::equiv_disp(const Grid::NrmVec& disp, std::vector<Grid::NrmVec>& alldisp) const
{
   std::set<NrmVec> result;
   NrmVec r;
   int nsym = sym_op.size();
   for(int i=0; i<nsym; i++)
   {
      sym_op[i](disp,r);
      result.insert(r);
   }
   alldisp.resize( result.size() );
   std::copy(result.begin(),result.end(),alldisp.begin());
}

void Grid::neighbors(int index, const std::vector<Grid::NrmVec>& alldisp, std::vector<int>& nindex) const
{
   NrmVec ni,nj;
   hnorml(index,ni);
   int twoL[3] = { 2*L[0], 2*L[1], 2*L[2] };
   int numdisp = alldisp.size();
   nindex.resize(2*numdisp);
   for(int j=0; j<numdisp; j++)
   { 
     nj = ni + alldisp[j];
     bool valid = apply_pbc(nj,twoL);
     nindex[j] = (valid)? this->index(nj) : -1; 
     nj = ni - alldisp[j];
     valid = apply_pbc(nj,twoL);
     nindex[numdisp+j] = (valid)? this->index(nj) : -1; 
   }
}

////////// Three-body calls

void Grid::equiv_three(const Vec<NrmVec,2>& disp, std::vector< Vec<NrmVec,2> >& alldisp) const
{
   alldisp.resize(1);
   alldisp[0][0] = disp[0];
   alldisp[0][1] = disp[1];
   equiv_three(alldisp);
}

void Grid::equiv_three(std::vector< Vec<NrmVec,2> >& alldisp) const
{
   int oldsize = alldisp.size();
   std::set< Vec<NrmVec,2> > result;
   Vec<NrmVec,2> r;
   int nsym = sym_op.size();
   for(int j=0; j<oldsize; j++)
   {
      for(int i=0; i<nsym; i++)
      {
         sym_op[i](alldisp[j][0],r[0]);
         sym_op[i](alldisp[j][1],r[1]);
         result.insert(r);
      }
   }
   alldisp.resize( result.size() );
   std::copy(result.begin(),result.end(),alldisp.begin());
// I haven't found a case where this found new displacement patterns
//std::cout << "oldsize=" << oldsize << " newsize=" << alldisp.size() << std::endl;
   if( oldsize<alldisp.size() ) equiv_three(alldisp);
};

// Find neighbors involved in three body interactions
void Grid::neighbors(int index, const std::vector< Vec<Grid::NrmVec,2> >& alldisp, std::vector<int>& nindex) const
{
   NrmVec ni,nj;
   hnorml(index,ni);
   int twoL[3] = { 2*L[0], 2*L[1], 2*L[2] };
   int numdisp = alldisp.size();
   nindex.resize(3*2*numdisp);
   // This is for (0, d0, d1)
   int index_off = 0;
   for(int j=0; j<numdisp; j++)
   { 
      for(int id=0; id<2; id++)
      {
         nj = ni + alldisp[j][id];
         bool valid = apply_pbc(nj,twoL);
         nindex[2*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
   // This is for (-d0,0,d1-d0)
   index_off += 2*numdisp;
   Vec< NrmVec,2 > tdisp;
   for(int j=0; j<numdisp; j++)
   { 
      tdisp[0] =              -alldisp[j][0];
      tdisp[1] = alldisp[j][1]-alldisp[j][0];
      for(int id=0; id<2; id++)
      {
         nj = ni + tdisp[id];
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+2*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
   // This is for (-d1,d0-d1,0)
   index_off += 2*numdisp;
   for(int j=0; j<numdisp; j++)
   { 
      tdisp[0] =              -alldisp[j][1];
      tdisp[1] = alldisp[j][0]-alldisp[j][1];
      for(int id=0; id<2; id++)
      {
         nj = ni + tdisp[id];
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+2*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
}

////////// Four-body calls

void Grid::equiv_four(const Vec<NrmVec,3>& disp, std::vector< Vec<NrmVec,3> >& alldisp) const
{
   alldisp.resize(1);
   alldisp[0][0] = disp[0];
   alldisp[0][1] = disp[1];
   alldisp[0][2] = disp[2];
   equiv_four(alldisp);
}

void Grid::equiv_four(std::vector< Vec<NrmVec,3> >& alldisp) const
{
   int oldsize = alldisp.size();
   std::set< Vec<NrmVec,3> > result;
   Vec<NrmVec,3> r;
   int nsym = sym_op.size();
   for(int j=0; j<oldsize; j++)
   {
      for(int i=0; i<nsym; i++)
      {
         sym_op[i](alldisp[j][0],r[0]);
         sym_op[i](alldisp[j][1],r[1]);
         sym_op[i](alldisp[j][2],r[2]);
         result.insert(r);
      }
   }
   alldisp.resize( result.size() );
   std::copy(result.begin(),result.end(),alldisp.begin());
// I haven't found a case where this found new displacement patterns
//std::cout << "equiv_four oldsize=" << oldsize << " newsize=" << alldisp.size() << std::endl;
   if( oldsize<alldisp.size() ) equiv_four(alldisp);
};

// Find neighbors involved in four body interactions
void Grid::neighbors(int index, const std::vector< Vec<Grid::NrmVec,3> >& alldisp, std::vector<int>& nindex) const
{
   const int NBody  = 4;
   const int NOther = NBody-1;
   NrmVec ni,nj;
   hnorml(index,ni);
   int twoL[3] = { 2*L[0], 2*L[1], 2*L[2] };
   int numdisp = alldisp.size();
   nindex.resize(NBody*NOther*numdisp);
   // This is for (0,d0,d1,d2)
   int index_off = 0;
   for(int j=0; j<numdisp; j++)
   { 
      for(int id=0; id<NOther; id++)                  // 3=NBody-1, sites related to "index"
      {
         nj = ni + alldisp[j][id];
         // apply periodic boundary conditions
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+NOther*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
   // This is for (-d0,0,d1-d0,d2-d0)
   index_off += NOther*numdisp;
   Vec<NrmVec,NOther> tdisp;
   for(int j=0; j<numdisp; j++)
   { 
      tdisp[0] =              -alldisp[j][0];
      tdisp[1] = alldisp[j][1]-alldisp[j][0];
      tdisp[2] = alldisp[j][2]-alldisp[j][0];
      for(int id=0; id<NOther; id++)
      {
         nj = ni + tdisp[id];
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+NOther*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
   // This is for (-d1,d0-d1,0,d2-d1)
   index_off += NOther*numdisp;
   for(int j=0; j<numdisp; j++)
   { 
      tdisp[0] =              -alldisp[j][1];
      tdisp[1] = alldisp[j][0]-alldisp[j][1];
      tdisp[2] = alldisp[j][2]-alldisp[j][1];
      for(int id=0; id<NOther; id++)
      {
         nj = ni + tdisp[id];
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+NOther*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
   // This is for (-d2,d0-d2,d1-d2,0)
   index_off += NOther*numdisp;
   for(int j=0; j<numdisp; j++)
   { 
      tdisp[0] =              -alldisp[j][2];
      tdisp[1] = alldisp[j][0]-alldisp[j][2];
      tdisp[2] = alldisp[j][1]-alldisp[j][2];
      for(int id=0; id<NOther; id++)
      {
         nj = ni + tdisp[id];
         bool valid = apply_pbc(nj,twoL);
         nindex[index_off+NOther*j+id] = (valid)? this->index(nj) : -1; 
      }
   }
}


///// Functions that use Grid


struct GridNghbr
{
public:
   Grid::NrmVec              disp;     // Initial displacement
   int                       d2;       // distance squared
   int                       degen;    // degeneracy
   std::vector<Grid::NrmVec> alldisp;  // equivalent displacements
public:
   GridNghbr() {;}
   GridNghbr(const GridNghbr& gn) { copy(gn); }
   void operator=(const GridNghbr& gn) { copy(gn); }
   void copy(const GridNghbr& gn) { disp=gn.disp; d2=gn.d2; degen=gn.degen; alldisp=gn.alldisp; }
};


// assumes pbc in all directions
void FindNeighbors(const Grid& grid, std::vector<GridNghbr>& nclasses)
{
   std::multimap<int,GridNghbr> results;
   int numnode = grid.size();
   for(int inode=1; inode<numnode; inode++)
   {
      Grid::NrmVec disp;
      grid.hnorml(inode,disp);
      bool unclassed = true;
      for(std::multimap<int,GridNghbr>::const_iterator imap=results.begin(); imap!=results.end(); imap++)
      {
         const GridNghbr& oldrec(imap->second);
         for(int j=0; j<oldrec.alldisp.size() && unclassed; j++)
            if( disp==oldrec.alldisp[j] ) unclassed = false;
      }
      if( unclassed )
      {
         GridNghbr newclass;
         newclass.disp = disp;
         newclass.d2 = disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2];
         grid.equiv_disp(disp,newclass.alldisp);
         newclass.degen = newclass.alldisp.size();
         results.insert( std::pair<int,GridNghbr>(newclass.d2,newclass) );
      }
   }
   nclasses.resize(results.size());
   int iclass = 0;
   for(std::multimap<int,GridNghbr>::const_iterator imap=results.begin(); imap!=results.end(); imap++)
      nclasses[iclass++] = imap->second;
}

// Find the details about the nearest neighbors
GridNghbr FindNearNeighbors(const Grid& grid)
{
  GridNghbr newclass;
  std::vector< Grid::NrmVec > basis = Sym::BasisVectors(grid.get_center());
  newclass.disp = Grid::NrmVec(2,2,2);
  if( basis.size()>2 ) newclass.disp = basis[1]; // basis[0]=(0,0,0)
  newclass.d2 = newclass.disp[0]*newclass.disp[0]+newclass.disp[1]*newclass.disp[1]+newclass.disp[2]*newclass.disp[2];
  grid.equiv_disp(newclass.disp,newclass.alldisp);
  newclass.degen = newclass.alldisp.size();
  return newclass;
}


/*
// Find the centering from the primitive vectors
// numbering is defined in BravaisLattice.hpp
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
      std::cout << "DeduceCenter(p0,p1,p2) only classified " << check << " vectors" << std::endl
                << "p0 = " << p0[0] << " " << p0[1] << " " << p0[2] << " |p0|^2=" << mag0 << std::endl
                << "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << " |p1|^2=" << mag1 << std::endl
                << "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << " |p2|^2=" << mag2 << std::endl
                << "axial = " << axial << " face = " << face << " body = " << body << std::endl;
   }
   if(axial==3) return Sym::PRIMITIVE;
   if(face ==3) return Sym::FACECENTERED;
   if(body ==3) return Sym::BODYCENTERED;
   if(face ==2) return Sym::FACECENTERED;
   std::cout << "DeduceCenter(p0,p1,p2) could not classify" << std::endl
             << "p0 = " << p0[0] << " " << p0[1] << " " << p0[2] << " |p0|^2=" << mag0 << std::endl
             << "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << " |p1|^2=" << mag1 << std::endl
             << "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << " |p2|^2=" << mag2 << std::endl
             << "axial = " << axial << " face = " << face << " body = " << body << std::endl;
}
*/




#endif  // GRID_HPP
