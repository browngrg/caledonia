#ifndef NBODYTABLE_HPP_
#define NBODYTABLE_HPP_

// Data structure for encoding NBody interactions
//
// April 23, 2014
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)

// TODO: Test the build_ave methods, as well as the concept

#include"Vec.hpp"
#include"Grid.hpp"

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>
#include<map>
#include<string>
#include<algorithm>

// "Multicomponent multisublattice alloys, nonconfigurational entropy
// and other additions to the Alloy Theoretic Automated Toolkit,"
// A. van de Walle   http://arxiv.org/abs/0906.1608v1
// June 8, 2009
//
// q(sigma) = sum_a m_a J_a <Gamma_a'(sigma)>_a
//     a (alpha) is a cluster.
//         for a binary alloy, a_i specifies if site i in cluster or not
//         for M_i>2 can take values 0,1,...,M_i-1; 0 == not part of cluster
//     the average <...>_a is over all clusters a' that are symmetrically equivalent
//     Gamma_a'(sigma) are cluster functions
//        Gamma_a(sigma) = prod_i gamma_{ai,Mi}(sigma_i)
//        gamma described in papers. The choice is (apparently for homogeneous M)
//                                              1               if ai=0
//           gamma_{ai,Mi} =         -cos(2pi*ceil(ai/2)*sigma_i/M) if ai>0 and odd
//                                   -sin(2pi*ceil(ai/2)*sigma_i/M) if ai>0 and even
//
//     m_a is the multiplicity (degeneracy), number of symmetry equivalent clusters
//     J_a is an "Effective Cluster Interaction," these are the expansion coefficients
//     
// The number of nonzero sites in the cluster is the degree of the NBody interaction.
//
// Default ATAT choice as described by van de Walle 2009.

// This is the gamma function defined in the reference above  (van de Walle 2009)
// gamma(alpha_i) = clusterfun(alpha_i-1)  [check, old calls may send alpha_i-1]
double gamma_ATAT(int alphai, int M, int sigmai)
{
   // const long double PI = 3.141592653589793238L;
   // const double PI = 3.141592653589793;
   // const float PI = 3.1415927;
   const double pi = 4.*atan(1.);
   const double twopi = 2.*pi;
   if(alphai==0) return 1;
   int rup = static_cast<int>((alphai+1)/2);   
   double x = static_cast<double>(sigmai)/static_cast<double>(M);
   double val = (alphai%2)? cos(twopi*rup*x) : sin(twopi*rup*x);
   if(fabs(val)<1.e-10) val=0;
   if(val>2) std::cout << "gamma wrong " << val << " " << alphai << " " << M << " " << sigmai << std::endl;
   return -val;
}


// This is the clsuster function appropriate to clusters.out
// They are defined in gcdhelp.txt as:
//  The cluster functions in a m-component system are defined as
//   function '0' : -cos(2*PI*  1  *s/m)
//   function '1' : -sin(2*PI*  1  *s/m)
//                             .
//                             .
//                             .
//                  -cos(2*PI*[m/2]*s/m)
//                  -sin(2*PI*[m/2]*s/m)   <--- the last sin( ) is omitted if m is even
//  where the occupation variable s can take any values in {0,...,m-1}
//  and [...] denotes the 'round down' operation.
//  Note that, these functions reduce to the single function (-1)^s in the binary case.
// 
double clusterfun(int funindex, int M, int sigmai)
{
   const double pi = 4.*atan(1.);
   const double twopi = 2.*pi;
   double x = static_cast<double>(sigmai)/static_cast<double>(M);
   double val = 1;
   switch(funindex)
   {
      case 0: val = -cos(twopi*x); break;       // 2pi/1
      case 1: val = -sin(twopi*x); break;
      case 2: val = -cos(   pi*x); break;       // 2pi/2
      case 3: val = -sin(   pi*x); break;
      case 4: val = -cos(twopi/3.*x); break;    // 2pi/3
      case 5: val = -sin(twopi/3.*x); break;
      case 6: val = -cos(   pi/2.*x); break;    // 2pi/4
      case 7: val = -sin(   pi/2.*x); break;
      case 8: val = -cos(twopi/5.*x); break;    // 2pi/5
      case 9: val = -sin(twopi/5.*x); break;
      default: std::cout << "clusterfun not implemented for funindex>9" << std::endl;
   }
   if(fabs(val)<1.e-10) val=0;
   if(val>2) std::cout << "clusterfun" << funindex << " wrong " << val << " ,M=" << M << " sigma=" << sigmai << std::endl;
   return val;
}


// Class NBodyTable is M x M ... = M^N table of interactions J(sigma_p,sigma_j,...) for cluster with N sites.
// N=0 is the vacuum energy (a constant for all confiurations)
// N=1 is the one-site energy (chemical potential, see literature for a better explaination)
// N=2 are pairwise interactions
// N=3 three body interactions, etc.
// Currently each N is a specialization, an arbitrary definition for N>=2 is possible and may emerge with time.
// But is a general definition necessary?
template<int N, int M> struct NBodyTable;

// This changes the calculated energy 5/6/2014 GPB
// Have to keep =false for now to get right energy
static const bool BUILDAVE=false;

// Constant for the whole system
template<int M>
struct NBodyTable<0,M>
{
public:
   typedef Vec<float,3> Vec3;
   float J;
public:
   NBodyTable(float Jval=0) { J=Jval; }
   NBodyTable(const NBodyTable<0,M>& orig) { J=orig.J; }
};

// One site energy J[sigmai]
template<int M>
struct NBodyTable<1,M>
{
public:
//   enum { M = MCFramework::MSpecies };  // Maximum number of species
   float J[M];                          // The table of energies
public:
   NBodyTable(float Jval=0) { for(int i=0; i<M; i++) J[i]=Jval; }
   NBodyTable(const NBodyTable<1,M>& orig) { for(int i=0; i<M; i++) J[i]=orig.J[i]; }
public:
   // Build the table from ATAT data
   // Use a list of length 1 to make interface compatible with N>1
   void build(const std::vector<int>& cfi, float Jeci)
   {
      for(int sigma=0; sigma<M; sigma++)
      {
         double Gamma = clusterfun(cfi[0],M,sigma);
         this->J[sigma] = Jeci*Gamma;
      }
   }
   // Add the contribution from geometrically-related cluster
   // The two clusters should only differ in cfi and Jeci (not an issue for N=1)
   void add(const std::vector<int>& cfi, float Jeci)
   {
      for(int sigma=0; sigma<M; sigma++)
      {
         double Gamma = clusterfun(cfi[0],M,sigma);
         this->J[sigma] += Jeci*Gamma;
      }
   }
};

// Pairwise energy J(sigmai[0],sigmai[1])
template<int M>
struct NBodyTable<2,M>
{
public:
// enum { M = MCFramework::MSpecies };
   typedef Vec<float,3> Vec3;            // Describes displacement in lattice normal units
   Vec3 disp[1];                         // List of displacements from cluster center, for N=1 pair displacement
   float J[M*M];                         // Table of interactions J(sigma0,sigma1) = J[sigma1+M*sigma0]
   int degen;                            // Number of disp related by space-group symmetry defined by grid
   std::vector<Grid::NrmVec> d;          // set of all equivalent displacements
   std::vector<int> nghbr;               // set of neighbor indices, depends on isite. Provided for (re)use by calling algorithms
   int nghbr_isite;                      // site that nghbr is calculated for, external algorithms need to avoid recalculating nghbr
public:
   NBodyTable(Vec3 displace=Vec3(1,0,0),float Jval=1)             // Default constructor is a Potts model interaction
   {
      disp[0] = displace;
      for(int i=0; i<M*M; i++) this->J[i]=0;
      for(int i=0; i<M;   i++) this->J[i+i*M]=Jval;
   }
   NBodyTable(const NBodyTable<2,M>& orig) { this->copy(orig); }   // Copy constructor calls the deep copy method
public:
   // Build one term (of alpha[]={ {0,0}, {1,0}, {1,1} }  as defined by site_energy.pdf
   // These are then combined using the this->add(cluster) method defined below
   // The energy mJeci is the multiplicity times the cluster expansion coefficient mi*Jeci
   // Currently the factor 1/NN is included at calculation in MCMethod::full_calc
   void build(Grid& grid, std::vector<Vec3>& displace, const std::vector<int>& cfi, float mJeci)
   {
      //  The neighbor displacement
      disp[0] = displace[1]-displace[0]; 
      Grid::NrmVec hdisp(2*disp[0][0],2*disp[0][1],2*disp[0][2]);
      grid.equiv_disp(hdisp,this->d);
      degen = this->d.size();
      if(BUILDAVE) { this->build_ave(cfi,mJeci); return; }
      mJeci /= static_cast<float>(degen);               // The factor 1/N from site_energy.pdf
      for(int i=0; i<M*M; i++) this->J[i]=0;
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         double Gamma0 = clusterfun(cfi[0],M,sigma0);
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            double Gamma1 = clusterfun(cfi[1],M,sigma1);
            this->J[sigma1+M*sigma0] = mJeci*Gamma0*Gamma1;
//          this->J[sigma1+M*sigma0] += mJeci*Gamma0*Gamma1/2.;
//          this->J[sigma0+M*sigma1] += mJeci*Gamma0*Gamma1/2.;
         }
      }
   }
   // Build one term assuming that <Gamma_a> = 0.5*( Gamma_(alpha0,alpha1) + Gamma(alhpa1,alpha0) )
   // This gives a symmetric table J
   void build_ave(const std::vector<int>& cfi, float Jeci)
   {
      Jeci /= static_cast<float>(degen);                // The factor 1/N from site_energy.pdf
      for(int i=0; i<M*M; i++) this->J[i] = 0;
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            double Gamma0 = clusterfun(cfi[0],M,sigma0);
            double Gamma1 = clusterfun(cfi[1],M,sigma1);
            this->J[sigma1+M*sigma0] += Jeci*Gamma0*Gamma1;
            Gamma0 = clusterfun(cfi[1],M,sigma0);
            Gamma1 = clusterfun(cfi[0],M,sigma1);
            this->J[sigma1+M*sigma0] += Jeci*Gamma0*Gamma1;
         }
      }
      for(int i=0; i<M*M; i++) this->J[i] /= 2;
   }
   // Add the contribution of another cluster to this one
   // Used to combine terms created by method this->build()
   void add(const NBodyTable<2,M>& clust)
   {
      for(int i=0; i<M*M; i++) this->J[i] += clust.J[i];
   }
   // deep copy of cluster information
   void copy(const NBodyTable<2,M>& clust)
   {
      disp[0] = clust.disp[0];
      degen   = clust.degen;
      for(int i=0; i<M*M; i++) this->J[i] = clust.J[i];
      d.resize(clust.d.size()); std::copy(clust.d.begin(),clust.d.end(),d.begin());
      nghbr.resize(clust.nghbr.size()); std::copy(clust.nghbr.begin(),clust.nghbr.end(),nghbr.begin());    
      nghbr_isite = clust.nghbr_isite;
   }
   bool operator==(const NBodyTable<2,M>& other)
   {
     return (this->disp[0])==(other.disp[0]);
   }
   bool operator!=(const NBodyTable<2,M>& other)
   {
     return (this->disp[0])!=(other.disp[0]);
   }
};


// Partial specialization for three-body interactions
// See comments on two-body term NBodyTable<2> for explanation of methods
template<int M>
struct NBodyTable<3,M>
{
public:
// enum { M = MCFramework::MSpecies };
   typedef Vec<float,3> Vec3;
   Vec3 disp[2];
   float J[M*M*M];                         // J(s0,s1,s2) = J[s2+M*s1+M*M*s0]
   int degen;
   std::vector< Vec<Grid::NrmVec,2> > d;   // set of all equivalent displacements
   std::vector<int> nghbr;                 // set of neighbor indices
   int nghbr_isite;                        // site that nghbr is calculated for
public:
   NBodyTable(Vec3 displace0=Vec3(1,0,0), Vec3 displace1=Vec3(1,1,0), float Jval=1)
   {
      disp[0] = displace0;
      disp[1] = displace1;
      for(int i=0; i<M*M*M; i++) this->J[i]=0;
      for(int i=0; i<M;     i++) this->J[i+M*(i+M*i)]=Jval;
      nghbr_isite = -1;
   }
   NBodyTable(const NBodyTable<3,M>& orig) { this->copy(orig); }
public:
   void build(Grid& grid, std::vector<Vec3>& displace, const std::vector<int>& cfi, float Jeci)
   {
      this->disp[0] = displace[1]-displace[0];
      this->disp[1] = displace[2]-displace[0];
      Vec<Grid::NrmVec,2> hdisp = Vec<Grid::NrmVec,2>( Grid::NrmVec(2*disp[0][0],2*disp[0][1],2*disp[0][2]),
                                                       Grid::NrmVec(2*disp[1][0],2*disp[1][1],2*disp[1][2]));
      grid.equiv_three(hdisp,this->d);
      degen = this->d.size();
      if(BUILDAVE) { this->build_ave(cfi,Jeci); return; }   // Use _ave method to build a symmetric table
      Jeci /= static_cast<float>(degen);                // The factor 1/N from site_energy.pdf
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         double Gamma0 = clusterfun(cfi[0],M,sigma0);
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            double Gamma1 = clusterfun(cfi[1],M,sigma1);
            for(int sigma2=0; sigma2<M; sigma2++)
            {
               double Gamma2 = clusterfun(cfi[2],M,sigma2);
               this->J[sigma2+M*(sigma1+M*sigma0)] = Jeci*Gamma0*Gamma1*Gamma2;
            }
         }
      }
   }
   // Build assuming <Gamma> = 0.5*( Gamma(a0,a1,a2)+Gamma(a2,a1,a0) )
   // which means labelling the sites left to right and right to left
   void build_ave(const std::vector<int>& cfi, float Jeci)
   {
      Jeci /= static_cast<float>(degen);                // The factor 1/N from site_energy.pdf
      for(int i=0; i<M*M*M; i++) this->J[i]=0;
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            for(int sigma2=0; sigma2<M; sigma2++)
            {
               int idx = sigma2+M*(sigma1+M*sigma0);
#if 0
               double Gamma0 = clusterfun(cfi[0],M,sigma0);
               double Gamma1 = clusterfun(cfi[1],M,sigma1);
               double Gamma2 = clusterfun(cfi[2],M,sigma2);
               //this->J[sigma2+M*(sigma1+M*sigma0)] += Jeci*Gamma0*Gamma1*Gamma2;
               this->J[idx] += Jeci*Gamma0*Gamma1*Gamma2;
               Gamma0 = clusterfun(cfi[2],M,sigma0);
               Gamma1 = clusterfun(cfi[1],M,sigma1);
               Gamma2 = clusterfun(cfi[0],M,sigma2);
               //this->J[sigma2+M*(sigma1+M*sigma0)] += Jeci*Gamma0*Gamma1*Gamma2;
               this->J[idx] += Jeci*Gamma0*Gamma1*Gamma2;
               this->J[idx] /= 2;
#endif
               static int perm[6][3] = { {0,1,2},
                                         {0,2,1},
                                         {1,0,2},
                                         {1,2,0},
                                         {2,0,1},
                                         {2,1,0} };
               for(int ip=0; ip<6; ip++)
               {
                  double Gamma0 = clusterfun(cfi[perm[ip][0]],M,sigma0);
                  double Gamma1 = clusterfun(cfi[perm[ip][1]],M,sigma1);
                  double Gamma2 = clusterfun(cfi[perm[ip][2]],M,sigma2);
                  this->J[idx] += Jeci*Gamma0*Gamma1*Gamma2;
               }
               this->J[idx] /= 6.;
            }
         }
      }
   }
   void add(const NBodyTable<3,M>& clust)
   {
      for(int i=0; i<M*M*M; i++) 
      {
         this->J[i] += clust.J[i];
      }
   }
   void copy(const NBodyTable<3,M>& clust)
   {
      disp[0] = clust.disp[0];
      disp[1] = clust.disp[1];
      degen   = clust.degen;
      for(int i=0; i<M*M*M; i++)
      {
         this->J[i] = clust.J[i];
      }
      d.resize(clust.d.size()); std::copy(clust.d.begin(),clust.d.end(),d.begin());
      nghbr.resize(clust.nghbr.size()); std::copy(clust.nghbr.begin(),clust.nghbr.end(),nghbr.begin());    
      nghbr_isite = clust.nghbr_isite;
   }
   bool operator==(const NBodyTable<3,M>& other)
   {
     return ( (this->disp[0])==(other.disp[0]) && (this->disp[1])==(other.disp[1]) );
   }
   bool operator!=(const NBodyTable<3,M>& other)
   {
     return ( (this->disp[0])!=(other.disp[0]) && (this->disp[1])!=(other.disp[1]) );
   }
};

template<int M>
struct NBodyTable<4,M>
{
public:
// enum { M = MCFramework::MSpecies };
   typedef Vec<float,3> Vec3;
   Vec3 disp[3];
   float J[M*M*M*M];
   int degen;
   std::vector< Vec<Grid::NrmVec,3> > d;   // set of all equivalent displacements
   std::vector<int> nghbr;                 // set of neighbor indices
   int nghbr_isite;                        // site that nghbr is calculated for
public:
   NBodyTable(Vec3 displace0=Vec3(1,0,0), Vec3 displace1=Vec3(1,1,0), Vec3 displace2=Vec3(1,1,1), float Jval=1)
   {
      disp[0] = displace0;
      disp[1] = displace1;
      disp[2] = displace2;
      for(int i=0; i<M*M*M*M; i++) this->J[i]=0;
      for(int i=0; i<M;       i++) this->J[i+M*(i+M*(i+M*i))]=Jval;
   }
   NBodyTable(const NBodyTable<4,M>& orig) { this->copy(orig); }
public:
   void build(Grid& grid, std::vector<Vec3>& displace, const std::vector<int>& cfi, float Jeci)
   {
      this->disp[0] = displace[1]-displace[0];
      this->disp[1] = displace[2]-displace[0];
      this->disp[2] = displace[3]-displace[0];
      Vec<Grid::NrmVec,3> hdisp = Vec<Grid::NrmVec,3>( Grid::NrmVec(2*disp[0][0],2*disp[0][1],2*disp[0][2]),
                                                       Grid::NrmVec(2*disp[1][0],2*disp[1][1],2*disp[1][2]),
                                                       Grid::NrmVec(2*disp[2][0],2*disp[2][1],2*disp[2][2]));
      grid.equiv_four(hdisp,this->d);
      degen = this->d.size();
      if(BUILDAVE) { this->build_ave(cfi,Jeci); return; }
      Jeci /= static_cast<float>(degen);                // The factor 1/N from site_energy.pdf
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         double Gamma0 = clusterfun(cfi[0],M,sigma0);
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            double Gamma1 = clusterfun(cfi[1],M,sigma1);
            for(int sigma2=0; sigma2<M; sigma2++)
            {
               double Gamma2 = clusterfun(cfi[2],M,sigma2);
               for(int sigma3=0; sigma3<M; sigma3++)
               {
                  double Gamma3 = clusterfun(cfi[3],M,sigma3);
                  this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] = Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
               }
            }
         }
      }
   }
   // Build assuming <Gamma> = 1/4( Gamma(a0,a1,a2,a3) + Gamma(a1,a0,a3,a2) + ... )
   // These are the permutations of alpha appropriate for the four-body cluster used
   // in the AlNiFe clusters.out expansion.
   // Not at all sure how to do for general space symmetry, this could be general
   // It groups the possible values of alpha as
   // { { {0,0,0,0} },  
   //   { {1,0,0,0}, {0,0,0,1}, {0,1,0,0}, {0,0,1,0} },  
   //   { {1,1,0,0}, {0,0,1,1}, {0,1,0,1}, {1,0,1,0} },
   //   { {0,1,1,0}, {1,0,0,1},
   //   { {1,1,1,0}, {0,1,1,1}, {1,1,0,1}, {1,0,1,1} },
   //   { {1,1,1,1} }
   // Using the first element of each subset, the rest of the subset are the
   // unique subsets generated by applying the set of permutations
   // S[0123], S[1032], S[2301], S[3120]
   // This partitions the 4^2 distinct sets {w,x,y,z}.
   // Everything would be completely different for M>3!!
   // GPB 5/2/2014
   void build_ave(const std::vector<int>& cfi, float Jeci)
   {
      Jeci /= static_cast<float>(degen);                // The factor 1/N from site_energy.pdf
      for(int i=0; i<M*M*M*M; i++) this->J[i]=0;
      for(int sigma0=0; sigma0<M; sigma0++)
      {
         for(int sigma1=0; sigma1<M; sigma1++)
         {
            for(int sigma2=0; sigma2<M; sigma2++)
            {
               for(int sigma3=0; sigma3<M; sigma3++)
               {
#if 0
                  double Gamma0 = clusterfun(cfi[0],M,sigma0);
                  double Gamma1 = clusterfun(cfi[1],M,sigma1);
                  double Gamma2 = clusterfun(cfi[2],M,sigma2);
                  double Gamma3 = clusterfun(cfi[3],M,sigma3);
                  this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] += Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
                  Gamma0 = clusterfun(cfi[1],M,sigma0);
                  Gamma1 = clusterfun(cfi[0],M,sigma1);
                  Gamma2 = clusterfun(cfi[3],M,sigma2);
                  Gamma3 = clusterfun(cfi[2],M,sigma3);
                  this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] += Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
                  Gamma0 = clusterfun(cfi[2],M,sigma0);
                  Gamma1 = clusterfun(cfi[3],M,sigma1);
                  Gamma2 = clusterfun(cfi[0],M,sigma2);
                  Gamma3 = clusterfun(cfi[1],M,sigma3);
                  this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] += Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
                  Gamma0 = clusterfun(cfi[3],M,sigma0);
                  Gamma1 = clusterfun(cfi[1],M,sigma1);
                  Gamma2 = clusterfun(cfi[2],M,sigma2);
                  Gamma3 = clusterfun(cfi[0],M,sigma3);
                  this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] += Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
               }
            }
         }
      }
      for(int i=0; i<M*M*M*M; i++) this->J[i] /= 4;
#endif
                  static int perm[24][4] = { {0,1,2,3},
                                             {0,1,3,2},
                                             {0,2,1,3},
                                             {0,2,3,1},
                                             {0,3,1,2},
                                             {0,3,2,1},
                                             {1,0,2,3},
                                             {1,0,3,2},
                                             {1,2,0,3},
                                             {1,2,3,0},
                                             {1,3,0,2},
                                             {1,3,2,0},
                                             {2,0,1,3},
                                             {2,0,3,1},
                                             {2,1,0,3},
                                             {2,1,3,0},
                                             {2,3,0,1},
                                             {2,3,1,0},
                                             {3,0,1,2},
                                             {3,0,2,1},
                                             {3,1,0,2},
                                             {3,1,2,0},
                                             {3,2,0,1},
                                             {3,2,1,0} };
                  for(int ip=0; ip<24; ip++)
                  {
                     double Gamma0 = clusterfun(cfi[perm[ip][0]],M,sigma0);
                     double Gamma1 = clusterfun(cfi[perm[ip][1]],M,sigma1);
                     double Gamma2 = clusterfun(cfi[perm[ip][2]],M,sigma2);
                     double Gamma3 = clusterfun(cfi[perm[ip][3]],M,sigma3);
                     this->J[sigma3+M*(sigma2+M*(sigma1+M*sigma0))] += Jeci*Gamma0*Gamma1*Gamma2*Gamma3;
                  }
               }
            }
         }
      }
      for(int i=0; i<M*M*M*M; i++) this->J[i] /= 24;
   }
   void add(const NBodyTable<4,M>& clust)
   {
      for(int i=0; i<M*M*M; i++) 
      {
         this->J[i] += clust.J[i];
      }
   }
   void copy(const NBodyTable<4,M>& clust)
   {
      disp[0] = clust.disp[0];
      disp[1] = clust.disp[1];
      disp[2] = clust.disp[2];
      degen   = clust.degen;
      for(int i=0; i<M*M*M*M; i++)
      {
         this->J[i] = clust.J[i];
      }
      d.resize(clust.d.size()); std::copy(clust.d.begin(),clust.d.end(),d.begin());
      nghbr.resize(clust.nghbr.size()); std::copy(clust.nghbr.begin(),clust.nghbr.end(),nghbr.begin());    
      nghbr_isite = clust.nghbr_isite;
   }
   bool operator==(const NBodyTable<4,M>& other)
   {
     return ( (this->disp[0])==(other.disp[0]) 
           && (this->disp[1])==(other.disp[1])
           && (this->disp[2])==(other.disp[2]) );
   }
   bool operator!=(const NBodyTable<4,M>& other)
   {
     return ( (this->disp[0])!=(other.disp[0])
           && (this->disp[1])!=(other.disp[1])
           && (this->disp[2])!=(other.disp[2]) );
   }
};


template<int N, int M>
void merge(std::vector< NBodyTable<N,M> >& nbodyv, const NBodyTable<N,M>& cluster)
{
   bool newc=true;
   int iold=0;
   while(iold<nbodyv.size() && nbodyv[iold]!=cluster) iold++;
   if( iold<nbodyv.size() )
      nbodyv[iold].add(cluster);
   else
      nbodyv.push_back(cluster);
}



#endif   // NBODYTABLE_HPP_
