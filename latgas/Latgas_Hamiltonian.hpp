#ifndef LATGAS_HAMILTONIAN_HPP_
#define LATGAS_HAMILTONIAN_HPP_

#include"NBodyTable.hpp"
#include"Vec.hpp"
#include"Grid.hpp"
#include"Symmetry.hpp"

#include<vector>
#include<algorithm>
#include<string>
#include<cmath>

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>
#include<map>
#include<string>

#undef LATGAS_FULL_CALC_DEBUG

struct Latgas_Observables
{
public:
   enum { MMax = 4 };    // Must be >= Latgas_Hamiltonian::M
   double E,X;           // Required: Energy and the order parameter
   double V;             // The volume of the system (number of lattice sites)
   double ebody[5];      // Energy up to four-body interactions
   double theta[MMax];   // Coverage
public:
   Latgas_Observables() { E=0; X=0; }
   Latgas_Observables(const Latgas_Observables& orig) { this->copy(orig); }
   void copy(const Latgas_Observables& orig)
   {
      V = orig.V;
      E = orig.E;
      X = orig.X;
      for(int i=0; i<5; i++) ebody[i] = orig.ebody[i]; 
      for(int i=0; i<MMax; i++) theta[i] = orig.theta[i]; 
   }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The model class definition: NBody interactions in a lattice-gas
class Latgas_Hamiltonian
{
public:
   
   enum { M = 4 };                         // Number of species 
   enum { KEEPEQUIV = true  };             // Keep list of equivalent displacements beteween calls, memory intensive
   enum { SPINFLIP  = false };             // The type of Monte Carlo move

   typedef Latgas_Observables Observables; 
   typedef std::vector<char>  Config;

   typedef Vec<float,3> Vec3;              // The basic vector type

   // Description of the NBody interactions
   NBodyTable<0,M>                nbody0;    // Zero-body: constant off-set for the energy
   NBodyTable<1,M>                nbody1;    // One-body: the chemical potential of each species
   std::vector< NBodyTable<2,M> > nbody2;    // Two-body: pairwise interactions
   std::vector< NBodyTable<3,M> > nbody3;
   std::vector< NBodyTable<4,M> > nbody4;

public:
  
   Latgas_Hamiltonian(int L=10, std::string model="Potts");

   Latgas_Hamiltonian(std::string filename) { this->Read(filename); }

   std::string header() const;

   // Calculate all macroscopic variables based on this sigma
   void calc_observable(Config& sigma, Latgas_Observables& macro);

   // Change the value of sigma_i at one grid point, update observables
   void change(Config& sigma, Latgas_Observables& macro, int igrid, char new_sigma);

   // Make an attempted Monte Carlo step
   template<typename MCWalker, typename URNG> void mc_step(MCWalker& walker, URNG& urng);

   // Some code will want this, Latgas ignores it
   double step_size;
   // Monte Carlo moves
   template<typename MCWalker, typename URNG> void spin_flip(MCWalker& walker, URNG& urng);
   template<typename MCWalker, typename URNG> void spin_exch(MCWalker& walker, URNG& urng);

public:   

   // Lattice details
   Grid grid;                              // The bravais lattice of nodes 
   float a[3];                             // lattice constants
   float alpha[3];                         // lattice angles alpha,beta,gamma
   std::map<Vec3,std::string> sublattice;  // Designate species allowed on each sublattice
   std::vector<std::string>   species;     // This maps index to species label
   std::vector<Grid::NrmVec>  exch_disp;   // Allowed neighbor relationships for spin-exchange moves

public:

   void clear();
   void ReadATAT();
   void Read(std::string fname="Latgas_ModelIn.txt");
   void Write(std::string fname="Latgas_ModelOut.txt");

   Latgas_Observables full_calc(Config& s);
   Latgas_Observables site_calc(Config& s, int isite);
   Latgas_Observables full_calc_bysite(Config& s);

};  // class Latgas_Hamiltonian



std::string Latgas_Hamiltonian::header() const
{
   std::ostringstream buffer;
   buffer << "# Latgas Hamiltonian" << std::endl;
   buffer << "# D=" << 3 << " L0=" << grid.L[0] << " L1=" << grid.L[1] << " L2=" << grid.L[2] << " N=" << grid.size() << std::endl;
   return buffer.str();
}


void Latgas_Hamiltonian::clear()
{
   sublattice.clear();
   species.resize(0);
   nbody0.J = 0;
   for(int i=0; i<M; i++) nbody1.J[i] = 0;
   nbody2.resize(0);
   nbody3.resize(0);
   nbody4.resize(0);
}

// Read ATAT output to caledonia 
//
// April 24, 2014
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)

// "Multicomponent multisublattice alloys, nonconfigurational entropy
// and other additions to the Alloy Theoretic Automated Toolkit,"
// A. van de Walle   http://arxiv.org/abs/0906.1608v1
// June 8, 2009
void Latgas_Hamiltonian::ReadATAT()
{
   using namespace std;
   this->clear();                     // clear all nbody terms
   grid.pbc[0] = grid.pbc[1] = grid.pbc[2] = true;
   //////////
   // Read lattice information from lat.in
   string buffer;
   Vec3 p[3];                         // primitive unit cell
   ifstream latin("lat.in");
   // Line 1
   latin >> a[0] >> a[1] >> a[2];
   latin >> alpha[0] >> alpha[1] >> alpha[2];
   getline(latin,buffer);
   if(false) cout << "Lattice constants: " << a[0] << " " <<  a[1] << " " <<  a[2] << endl;
   // Line 2 - 4  primitive vectors
   latin >> p[0][0] >> p[0][1] >> p[0][2]; getline(latin,buffer);
   latin >> p[1][0] >> p[1][1] >> p[1][2]; getline(latin,buffer);
   latin >> p[2][0] >> p[2][1] >> p[2][2]; getline(latin,buffer);
   // rest of the file is sublattice occupations
   getline(latin,buffer);
   std::map<Vec3,std::string>& sublattice(sublattice);
   while( !latin.eof() )
   {
      if(false) cout << buffer << endl;;
      istringstream input(buffer);          
      Vec3 basis;
      input >> basis[0] >> basis[1] >> basis[2];
      input >> sublattice[basis];
      if(false) cout << basis[0] << " " << basis[1] << " " << basis[2] << " " << sublattice[basis] << endl;
      getline(latin,buffer);
   }
   int lsystm = Sym::DeduceLatticeSystem(a,alpha);
   int center = Sym::DeduceCentering(p[0],p[1],p[2]);
   grid.set_bravais(lsystm,center);
   if(true) cout << "deduced space group = " << grid.get_spacegroup() << endl;
   if( sublattice.size()<1 ) cout << "lat.in not parsed properly" << endl;
   if( sublattice.size()>1 ) cout << "lat.in violates current assumption of only one sublattice" << endl;
   // change species list
   if(true)                                       // more an anonymous scope than code to turn off
   {
      string s(sublattice[Vec3(0,0,0)]);
      size_t i = s.find_first_not_of(" \t");
      s = s.substr(i);
      while(s.size()>1 && species.size()<10 )
      {
         i = s.find_first_of(",\n",i);
         species.push_back(s.substr(0,i));

         if( s[i]==',')
            s = s.substr(i+1);
         else
            s = "";
         if(false) cout << species.size() << " " <<  species.back() << " " << s << endl;
      }
   } 
   if( M<species.size() ) cout << "lat.in has " << species.size() << " species, code assumes " << M << endl;
   if( M>species.size() ) species.resize(M,"");
   if(true)
   {
      cout << "Species table deduced form lat.in" << endl;
      for(size_t i=0; i<species.size(); i++)
         cout << i << " " << species[i] << endl;
   }
   ///////////////////////
   // Read cluster expasion coefficients J_alpha from eci.out
   vector<double> eci;
   if(true)
   {
      ifstream ein("eci.out");
      while( !ein.eof() )           // a paranoid way to read the file
      {                             // but not thoroughly debugged
         string buffer;
         double J;
         ein >> J;
         getline(ein,buffer);
         if( !ein.eof() )
            eci.push_back(J);
      }
      if(false)
      {
         cout << "Interaction values (" << eci.size() << "):" << endl;
         for(size_t i=0; i<eci.size(); i++)
            cout << eci[i] << endl;
      }
   } 
   ///////////////////////
   // Read cluster expansion from clusters.out, the cluster descriptions parallel those in eci.out
   // The number of nonzero sites in the cluster is the degree of the NBody interaction.
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
   //     m_a is the multiplicity, number of symmetry equivalent clusters
   //     J_a is an "Effective Cluster Interaction," these are the expansion coefficients
   //     
   if( true )                           // provides a local scope, should not be set false
   {
      NBodyTable<2,M> clust2;             // clusters objects that can be reused
      NBodyTable<3,M> clust3;
      NBodyTable<4,M> clust4;
      ifstream fin("clusters.out");
      int icluster = 0;                 // index that allows parallel access of eci.out data
      while( !fin.eof() )
      {
         int mult,nb,Mm2;               // Multiplicity, num sites in cluster, M-2
         float mdist;                   // maximum distance in cluster
         fin >> mult >> mdist >> nb; 
         vector<Vec3> disp(nb);         // displacements of cluster
         vector<int>  cfi(nb);          // index of the cluster function, cfi=alphai-1 for fun_cfi(sigma_i) = gamma_alphai(sigma_i)
         for(size_t i=0; i<nb; i++)
         {
            fin >> disp[i][0] >> disp[i][1] >> disp[i][2] >> Mm2 >> cfi[i];
         }
         if( fin.eof() ) continue;
         switch(nb)
         {
         case 0:
            nbody0.J += eci[icluster];               // the vacuum energy term
            break;
         case 1:
            nbody1.add(cfi,eci[icluster]);           // the single site energies
            break;
         case 2:
            clust2.build(grid,disp,cfi,mult*eci[icluster]);    // build a cluster description
            merge(nbody2,clust2);                              // merge with known clusters
            break; 
         case 3:
            clust3.build(grid,disp,cfi,mult*eci[icluster]);
            merge(nbody3,clust3);
            break; 
         case 4:
            clust4.build(grid,disp,cfi,mult*eci[icluster]);
            merge(nbody4,clust4);
            break; 
         default:
            cout << "NBody cannot handle " << nb << " body interactions" << endl;
         }
         if(false) cout << icluster << " " << mult << " " << nb << " " << eci[icluster] << endl;
         icluster++;
      }
   }
   if(false)                                         // very verbose output
   {
      cout << "nbody0 = " << nbody0.J << endl;
      cout << "nbody1 = "; 
      for(int i=0; i<M; i++) cout << " " << nbody1.J[i];
      cout << endl;
      for(int icluster=0; icluster<nbody2.size(); icluster++)
      {
         Vec3 d = nbody2[icluster].disp[0];
         cout << endl << "nbody2, cluster = " << icluster;
         cout << " disp= " << d[0] << " " << d[1] << " " << d[2] << endl;
         for(int i=0; i<M; i++)
         {
            for(int j=0; j<M; j++) cout << "   " << nbody2[icluster].J[j+M*i];
            cout << endl;
         }
      }
      for(int icluster=0; icluster<nbody3.size(); icluster++)
      {
         Vec3 d0 = nbody3[icluster].disp[0];
         Vec3 d1 = nbody3[icluster].disp[1];
         cout << endl << "nbody3, cluster = " << icluster;
         cout << " disp0= " << d0[0] << " " << d0[1] << " " << d0[2];
         cout << " disp1= " << d1[0] << " " << d1[1] << " " << d1[2] << endl;
         for(int i=0; i<M; i++)
         {
            for(int j=0; j<M; j++)
            {
               for(int k=0; k<M; k++) cout << "   " << nbody3[icluster].J[k+M*(j+M*i)];
               cout << "\t\t";
            }
            cout << endl;
         }
      }
      for(int icluster=0; icluster<nbody3.size(); icluster++)
      {
         Vec3 d0 = nbody4[icluster].disp[0];
         Vec3 d1 = nbody4[icluster].disp[1];
         Vec3 d2 = nbody4[icluster].disp[2];
         cout << endl << "nbody4, cluster = " << icluster;
         cout << " disp0= " << d0[0] << " " << d0[1] << " " << d0[2];
         cout << " disp1= " << d1[0] << " " << d1[1] << " " << d1[2];
         cout << " disp2= " << d2[0] << " " << d2[1] << " " << d2[2] << endl;
      }
   }
}


void Latgas_Hamiltonian::Write(std::string fname)
{
   using namespace std;
   ofstream fout(fname.c_str());
   const int version = 1;
   if( version==1 )
   { 
      fout << "Latgas_Model 1" << endl;
      fout << "grid " << grid.to_string() << endl;
      fout << "a " << a[0] << " " << a[1] << " " << a[2] << endl;
      fout << "alpha " << alpha[0] << " " << alpha[1] << " " << alpha[2] << endl;
      fout << "L " << grid.L[0] << " " << grid.L[1] << " " << grid.L[2] << endl;
      fout << "pbc " << grid.pbc[0] << " " << grid.pbc[1] << " " << grid.pbc[2] << endl;
      fout << "M " << M;
      for(size_t i=0; i<species.size(); i++) fout << " " << species[i];
      fout << endl;
      fout << "nbody0 " << nbody0.J << endl;
      fout << "nbody1"; 
      for(int i=0; i<M; i++) fout << " " << nbody1.J[i];
      fout << endl;
      fout << "nbody2 " << nbody2.size() << endl;
      for(int icluster=0; icluster<nbody2.size(); icluster++)
      {
         Vec3 d = nbody2[icluster].disp[0];
         fout << d[0] << " " << d[1] << " " << d[2]<< endl;
         for(int i=0; i<M*M; i++) fout << " " << nbody2[icluster].J[i];
         fout << endl;
      }
      fout << "nbody3 " << nbody3.size() << endl;
      for(int icluster=0; icluster<nbody3.size(); icluster++)
      {
         Vec3 d0 = nbody3[icluster].disp[0];
         Vec3 d1 = nbody3[icluster].disp[1];
         fout << d0[0] << " " << d0[1] << " " << d0[2]<< endl;
         fout << d1[0] << " " << d1[1] << " " << d1[2]<< endl;
         for(int i=0; i<M*M*M; i++) fout << " " << nbody3[icluster].J[i];
         fout << endl;
      }
      fout << "nbody4 " << nbody4.size() << endl;
      for(int icluster=0; icluster<nbody4.size(); icluster++)
      {
         Vec3 d0 = nbody4[icluster].disp[0];
         Vec3 d1 = nbody4[icluster].disp[1];
         Vec3 d2 = nbody4[icluster].disp[2];
         fout << d0[0] << " " << d0[1] << " " << d0[2]<< endl;
         fout << d1[0] << " " << d1[1] << " " << d1[2]<< endl;
         fout << d2[0] << " " << d2[1] << " " << d2[2]<< endl;
         for(int i=0; i<M*M*M*M; i++) fout << " " << nbody4[icluster].J[i];
         fout << endl;
      }
   }
 
}
    
void Latgas_Hamiltonian::Read(std::string fname)
{
// Something set inside ATAT not set here that is causing a seg fault
   using namespace std;
   this->clear();
   ifstream fin(fname.c_str());
   string label,buffer;
   // Get version number
   int version = 1;
   fin >> label >> version; getline(fin,buffer); 
   if( label.compare("Latgas_Model")!=0 ) cout << "Latgas_Hamiltonian::Read() unrecognized first line \"" << label << " " << version << "\" in file \"" << fname << "\"" << endl;
   fin >> label; getline(fin,buffer);  grid.from_string(buffer);
   if( version==1 )
   { 
      int MSpecies,  NTwo, NThree, NFour;
      fin >> label >> a[0] >> a[1] >> a[2]; getline(fin,buffer);
      fin >> label >> alpha[0] >> alpha[1] >> alpha[2]; getline(fin,buffer);
      fin >> label >> grid.L[0] >> grid.L[1] >> grid.L[2]; getline(fin,buffer);
      fin >> label >> grid.pbc[0] >> grid.pbc[1] >> grid.pbc[2]; getline(fin,buffer);
      fin >> label >> MSpecies;
      if(MSpecies>M) { MSpecies=M; cout << "Latgas_Hamiltonian::Read() compiled for " << M << " species not " << MSpecies << endl; }
      species.resize(MSpecies);
      for(size_t i=0; i<species.size(); i++) fin >> species[i];
      getline(fin,buffer);
      fin >> label >> nbody0.J; getline(fin,buffer);
      fin >> label;
      for(int i=0; i<MSpecies; i++) fin >> nbody1.J[i];
      getline(fin,buffer);
      fin >> label >> NTwo; getline(fin,buffer);
      nbody2.resize(NTwo);
      for(int icluster=0; icluster<nbody2.size(); icluster++)
      {
         Vec3& d(nbody2[icluster].disp[0]);
         fin >> d[0] >> d[1] >> d[2]; 
         getline(fin,buffer);
	 std::vector<double> Jin(MSpecies*MSpecies);
         for(int i=0; i<Jin.size(); i++) fin >> Jin[i];
         getline(fin,buffer);
	 for(int i=0; i<M*M; i++) nbody2[icluster].J[i] = 0;
         for(int i=0; i<MSpecies; i++)
            for(int j=0; j<MSpecies; j++)
	       nbody2[icluster].J[j+M*i] = Jin[j+MSpecies*i];
         Grid::NrmVec hdisp(2*d[0],2*d[1],2*d[2]);
         grid.equiv_disp(hdisp,nbody2[icluster].d);
         nbody2[icluster].degen = nbody2[icluster].d.size();
      }
      fin >> label >> NThree; getline(fin,buffer);
      nbody3.resize(NThree);
      for(int icluster=0; icluster<nbody3.size(); icluster++)
      {
         Vec3& d0(nbody3[icluster].disp[0]);
         Vec3& d1(nbody3[icluster].disp[1]);
         fin >> d0[0] >> d0[1] >> d0[2]; getline(fin,buffer);
         fin >> d1[0] >> d1[1] >> d1[2]; getline(fin,buffer);
	 int M3 = MSpecies*MSpecies*MSpecies;
	 std::vector<double> Jin(M3);
         for(int i=0; i<M3; i++) fin >> Jin[i];
         getline(fin,buffer);
         for(int i=0; i<M*M*M; i++) nbody3[icluster].J[i] = 0;
         for(int i=0; i<MSpecies; i++)
            for(int j=0; j<MSpecies; j++)
               for(int k=0; k<MSpecies; k++)
                  nbody3[icluster].J[k+M*(j+M*i)] = Jin[k+MSpecies*(j+MSpecies*i)];
         Vec<Grid::NrmVec,2> hdisp = Vec<Grid::NrmVec,2>( Grid::NrmVec(2*d0[0],2*d0[1],2*d0[2]),
                                                          Grid::NrmVec(2*d1[0],2*d1[1],2*d1[2]));
         grid.equiv_three(hdisp,nbody3[icluster].d);
         nbody3[icluster].degen = nbody3[icluster].d.size();
      }
      fin >> label >> NFour; getline(fin,buffer);
      nbody4.resize(NFour);
      for(int icluster=0; icluster<nbody4.size(); icluster++)
      {
         Vec3& d0(nbody4[icluster].disp[0]);
         Vec3& d1(nbody4[icluster].disp[1]);
         Vec3& d2(nbody4[icluster].disp[2]);
         fin >> d0[0] >> d0[1] >> d0[2]; getline(fin,buffer);
         fin >> d1[0] >> d1[1] >> d1[2]; getline(fin,buffer);
         fin >> d2[0] >> d2[1] >> d2[2]; getline(fin,buffer);
	 int M4 = MSpecies*MSpecies*MSpecies*MSpecies;
	 std::vector<double> Jin(M4);
         for(int i=0; i<M4; i++) fin >> Jin[i];
         getline(fin,buffer);
         for(int i=0; i<M4; i++) nbody4[icluster].J[i] = 0;
	 for(int i=0; i<MSpecies; i++)
            for(int j=0; j<MSpecies; j++)
               for(int k=0; k<MSpecies; k++)
                  for(int l=0; l<MSpecies; l++)
                     nbody4[icluster].J[l+M*(k+M*(j+M*i))] = Jin[l+MSpecies*(k+MSpecies*(j+MSpecies*i))];
         Vec<Grid::NrmVec,3> hdisp = Vec<Grid::NrmVec,3>( Grid::NrmVec(2*d0[0],2*d0[1],2*d0[2]),
                                                          Grid::NrmVec(2*d1[0],2*d1[1],2*d1[2]),
                                                          Grid::NrmVec(2*d2[0],2*d2[1],2*d2[2]));
         grid.equiv_four(hdisp,nbody4[icluster].d);
         nbody4[icluster].degen = nbody4[icluster].d.size();
      }
      return;
   }
   cout << "Latgas_Hamiltonian::Read() unrecognized version " << version << " in file \"" << fname << "\"" << endl;
}
    


Latgas_Hamiltonian::Latgas_Hamiltonian(int Lval, std::string model)
{
   species.resize(M);
   for(int i=0; i<M; i++) species[i] = static_cast<char>('A'+i);
   grid.L[0] = Lval;
   grid.L[1] = Lval;
   grid.L[2] = 1;
   grid.set_bravais(Sym::CUBIC,Sym::PRIMITIVE);
   bool set = false;
   if( model.compare("Potts")!=0 )
   {
      set = true;
      nbody0.J = 0;
      for(int i=0; i<M; i++) nbody1.J[i] = 0;
      nbody2.resize(1);
      NBodyTable<2,M>& nn(nbody2[0]);
      nn.disp[0] = Vec3(1.,1.,1.);
      for(int i=0; i<M*M; i++) nn.J[i] = 0;
      for(int i=0; i<M;   i++) nn.J[i+i*M] = 1;
   }
   nbody3.resize(0);
   nbody4.resize(0);
}

// calculate energy and occupation numbers from scratch
Latgas_Observables Latgas_Hamiltonian::full_calc(Latgas_Hamiltonian::Config& s)
{
#undef LATAG_FULL_CALC_DEBUG
#undef MCMODEL_RECALC_EQUIV
   int nsite = s.size();
   // debugging info
#  ifdef LATAG_FULL_CALC_DEBUG
   std::vector<float> site_energy(nsite,nbody0.J);
   std::vector<float> bysite_energy(nsite,nbody0.J); 
   std::vector<float> bysite_energy2(nsite,nbody0.J); 
#  endif
   // single-site calculations
   Latgas_Observables macro;
   for(int i=0; i<M; i++) macro.theta[i] = 0;
   double e1 = 0;
   for(int isite=0; isite<nsite; isite++) 
   {
      macro.theta[s[isite]] += 1;
      e1 += nbody1.J[s[isite]];
#     ifdef LATAG_FULL_CALC_DEBUG
      site_energy[isite] += nbody1.J[s[isite]];
      bysite_energy[isite] += nbody1.J[s[isite]];
      bysite_energy2[isite] += nbody1.J[s[isite]];
#     endif
   }
   if(false) std::cout << "full_calc e1=" << e1 << std::endl;
   // two-site calculations
   double e2 = 0;
   for(int icluster=0; icluster<nbody2.size(); icluster++)
   {
      NBodyTable<2,M>& clust(nbody2[icluster]);
#     ifdef MCMODEL_RECALC_EQUIV
      Grid::NrmVec hdisp(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]);
      grid.equiv_disp(hdisp,clust.d);
#     endif
      for(int isite=0; isite<nsite; isite++)
      {
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
         int numnghbr = clust.degen;
         for(int j=0; j<numnghbr; j++)
         {
            int jsite= clust.nghbr[j];
            if( jsite!=-1 )                                         // This test not needed if known PBC=true
            {
               e2 += clust.J[s[jsite]+M*s[isite]];                 // Might be able to use #Define to speed up execution
#              ifdef LATAG_FULL_CALC_DEBUG
               site_energy[isite] += clust.J[s[jsite]+M*s[isite]];       
               bysite_energy[isite] += clust.J[s[jsite]+M*s[isite]];       
               bysite_energy[jsite] += clust.J[s[jsite]+M*s[isite]];       
               bysite_energy2[isite] += clust.J[s[jsite]+M*s[isite]]/2.; 
               bysite_energy2[jsite] += clust.J[s[jsite]+M*s[isite]]/2.;       
#              endif
            }
         }
      }
      clust.nghbr_isite = nsite-1;
   }
   if(false) std::cout << "full_calc e2=" << e2 << std::endl;
   // three-site calculations
   double e3 = 0;
   for(int icluster=0; icluster<nbody3.size(); icluster++)
   {
      NBodyTable<3,M>& clust(nbody3[icluster]);
#     ifdef MCMODEL_RECALC_EQUIV
      Vec<Grid::NrmVec,2> hdisp = Vec<Grid::NrmVec,2>( Grid::NrmVec(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]),
                                                       Grid::NrmVec(2*clust.disp[1][0],2*clust.disp[1][1],2*clust.disp[1][2]));
      grid.equiv_three(hdisp,clust.d);
#     endif
      for(int isite=0; isite<nsite; isite++)
      {
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
         int numnghbr = 2*clust.degen;                             // Two other atoms per cluster
         for(int j=0; j<numnghbr; j+=2)
         {
            int jsite = clust.nghbr[j];
            int ksite = clust.nghbr[j+1];
            if( jsite!=-1 && ksite!=-1 )
            {
                e3 += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
#               ifdef LATAG_FULL_CALC_DEBUG
                site_energy[isite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
                bysite_energy[isite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
                bysite_energy[jsite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
                bysite_energy[ksite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
                bysite_energy2[isite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])]/3.;
                bysite_energy2[jsite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])]/3.;
                bysite_energy2[ksite] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])]/3.;
#               endif
            }
         }
      }
      clust.nghbr_isite = nsite-1;
   }
   if(false) std::cout << "full_calc e3=" << e3 << std::endl;
   // four-site calculations
   double e4 = 0;
   for(int icluster=0; icluster<nbody4.size(); icluster++)
   {
      NBodyTable<4,M>& clust(nbody4[icluster]);
#     ifdef MCMODEL_RECALC_EQUIV
      Vec<Grid::NrmVec,3> hdisp = Vec<Grid::NrmVec,3>( Grid::NrmVec(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]),
                                                       Grid::NrmVec(2*clust.disp[1][0],2*clust.disp[1][1],2*clust.disp[1][2]),
                                                       Grid::NrmVec(2*clust.disp[2][0],2*clust.disp[2][1],2*clust.disp[2][2]));
      grid.equiv_four(hdisp,clust.d);
#     endif
      for(int isite=0; isite<nsite; isite++)
      {
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
         int numnghbr = 3*clust.degen;                             // Three other atoms per cluster
         for(int j=0; j<numnghbr; j+=3)
         {
            int jsite = clust.nghbr[j];
            int ksite = clust.nghbr[j+1];
            int lsite = clust.nghbr[j+2];
            if( jsite!=-1 && ksite!=-1 && lsite!=-1 )
            {
                e4 += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
#               ifdef LATAG_FULL_CALC_DEBUG
                site_energy[isite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
                bysite_energy[isite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
                bysite_energy[jsite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
                bysite_energy[ksite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
                bysite_energy[lsite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
                bysite_energy2[isite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))]/4.;
                bysite_energy2[jsite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))]/4.;
                bysite_energy2[ksite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))]/4.;
                bysite_energy2[lsite] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))]/4.;
#               endif
            }
         }
      }
      clust.nghbr_isite = nsite-1;
   }
   if(false) std::cout << "full_calc e4=" << e4 << std::endl;
   //
   macro.E = nsite*nbody0.J + e1 + e2 + e3 + e4;
   macro.ebody[0] = nsite*nbody0.J;
   macro.ebody[1] = e1;
   macro.ebody[2] = e2;
   macro.ebody[3] = e3;
   macro.ebody[4] = e4;
   if(false)
   {
      std::cout << "Latgas_Hamiltonian::full_calc results:" << std::endl;
      for(int i=0; i<M; i++)
         std::cout << "Num " << species[i] << " = " << macro.theta[i] << std::endl;
      std::cout << "E NBody 0  = " << nsite*nbody0.J << std::endl; 
      std::cout << "E NBody 1  = " << e1 << std::endl; 
      std::cout << "E NBody 2  = " << e2 << std::endl; 
      std::cout << "E NBody 3  = " << e3 << std::endl; 
      std::cout << "E NBody 4  = " << e4 << std::endl; 
      std::cout << "macro.E    = " << macro.E << "  " << macro.E/nsite << std::endl;
#     ifdef LATAG_FULL_CALC_DEBUG
      float tot = 0; for(int i=0; i<nsite; i++) tot += site_energy[i];
      float tot2 = 0; for(int i=0; i<nsite; i++) tot2 += bysite_energy2[i];
      std::cout << "sum_site   = " << tot     << "  " << tot/nsite << std::endl;
      std::cout << "sum_bysite = " << tot2    << "  " << tot2/nsite << std::endl;
      std::ofstream sout("site_energy.out");
      sout <<"site, sigma_site, E_site, E_site/nsite, E_bysite2, E_bysite" << std::endl;
      for(int i=0; i<nsite; i++) sout << i << " " << static_cast<int>(s[i]) << " " << site_energy[i] << " " << site_energy[i]/nsite  << " " << bysite_energy2[i] << " " << bysite_energy[i] << std::endl;
#     endif
   }
   const bool test_bysite = false;
   if(test_bysite)
   {
      this->full_calc_bysite(s);
   } 
   return macro;
}


// calculate energy and occupation numbers from scratch
Latgas_Observables Latgas_Hamiltonian::site_calc(Latgas_Hamiltonian::Config& s, int isite)
{
   Latgas_Observables macro;
   int nsite = s.size();
   // vacuum energy
   macro.ebody[0] = nbody0.J;
   // single-site calculations
   macro.ebody[1] = nbody1.J[s[isite]];
   // two-site calculations
   macro.ebody[2] = 0;
   for(int icluster=0; icluster<nbody2.size(); icluster++)
   {
      NBodyTable<2,M>& clust(nbody2[icluster]);
      if(!KEEPEQUIV) 
      {
         Grid::NrmVec hdisp(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]);
         grid.equiv_disp(hdisp,clust.d);
      }
      if( clust.nghbr_isite != isite )
      {
         clust.nghbr_isite = isite;
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
      }
      // This gives the right answer in full_calc
      //       bysite_energy2[isite] += clust.J[s[jsite]+M*s[isite]]/2.; 
      //       bysite_energy2[jsite] += clust.J[s[jsite]+M*s[isite]]/2.;       
      int numnghbr = clust.degen;
      for(int j=0; j<numnghbr; j++)
      {
         int jsite= clust.nghbr[j];
         if( jsite!=-1 ) macro.ebody[2] += clust.J[s[jsite]+M*s[isite]];
      }
      for(int j=numnghbr; j<2*numnghbr; j++)
      {
         int jsite= clust.nghbr[j];
         if( jsite!=-1 ) macro.ebody[2] += clust.J[s[isite]+M*s[jsite]];
      }
   }
   // three-site calculations
   macro.ebody[3] = 0;
   for(int icluster=0; icluster<nbody3.size(); icluster++)
   {
      NBodyTable<3,M>& clust(nbody3[icluster]);
      if(!KEEPEQUIV) 
      {
         Vec<Grid::NrmVec,2> hdisp = Vec<Grid::NrmVec,2>( Grid::NrmVec(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]),
                                                          Grid::NrmVec(2*clust.disp[1][0],2*clust.disp[1][1],2*clust.disp[1][2]));
         grid.equiv_three(hdisp,clust.d);
      }
      if( clust.nghbr_isite != isite )
      {
         clust.nghbr_isite = isite;
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
      }
      int numnghbr = 2*clust.degen;
      for(int j=0; j<numnghbr; j+=2)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         if( jsite!=-1 && ksite!=-1 )
            macro.ebody[3] += clust.J[s[ksite]+M*(s[jsite]+M*s[isite])];
      }
      for(int j=numnghbr; j<2*numnghbr; j+=2)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         if( jsite!=-1 && ksite!=-1 )
            macro.ebody[3] += clust.J[s[ksite]+M*(s[isite]+M*s[jsite])];
      }
      for(int j=2*numnghbr; j<3*numnghbr; j+=2)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         if( jsite!=-1 && ksite!=-1 )
            macro.ebody[3] += clust.J[s[isite]+M*(s[ksite]+M*s[jsite])];
      }
   }
   // four-site calculations
   macro.ebody[4] = 0;
   for(int icluster=0; icluster<nbody3.size(); icluster++)
   {
      NBodyTable<4,M>& clust(nbody4[icluster]);
      if(!KEEPEQUIV) 
      {
         Vec<Grid::NrmVec,3> hdisp = Vec<Grid::NrmVec,3>( Grid::NrmVec(2*clust.disp[0][0],2*clust.disp[0][1],2*clust.disp[0][2]),
                                                          Grid::NrmVec(2*clust.disp[1][0],2*clust.disp[1][1],2*clust.disp[1][2]),
                                                          Grid::NrmVec(2*clust.disp[2][0],2*clust.disp[2][1],2*clust.disp[2][2]));
         grid.equiv_four(hdisp,clust.d);
      }
      if( clust.nghbr_isite != isite )
      {
         clust.nghbr_isite = isite;
         grid.neighbors(isite,clust.d,clust.nghbr);   // generate neighbor list
      }
      int numnghbr = 3*clust.degen;
      for(int j=0; j<numnghbr; j+=3)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         int lsite= clust.nghbr[j+2];
         if( jsite!=-1 && ksite!=-1 && lsite!=-1 )
            macro.ebody[4] += clust.J[s[lsite]+M*(s[ksite]+M*(s[jsite]+M*s[isite]))];
      }
      for(int j=numnghbr; j<2*numnghbr; j+=3)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         int lsite= clust.nghbr[j+2];
         if( jsite!=-1 && ksite!=-1 && lsite!=-1 )
            macro.ebody[4] += clust.J[s[lsite]+M*(s[ksite]+M*(s[isite]+M*s[jsite]))];
      }
      for(int j=2*numnghbr; j<3*numnghbr; j+=3)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         int lsite= clust.nghbr[j+2];
         if( jsite!=-1 && ksite!=-1 && lsite!=-1 )
            macro.ebody[4] += clust.J[s[lsite]+M*(s[isite]+M*(s[ksite]+M*s[jsite]))];
      }
      for(int j=3*numnghbr; j<4*numnghbr; j+=3)
      {
         int jsite= clust.nghbr[j];
         int ksite= clust.nghbr[j+1];
         int lsite= clust.nghbr[j+2];
         if( jsite!=-1 && ksite!=-1 && lsite!=-1 )
            macro.ebody[4] += clust.J[s[isite]+M*(s[lsite]+M*(s[ksite]+M*s[jsite]))];
      }
   }
   //
   macro.E = macro.ebody[0];
   for(int i=1; i<5; i++) macro.E += macro.ebody[i];
   return macro;
}
    

Latgas_Observables Latgas_Hamiltonian::full_calc_bysite(Latgas_Hamiltonian::Config& s)
{
   Latgas_Observables macro;
   for(int i=0; i<M; i++) macro.theta[i] = 0;
   for(int i=0; i<5; i++) macro.ebody[i] = 0;
   int nsite = s.size();
#  ifdef LATAG_FULL_CALC_DEBUG
   std::vector<float> site_energy(nsite,0);
#  endif
   for(int isite=0; isite<nsite; isite++)
   {
      macro.theta[s[isite]] += 1;
      Latgas_Observables macro_site = site_calc(s,isite);
      for(int i=0; i<5; i++) macro.ebody[i] += macro_site.ebody[i];
#     ifdef LATAG_FULL_CALC_DEBUG
      site_energy[isite] = macro_site.E;
#     endif
   }
   macro.E = macro.ebody[0];
   for(int i=1; i<5; i++) macro.E += macro.ebody[i]/static_cast<float>(i);
   if(false)
   {
      std::cout << __FILE__ << ":" << __LINE__ << " full_calc_bysite results:" << std::endl;
      for(int i=0; i<M; i++)
         std::cout << "Num " << species[i] << " = " << macro.theta[i] << std::endl;
      std::cout << "E NBody 0 = " << macro.ebody[0] << " " << macro.ebody[0]/1. << std::endl; 
      std::cout << "E NBody 1 = " << macro.ebody[1] << " " << macro.ebody[1]/1. << std::endl; 
      std::cout << "E NBody 2 = " << macro.ebody[2] << " " << macro.ebody[2]/2. << std::endl; 
      std::cout << "E NBody 3 = " << macro.ebody[3] << " " << macro.ebody[3]/3. << std::endl; 
      std::cout << "E NBody 4 = " << macro.ebody[4] << " " << macro.ebody[4]/4. << std::endl; 
      std::cout << "macro.E  = " << macro.E << "  " << macro.E/nsite << std::endl;
#     ifdef LATAG_FULL_CALC_DEBUG
      std::ofstream sout("site_energy_bysite.out");
      sout <<"site, sigma_site, E_site, E_site/nsite" << std::endl;
      for(int i=0; i<nsite; i++) sout << i << " " << static_cast<int>(s[i]) << " " << site_energy[i] << " " << site_energy[i]/nsite << std::endl;
#     endif
   }
   return macro;
}

// Calculate all macroscopic variables based on this sigma
void Latgas_Hamiltonian::calc_observable(Latgas_Hamiltonian::Config& s, Latgas_Observables& macro)
{
   for(int i=0; i<M; i++) macro.theta[i] = 0;
   for(int i=0; i<5; i++) macro.ebody[i] = 0;
   int nsite = s.size();
   macro.V = nsite;
#  ifdef LATGAS_FULL_CALC_DEBUG
   std::vector<float> site_energy(nsite,0);
#  endif
   for(int isite=0; isite<nsite; isite++)
   {
      macro.theta[s[isite]] += 1;
      Latgas_Observables macro_site = site_calc(s,isite);
      for(int i=0; i<5; i++) macro.ebody[i] += macro_site.ebody[i];
#     ifdef LATGAS_FULL_CALC_DEBUG
      site_energy[isite] = macro_site.E;
#     endif
   }
   macro.E = macro.ebody[0];
   for(int i=1; i<5; i++) macro.E += macro.ebody[i]/static_cast<float>(i);
   if(false)
   {
      std::cout << "Latgas_Hamiltonian::full_calc_bysite results:" << std::endl;
      for(int i=0; i<M; i++)
         std::cout << "Num " << species[i] << " = " << macro.theta[i] << std::endl;
      std::cout << "E NBody 0 = " << macro.ebody[0] << " " << macro.ebody[0]/1. << std::endl; 
      std::cout << "E NBody 1 = " << macro.ebody[1] << " " << macro.ebody[1]/1. << std::endl; 
      std::cout << "E NBody 2 = " << macro.ebody[2] << " " << macro.ebody[2]/2. << std::endl; 
      std::cout << "E NBody 3 = " << macro.ebody[3] << " " << macro.ebody[3]/3. << std::endl; 
      std::cout << "E NBody 4 = " << macro.ebody[4] << " " << macro.ebody[4]/4. << std::endl; 
      std::cout << "macro.E  = " << macro.E << "  " << macro.E/nsite << std::endl;
#     ifdef LATGAS_FULL_CALC_DEBUG
      std::ofstream sout("site_energy_bysite.out");
      sout <<"site, sigma_site, E_site, E_site/nsite" << std::endl;
      for(int i=0; i<nsite; i++) sout << i << " " << static_cast<int>(s[i]) << " " << site_energy[i] << " " << site_energy[i]/nsite << std::endl;
#     endif
   }
   return;
}


// Change the value of sigma_i at one grid point, update observables
void Latgas_Hamiltonian::change(Latgas_Hamiltonian::Config& sigma, Latgas_Observables& macro, int igrid, char new_sigma)
{
   char old_sigma = sigma[igrid];
   // Use full_calc() for debugging 
   if(false)
   {
      sigma[igrid] = new_sigma;
      macro = this->site_calc(sigma,igrid);
      return;
   }
   // calculate changes
   Latgas_Observables macro0 = site_calc(sigma,igrid);
   sigma[igrid] = new_sigma;
   Latgas_Observables macro1 = site_calc(sigma,igrid);
   macro.theta[old_sigma]--;
   macro.theta[new_sigma]++;
   macro.E += macro1.E-macro0.E;
   for(int i=0; i<5; i++) macro.ebody[i] += macro1.ebody[i]-macro0.ebody[i];
}


// propose a Monte Carlo move
template<typename MCWalker, typename URNG>
void Latgas_Hamiltonian::spin_flip(MCWalker& walker, URNG& urng)
{
   // Randomly choose a move
   int NSite = walker.sigma.size();
   int isite = static_cast<int>(NSite*urng());
   while( isite<0 || isite>=NSite )
   {
      std::cerr << __FILE__ << ":" << __LINE__ << " markov_propose bad isite=" << isite << "/" << NSite << std::endl;
      isite = static_cast<int>(NSite*urng());
   }
   char sold = walker.sigma[isite];
   char snow = static_cast<char>(M*urng());
   while( snow>=M && snow==sold ) snow = static_cast<char>(M*urng());
   walker.save_initial();
   walker.add_change(isite);
   this->change(walker.sigma,walker.now,isite,snow);
   if( false )
   {
      // Check that the method gets the right total energy
      double Ealgorithm = walker.now.E;
      walker.now = this->full_calc(walker.sigma);
      static int ialert = 0;
      if( std::fabs(Ealgorithm-walker.now.E)>1.e-5 && ialert<10 )
      {
         ialert++;
         std::cout << "Latgas_Hamiltonian::spin_flip gets wrong energy" << std::endl;
         std::cout << "Algorithm energy = " << Ealgorithm << std::endl;
         std::cout << "full_calc energy = " << walker.now.E << std::endl;
      } 
   }
}


// propose a Monte Carlo move
template<typename MCWalker, typename URNG>
void Latgas_Hamiltonian::spin_exch(MCWalker& walker, URNG& urng)
{
   // Randomly choose a move
   int NSite = walker.sigma.size();
   int isite = static_cast<int>(NSite*urng());
   while( isite<0 || isite>=NSite )
   {
      std::cerr << __FILE__ << ":" << __LINE__ << " markov_propose bad isite=" << isite << "/" << NSite << std::endl;
      isite = static_cast<int>(NSite*urng());
   }
   // Random choose a neighbor
   if( exch_disp.size()==0 )
   {  // By default use nearest neighbors
      GridNghbr nndata = FindNearNeighbors(grid);
      exch_disp = nndata.alldisp;
   } 
   int NumN = exch_disp.size();
   int jsite = -1;
   while( jsite<0 || jsite>=NSite )
   {
      int jdisp = static_cast<int>(NSite*urng());
      while( jdisp<0 || jdisp>=NumN )
         jdisp = static_cast<int>(NSite*urng());
      jsite = grid.neighbor(isite,exch_disp[jdisp]); 
   }
   // Do the move
   Latgas_Observables before_exchg(walker.now);   // save macrostate before composing invidual moves
   walker.old.copy(before_exchg);          // This is the old macrostate for the Monte Carlo step
   walker.save_initial();
   walker.add_change(isite);
   walker.add_change(jsite);
   char s_i = walker.sigma[isite];
   char s_j = walker.sigma[jsite];
   this->change(walker.sigma,walker.now,isite,s_j);      // spin_flip updates walker.old as build new state
   this->change(walker.sigma,walker.now,jsite,s_i);
   if( false )
   {
      // Check that the method gets the right total energy
      double Ealgorithm = walker.now.E;
      walker.now = this->full_calc(walker.sigma);
      static int ialert = 0;
      if( std::fabs(Ealgorithm-walker.now.E)>1.e-5 && ialert<10 )
      {
         ialert++;
         std::cout << "Latgas_Hamiltonian::spin_exchange gets wrong energy" << std::endl;
         std::cout << "Algorithm energy = " << Ealgorithm << std::endl;
         std::cout << "full_calc energy = " << walker.now.E << std::endl;
      } 
   }
}

// Make an attempted Monte Carlo step
template<typename MCWalker, typename URNG>
void Latgas_Hamiltonian::mc_step(MCWalker& walker, URNG& urng)
{
   if( SPINFLIP )
      this->spin_flip(walker,urng);
   else
      this->spin_exch(walker,urng);
}




#endif   /* MCMODEL_HPP */

