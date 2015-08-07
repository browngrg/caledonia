#ifndef MC_WALKER_HPP
#define MC_WALKER_HPP


#include<vector>

                                
template<typename Config, typename Observables>
class MC_Walker 
{
public:

   typedef Config ConfigType;
   typedef Observables ObserveType;

   Config sigma;
   Observables now,old;

   int       iwalk_global;      // unique, global identifier
   long long imcs;              // number of Monte Carlo updates, extensive Monte Carlo time
   double    MCSS;              // number of Monte Carlo Step per Site, intensive Monte Caro time

   // Save the state before building proposed move
   void save_initial() 
   { 
      old = now;
      undo.clear(); 
   }

   // Add one change in the chain building the proposed move
   void add_change(int igrid) { undo.push_back( UndoRec(igrid,sigma[igrid]) ); }

   // Undo all of the changes involved in the proposed move
   void restore_initial()
   {
      now = old;
      while( undo.size()>0 )
      {
         sigma[ undo.back().igrid ] = undo.back().sigmai;
         undo.pop_back();
      }
   }

   // Deep copy of Walker
   void copy(const MC_Walker<Config,Observables>& orig)
   {
      sigma = orig.sigma;
      now = orig.now;
      old = orig.old;
      imcs = orig.imcs;
      MCSS = orig.MCSS;
   }

private:

   // History used to roll back (proposed) Mote Carlo moves

   struct UndoRec 
   { 
      int igrid; 
      typename Config::value_type  sigmai; 
      UndoRec(int i, double s) { igrid=i; sigmai=s; }
      UndoRec(const UndoRec& orig) { igrid=orig.igrid; sigmai=orig.sigmai; }
   };

   std::vector<UndoRec> undo;
 
};


#endif
