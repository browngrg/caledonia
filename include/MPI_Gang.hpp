#ifndef MPI_GANG_HPP
#define MPI_GANG_HPP

#include<iostream>
#include<vector>

#ifdef USE_MPI
#include<mpi.h>
#endif

struct MPI_Struct
{
public:
   int iproc,nproc;
#  ifdef USE_MPI
   MPI_Group group;
   MPI_Comm comm;
#  else
   int group;
   int comm;
#  endif
public:

   MPI_Struct() 
   { 
      iproc=0; 
      nproc=1; 
#     ifdef USE_MPI
      group = MPI_GROUP_NULL;
      comm = MPI_COMM_NULL;
#     else
      group = 0;
      comm = 0;
#     endif
   }

   void copy(const MPI_Struct& s) 
   { 
      iproc=s.iproc;
      nproc=s.nproc;
      group=s.group;
      comm =s.comm;
   }

   void operator=(const MPI_Struct& s) { this->copy(s); }

   bool in() const 
   { 
#  ifdef USE_MPI
      return comm!=MPI_COMM_NULL; 
#  else
      return true;
#  endif
   }


   static MPI_Struct world() 
   {
      MPI_Struct tmp;
#     ifdef USE_MPI
      tmp.comm = MPI_COMM_WORLD;
      MPI_Comm_rank(tmp.comm,&tmp.iproc);
      MPI_Comm_size(tmp.comm,&tmp.nproc);
#     else
      tmp.iproc = 0;
      tmp.nproc = 1;
#     endif
      return tmp;
   }

};


// Break a pool of MPI processes into Gangs that work together
struct MPI_Gang
{
private:

   enum  { verbose_debug = false };          // very verbose tracing

public:

   int igang, ngang;
   MPI_Struct pool, gang, lead;
   bool owner;

public:

   MPI_Gang() 
   { 
      igang = 0; 
      ngang = 1; 
      owner = false;
   }

   MPI_Gang(const MPI_Gang& orig) { this->copy(orig); }

   ~MPI_Gang()
   {
#     ifdef USE_MPI
      if( !owner ) return;
      int final_flag;
      MPI_Finalized(&final_flag);
      if( final_flag ) return;
      if( pool.group!=MPI_GROUP_NULL ) MPI_Group_free(&pool.group);
      if( gang.group!=MPI_GROUP_NULL ) MPI_Group_free(&gang.group);
      if( lead.group!=MPI_GROUP_NULL ) MPI_Group_free(&lead.group);
      if(  pool.comm!=MPI_COMM_NULL  ) MPI_Comm_free(&pool.comm);
      if(  gang.comm!=MPI_COMM_NULL  ) MPI_Comm_free(&gang.comm);
      if(  lead.comm!=MPI_COMM_NULL  ) MPI_Comm_free(&lead.comm);
#     endif
   }

   void operator=(const MPI_Gang& orig) { this->copy(orig); }

   void copy(const MPI_Gang& orig)
   {
      igang = orig.igang;
      ngang = orig.ngang;
      pool  = orig.pool;
      gang  = orig.gang;
      lead  = orig.lead;
      owner = false;
   }

   
   void set_pool(const MPI_Struct& _pool)
   {
      pool.copy(_pool);
#     ifdef USE_MPI
      if(pool.in()) MPI_Comm_dup(_pool.comm,&pool.comm);
#     endif
   }


   void init(int NGang, const MPI_Struct& _pool)
   {
      this->set_pool(_pool);
      this->init(NGang);
   }


   void init(int NGang=1)
   {
#     ifndef USE_MPI
      std::cout << __FILE__ << ":" << __LINE__ << " MPI_Gang assumes -DUSE_MPI" << std::endl;
      return;
#     else
      owner = true;
      ngang = NGang;
      //if(pool.comm==MPI_COMM_NULL) this->set_pool(MPI_Struct::world());
      if(verbose_debug) std::cout << __FILE__ << ":" << __LINE__ << " NGang=" << ngang << std::endl;
      // Get information about the pool of processes
      if( !pool.in() ) return;
      MPI_Comm_rank(pool.comm, &pool.iproc);
      MPI_Comm_size(pool.comm, &pool.nproc);
      MPI_Comm_group(pool.comm,&pool.group);
      // Create a communicator for each Gang of workers
      int NPerGang   = pool.nproc/ngang;
      if( NPerGang<1 ) 
      {
         NPerGang = 1;
         ngang    = pool.nproc;
      }
      if( (pool.nproc % NPerGang) != 0 )
      {
         if( pool.iproc==0 ) 
            std::cout << __FILE__ << ":" << __LINE__ << "Can't evenly divide processes into gangs" << std::endl
                      << "ngang=" << ngang << " npool=" << pool.nproc << std::endl;
      }
      igang = pool.iproc/NPerGang;
      int gang_range[ngang][3];
      for(int ig=0; ig<ngang; ig++)
      {
         gang_range[ig][0] = ig*NPerGang;           // First process in gang
         gang_range[ig][1] = (ig+1)*NPerGang-1;     // Last process in gang
         gang_range[ig][2] = 1;                     // Stride through original group
         if(gang_range[ig][1]>=pool.nproc)
            gang_range[ig][1] = pool.nproc-1;
      }
      if(verbose_debug && pool.iproc==0)
      {
         std::cout << "range =";
         for(int i=0; i<(ngang*3); i++) std::cout << " " << gang_range[igang][i];
         std::cout << std::endl;
      }
      if(pool.in()) 
      {
         MPI_Group_range_incl(pool.group,1,&(gang_range[igang]),&gang.group);
         MPI_Comm_create(pool.comm,gang.group,&gang.comm);
         if( gang.in() )
         {
            MPI_Comm_rank(gang.comm,&gang.iproc);
            MPI_Comm_size(gang.comm,&gang.nproc);
         }
      }
      else
      {
         std::cout << "iproc=" << pool.iproc << " not assigned to a gang" << std::endl;
         return;
      }
      if(verbose_debug) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << pool.iproc << " igang=" << igang << "/" << ngang 
                                  << " iproc_gang=" << gang.iproc << "/" << gang.nproc << std::endl;
      // Create a communicator for lead-processes
      std::vector<int> lead_rank(ngang);
      for(int i=0; i<ngang; i++)
         lead_rank[i] = i*NPerGang;
      if( pool.in() ) 
      {
         MPI_Group_incl(pool.group,ngang,&(lead_rank[0]),&lead.group);
         MPI_Comm_create(pool.comm,lead.group,&lead.comm);
         if( lead.in() )
         {
            MPI_Comm_rank(lead.comm,&lead.iproc);
            MPI_Comm_size(lead.comm,&lead.nproc);
         }
      }
      if( gang.in() && gang.iproc==0 && !lead.in() )   // This is a paranoid check
         std::cout << __FILE__ << ":" << __LINE__ << "(" << pool.iproc << ") gang.iproc=" << gang.iproc << std::endl;
      if( lead.in() )
      {
         int ilead = -1;
         MPI_Comm_rank(lead.comm,&ilead);
         if( ilead!=igang ) 
         {
            std::cout << "ilead!=igang for iproc=" << pool.iproc << " ilead=" << ilead << " igang=" << igang << " iproc_gang=" << gang.iproc << std::endl;
         }
      }
#     endif
   }

};


#endif  // MPI_Gang.hpp
