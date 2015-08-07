#ifndef MPI_GANG_HPP
#define MPI_GANG_HPP

#include<iostream>
#include<vector>

#ifdef USE_MPI
#include<mpi.h>
#endif


// Break a pool of MPI processes into Gangs that work together
struct MPI_Gang
{
private:

   enum  { verbose_debug = false };          // very verbose tracing

public:
   int igang,ngang;
   int iproc_pool, nproc_pool;
   int iproc_gang, nproc_gang;
#  ifdef USE_MPI
   bool owner;
   MPI_Group group_pool, group_gang, group_lead;
   MPI_Comm  comm_pool, comm_gang, comm_lead;
#  endif

public:

   MPI_Gang() 
   { 
      igang      = 0; 
      ngang      = 1; 
      iproc_pool = 0; 
      nproc_pool = 1; 
      iproc_gang = 0;
      nproc_gang = 1;
#     ifdef USE_MPI
      owner = false;
      group_pool = MPI_GROUP_NULL;
      group_gang = MPI_GROUP_NULL;
      group_lead = MPI_GROUP_NULL;
      comm_pool  = MPI_COMM_NULL;
      comm_gang  = MPI_COMM_NULL;
      comm_lead  = MPI_COMM_NULL;
#     endif
   }

   MPI_Gang(const MPI_Gang& orig) { this->copy(orig); }

   ~MPI_Gang()
   {
#     ifdef USE_MPI
      if( !owner ) return;
      int final_flag;
      MPI_Finalized(&final_flag);
      if( final_flag ) return;
      if( group_pool!=MPI_GROUP_NULL ) MPI_Group_free(&group_pool);
      if( group_gang!=MPI_GROUP_NULL ) MPI_Group_free(&group_gang);
      if( group_lead!=MPI_GROUP_NULL ) MPI_Group_free(&group_lead);
      if(  comm_pool!=MPI_COMM_NULL  ) MPI_Comm_free(&comm_pool);
      if(  comm_gang!=MPI_COMM_NULL  ) MPI_Comm_free(&comm_gang);
      if(  comm_lead!=MPI_COMM_NULL  ) MPI_Comm_free(&comm_lead);
#     endif
   }

   void operator=(const MPI_Gang& orig) { this->copy(orig); }

   void copy(const MPI_Gang& orig)
   {
      igang = orig.igang;
      ngang = orig.ngang;
      iproc_pool = orig.iproc_pool;
      nproc_pool = orig.nproc_pool;
      iproc_gang = orig.iproc_gang;
      nproc_gang = orig.nproc_gang;
#     ifdef USE_MPI
      owner = false;
      group_pool = orig.group_pool;
      group_gang = orig.group_gang;
      group_lead = orig.group_lead;
      comm_pool  = orig.comm_pool;
      comm_gang  = orig.comm_gang;
      comm_lead  = orig.comm_lead;
#     endif
   }


   // Methods that set the pool of workers can
   // only be defined if MPI is available because
   // they need MPI types
#  ifdef USE_MPI
   void set_pool(const MPI_Comm& comm)
   {
      MPI_Comm_dup(comm,&comm_pool);
   }

   void init(int NGang, const MPI_Comm& comm)
   {
      this->set_pool(comm);
      this->init(NGang);
   }
#  endif


   void init(int NGang=1)
   {
#     ifndef USE_MPI
      std::cout << __FILE__ << ":" << __LINE__ << " MPI_Gang assumes -DUSE_MPI" << std::endl;
      return;
#     else
      owner = true;
      ngang = NGang;
      if(comm_pool==MPI_COMM_NULL) this->set_pool(MPI_COMM_WORLD);
      if(verbose_debug) std::cout << __FILE__ << ":" << __LINE__ << " NGang=" << ngang << std::endl;
      // Get information about the pool of processes
      MPI_Comm_rank(comm_pool,&iproc_pool);
      MPI_Comm_size(comm_pool,&nproc_pool);
      MPI_Comm_group(comm_pool,&group_pool);
      // Create a communicator for each Gang of workers
      int NPerGang   = nproc_pool/ngang;
      if( NPerGang<1 ) 
      {
         NPerGang = 1;
         ngang    = nproc_pool;
      }
      if( (nproc_pool % NPerGang) != 0 )
      {
         if( iproc_pool==0 ) 
            std::cout << __FILE__ << ":" << __LINE__ << "Can't evenly divide processes into gangs" << std::endl
                      << "ngang=" << ngang << " npool=" << nproc_pool << std::endl;
      }
      igang = iproc_pool/NPerGang;
      int gang_range[ngang][3];
      for(int ig=0; ig<ngang; ig++)
      {
         gang_range[ig][0] = ig*NPerGang;           // First process in gang
         gang_range[ig][1] = (ig+1)*NPerGang-1;     // Last process in gang
         gang_range[ig][2] = 1;                     // Stride through original group
         if(gang_range[ig][1]>=nproc_pool)
            gang_range[ig][1] = nproc_pool-1;
      }
      if(verbose_debug && iproc_pool==0)
      {
         std::cout << "range =";
         for(int i=0; i<(ngang*3); i++) std::cout << " " << gang_range[igang][i];
         std::cout << std::endl;
      }
      MPI_Group_range_incl(group_pool,1,&(gang_range[igang]),&group_gang);
      MPI_Comm_create(comm_pool,group_gang,&comm_gang);
      if( comm_gang!=MPI_COMM_NULL )
      {
         MPI_Comm_rank(comm_gang,&iproc_gang);
         MPI_Comm_size(comm_gang,&nproc_gang);
      }
      else
      {
         std::cout << "iproc=" << iproc_pool << " not assigned to a gang" << std::endl;
         return;
      }
      if(verbose_debug) std::cout << __FILE__ << ":" << __LINE__ << " iproc=" << iproc_pool << " igang=" << igang << "/" << ngang 
                                  << " iproc_gang=" << iproc_gang << "/" << nproc_gang << std::endl;
      // Create a communicator for lead-processes
      std::vector<int> lead_rank(ngang);
      for(int i=0; i<ngang; i++)
         lead_rank[i] = i*NPerGang;
      MPI_Group_incl(group_pool,ngang,&(lead_rank[0]),&group_lead);
      MPI_Comm_create(comm_pool,group_lead,&comm_lead);
      if( iproc_gang==0 )
      {
         int ilead = -1;
         MPI_Comm_rank(comm_lead,&ilead);
         if( ilead!=igang ) 
         {
            std::cout << "ilead!=igang for iproc=" << iproc_pool << " ilead=" << ilead << " igang=" << igang << " iproc_gang=" << iproc_gang << std::endl;
         }
      }
#     endif
   }

};


#endif  // MPI_Gang.hpp
