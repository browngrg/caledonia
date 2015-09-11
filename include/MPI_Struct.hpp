#ifndef MPI_STRUCT_HPP_
#define MPI_STRUCT_HPP_

#ifdef USE_MPI
#include<mpi.h>
#endif

// Combine the basic MPI info into one structure
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

   static MPI_Struct local()
   {
      MPI_Struct w = world();
      MPI_Struct tmp;
#     ifdef USE_MPI
      MPI_Comm_split(w.comm,w.iproc,0,&(tmp.comm));
      MPI_Comm_rank(tmp.comm,&tmp.iproc);
      MPI_Comm_size(tmp.comm,&tmp.nproc);
#     else
      tmp.iproc = 0;
      tmp.nproc = 1;
#     endif
      return tmp;
   }
};

#endif
