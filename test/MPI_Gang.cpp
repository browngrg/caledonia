

#include"MPI_Gang.hpp"

#include<iostream>


void test1(const MPI_Comm& pool)
{
   MPI_Gang mp;
   mp.init(5,pool);    // split processes into 5 gangs

   if(mp.iproc_pool==0) std::cout << "init successful" << std::endl;

   int ival = 0;
   for(int i=0; i<5; i++) ival += i*mp.nproc_gang;   // 10*mp.nproc_gang == 2*mp.nproc_pool

   // Sum over just leaders
   int isum = 0;
   if(mp.comm_lead!=MPI_COMM_NULL) MPI_Reduce( &(mp.iproc_pool), &(isum), 1, MPI_INT, MPI_SUM, 0, mp.comm_lead );
   if( mp.iproc_pool==0 )
      std::cout << "sum of leader ids = " << isum << std::endl;

   // Sum over all processes
   MPI_Reduce( &(mp.iproc_pool), &(isum), 1, MPI_INT, MPI_SUM, 0, mp.comm_pool );
   if( mp.iproc_pool==0  )
      std::cout << "sum of all ids = " << isum << std::endl;

   // Every processor get result from sum over leaders
   if(mp.comm_lead!=MPI_COMM_NULL) MPI_Allreduce( &(mp.iproc_pool), &(isum), 1, MPI_INT, MPI_SUM, mp.comm_lead );
   MPI_Bcast(&(isum),1,MPI_INT,0,mp.comm_gang);
   if( isum!=ival )
      std::cout << " Wrong result from Bcast iproc=" << mp.iproc_pool << " isum=" << isum << " expected=" << ival << std::endl;
}


int main(int argc, char* argv[])
{

   MPI_Init(&argc,&argv);
   test1(MPI_COMM_WORLD);
   MPI_Finalize();

}
