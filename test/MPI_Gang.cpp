

#include"MPI_Gang.hpp"

#include<iostream>



void test1(const MPI_Struct& pool)
{
   MPI_Gang mp;
   mp.init(5,pool);    // split processes into 5 gangs

   if(mp.pool.iproc==0) std::cout << "init successful" << std::endl;

   int ival = 0;
   for(int i=0; i<5; i++) ival += i*mp.gang.nproc;   // 10*mp.nproc_gang == 2*mp.nproc_pool

   // Sum over just leaders
   int isum = 0;
   if( mp.lead.in() ) MPI_Reduce( &(mp.pool.iproc), &(isum), 1, MPI_INT, MPI_SUM, 0, mp.lead.comm );
   if( mp.pool.in() && mp.pool.iproc==0 )
      std::cout << "sum of leader ids = " << isum << std::endl;

   // Sum over all processes
   if( mp.pool.in() ) MPI_Reduce( &(mp.pool.iproc), &(isum), 1, MPI_INT, MPI_SUM, 0, mp.pool.comm );
   if( mp.pool.in() && mp.pool.iproc==0  )
      std::cout << "sum of all ids = " << isum << std::endl;

   // Every processor get result from sum over leaders
   if( mp.lead.in() ) MPI_Allreduce( &(mp.pool.iproc), &(isum), 1, MPI_INT, MPI_SUM, mp.lead.comm );
   if( mp.gang.in() ) MPI_Bcast(&(isum),1,MPI_INT,0,mp.gang.comm);
   if( isum!=ival )
      std::cout << " Wrong result from Bcast iproc=" << mp.pool.iproc << " isum=" << isum << " expected=" << ival << std::endl;
}


void test2(const MPI_Struct& pool)
{
   MPI_Gang bigger;
   bigger.init(10,pool);
   if( bigger.lead.in() && bigger.lead.nproc<5 && bigger.lead.iproc==0 )
   {
      std::cout << __FILE__ << ":" << __LINE__ << "(" << bigger.pool.iproc << ") " << bigger.lead.nproc << std::endl;
      return;
   }
   MPI_Barrier(pool.comm);
   if( bigger.lead.in() ) test1(bigger.lead);
}


int main(int argc, char* argv[])
{

   MPI_Init(&argc,&argv);
   test1(MPI_Struct::world());
   test2(MPI_Struct::world());
   MPI_Finalize();

}
