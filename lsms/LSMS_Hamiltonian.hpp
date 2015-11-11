#ifndef CALEDONIA_LSMS_HAMILTONIAN_HPP
#define CALEDONIA_LSMS_HAMILTONIAN_HPP


#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <iostream>
#include <fstream>
#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL/VORPOL.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "Potential/PotentialShifter.hpp"
#include "calculateChemPot.hpp"
#include "lsmsClass.hpp"
#include "EvecGenerator.h"
#include "WangLandau.h"
#include "ExhaustiveIsing.h"
#include "WangLandau2d.h"

#include "Heisenberg_Hamiltonian.hpp"


// Occupation simulations in magnetic alloys demand evec
struct LSMS_Config
{
public:
   double energy;
   double band_energy;
   std::vector<double> evec;
   std::vector<int> occ;
   int init_steps;                            // This eventually becomes MC_Walker.time
   std::vector<double> vSpinShifts;
   std::vector<double> evecsAndSpinShifts;
};




class LSMS_Hamiltonian
{
public:

  int size, rank, world_rank, my_group, comm_size;
  int num_lsms;                 // number of parallel LSMS instances
  int size_lsms;                // number of atoms in a lsms instance
  int num_steps;                // number of energy calculations
  int initial_steps;            // number of steps before sampling starts
  int stepCount = 0;            // count the Monte Carlo steps executed
  double max_time;              // maximum walltime for this run in seconds
  bool restrict_time = false;   // was the maximum time specified?
  bool restrict_steps = false;  // or the max. numer of steps?
  int align;                    // alignment of lsms_instances

  double magnetization;
  double energy_accumulator;    // accumulates the enegy to calculate the mean
  int energies_accumulated;

  double walltime_0, walltime;

  double restartWriteFrequency = 30.0 * 60.0;
  double nextWriteTime = restartWriteFrequency;

  MPI_Comm local_comm;
  MPI_Status status;

  std::vector<int> lsms_rank0;

  char prefix[40];
  char i_lsms_name[64];
  char gWL_in_name[64], gWL_out_name[64];
  char mode_name[64];
  char energy_calculation_name[64];
  char stupid[37];
  std::vector<long> walkerSteps;

  char step_out_name[64];
  char wl_step_out_name[128];
  char *wl_stepf = NULL;
  bool step_out_flag = false;
  std::ofstream step_out_file;
  typedef enum {Constant, Random, WangLandau_1d, ExhaustiveIsing, WangLandau_2d} EvecGenerationMode;
  typedef enum {MagneticMoment, MagneticMomentZ, MagneticMomentX, MagneticMomentY} SecondDimension;
  typedef enum {OneStepEnergy, MultiStepEnergy, ScfEnergy} EnergyCalculationMode;

  bool isSpinSim = true;
  bool isOccupancySim = false;
  char config_space_name[64];
  AlloyMixingDesc alloyDesc;

  double ev0[3];

  bool return_moments_flag; // true -> return all magnetic moments from lsms run at each step.
  bool generator_needs_moment;

  EnergyCalculationMode energyCalculationMode;

  bool evec_random_flag;

  EvecGenerationMode evec_generation_mode;
  SecondDimension second_dimension;
  MoveChoice_t MoveChoice, prevMoveChoice;

  int energyIndex; // index for the return value to use for the MC step (0: total energy, 1: band energy)

   std::vector<double> recv_buffer;
   std::vector<double> r_values;
   std::vector<double> i_values;

   void calc_energy(LSMS& calc, LSMS_Config& sigma, MoveChoice_t local_choice);
   void set_config(LSMS& calc, LSMS_Config& sigma, MoveChoice_t local_choice);
   void send_num_sites(LSMS& calc);
   void send_shifter_info(LSMS& calc);

public:

  LSMS_Hamiltonian()
  {
     energy_accumulator = 0.0;
     energies_accumulated = 0;
     evec_generation_mode = Constant;
     evec_random_flag = false;
     second_dimension = MagneticMoment;
     MoveChoice, prevMoveChoice;
     return_moments_flag = true;
     generator_needs_moment = false;
     energyCalculationMode = OneStepEnergy;
     energyIndex = 1;
     ev0[0] = ev0[1] = 0.0; ev0[2] = 1.0;
     // size has to be align + size_lsms*num_lsms
     align = 1;
     num_lsms = 1;
     size_lsms = -1;
     my_group = -1;
     num_steps = -1;
     initial_steps = 0;
     walltime = -1;
     sprintf(i_lsms_name, "i_lsms");
     gWL_in_name[0] = gWL_out_name[0] = 0;
     mode_name[0] = 0;
     energy_calculation_name[0] = 0;
     config_space_name[0] = 0;
  }

public:

  // check command line arguments
  void parse_command_line(int argc, char* argv[]);
  template<typename OPTIONS> void add_options(OPTIONS& options);

  void init();
  void old_dosim();

private:

   enum { R_VALUE_OFFSET = 3 };

private:

   void set_modes();
   void write_header();
   void do_slave();
   void do_master(); 

   int  get_op(LSMS_Config& sigma, std::vector<double>& recv_buffer);
   void send_alloy_description(LSMS& lsms_calc);
   void recieve_alloy_description(EvecGenerator* generator);

   void send_evec(LSMS_Config& sigma, int itarget, bool vSpinShiftFlag);
   void send_occ (LSMS_Config& sigma, int itarget);

};


void LSMS_Hamiltonian::parse_command_line(int argc, char* argv[])
{
  for (int i=0; i<argc; i++)
  {
    if (!strcmp("-num_lsms", argv[i])) num_lsms = atoi(argv[++i]);
    if (!strcmp("-size_lsms", argv[i])) size_lsms = atoi(argv[++i]);
    if (!strcmp("-align", argv[i])) align = atoi(argv[++i]);
    if (!strcmp("-num_steps", argv[i])) {
      num_steps = atoi(argv[++i]);
      restrict_steps = true;
    }
    if (!strcmp("-initial_steps", argv[i])) initial_steps = atoi(argv[++i]); 
    if (!strcmp("-walltime", argv[i])) {
      max_time = 60.0*atof(argv[++i]);
      restrict_time = true;
    }
    if (!strcmp("-i", argv[i])) strncpy(i_lsms_name, argv[++i], 64);
    if (!strcmp("-random_dir", argv[i])) evec_generation_mode = Random;
    if (!strcmp("-step_out", argv[i])) {
      strncpy(step_out_name, argv[++i],64);
      step_out_flag = true;
      return_moments_flag = true;
    }
    if (!strcmp("-wl_out", argv[i])) strncpy(gWL_out_name, argv[++i], 64);
    if (!strcmp("-wl_in", argv[i])) strncpy(gWL_in_name, argv[++i], 64);
    if (!strcmp("-mode", argv[i])) strncpy(mode_name, argv[++i], 64);
    if (!strcmp("-energy_calculation", argv[i])) strncpy(energy_calculation_name, argv[++i], 64);
    if (!strcmp("-config_space", argv[i])) strncpy(config_space_name, argv[++i], 64);
  }
  if (!(restrict_steps || restrict_time)) restrict_steps = true;
}


template<typename OPTIONS> 
void LSMS_Hamiltonian::add_options(OPTIONS& options)
{
  options.add_option("num_lsms",      "Number of lsms calculators",     ' ', &(num_lsms) );
  options.add_option("size_lsms",     "Number of sites in calculation", ' ', &(size_lsms) );
  options.add_option("align",         "",                               ' ', &(align) );
  options.add_option("num_steps",     "Number Wang-Landau steps",       ' ', &(num_steps) ); 
  options.add_option("initial_steps", "",                               ' ', &(initial_steps) ); 
  options.add_option("walltime",      "Number minutes to run",          ' ', &(walltime) ); 
  options.add_option("ilsms",         "lsms input file",                ' ', i_lsms_name );
  options.add_option("random_dir",    "random direction for evecs",     ' ', &(evec_random_flag) );
  options.add_option("step_out",      "specify step output file",       ' ', step_out_name );
  options.add_option("wl_out",        "specify wanglandau output file", ' ', gWL_out_name );
  options.add_option("wl_int",        "specify wanglandau input file",  ' ', gWL_in_name );
  options.add_option("mode",          "specify mode by name",           ' ', mode_name );
  options.add_option("energy_calc",   "specify calculation by name",    ' ', energy_calculation_name );
  options.add_option("config_space",  "specify configuration by name",  ' ', config_space_name );
}


void LSMS_Hamiltonian::set_modes()
{
  // determine whether we perform simulation over occupancy or spin variables
  // in principle both can be done, but this has not been tested
  if( config_space_name[0] != 0 ) { 
    if (config_space_name[0] == 'o') {
      isOccupancySim = true;
      isSpinSim = false;
    }
  }
  if (mode_name[0] != 0)
  {
    if (!strcmp("constant", mode_name)) evec_generation_mode = Constant;
    if (!strcmp("random", mode_name)) evec_generation_mode = Random;
    if (!strcmp("1d", mode_name)) evec_generation_mode = WangLandau_1d;
    if (!strcmp("ising", mode_name)) evec_generation_mode = ExhaustiveIsing;
    if (!strcmp("2d", mode_name)) evec_generation_mode = WangLandau_2d;
    if (!strcmp("2d-m", mode_name)) 
    {
      evec_generation_mode = WangLandau_2d;
      second_dimension = MagneticMoment;
    }
    if (!strcmp("2d-x", mode_name)) 
    {
      evec_generation_mode = WangLandau_2d;
      second_dimension = MagneticMomentX;
    }
    if (!strcmp("2d-y", mode_name)) 
    {
      evec_generation_mode = WangLandau_2d;
      second_dimension = MagneticMomentY;
    }
    if (!strcmp("2d-z", mode_name)) 
    {
      evec_generation_mode = WangLandau_2d;
      second_dimension=MagneticMomentZ;
    }
  }
  if (energy_calculation_name[0] != 0)
  {
    if(energy_calculation_name[0] == 'o') 
    {
      energyCalculationMode = OneStepEnergy;
      energyIndex = 1; 
    }
    if(energy_calculation_name[0] == 'm')
    {
      energyCalculationMode = MultiStepEnergy;
      energyIndex = 1;
    }
    if(energy_calculation_name[0] == 's')
    {
      energyCalculationMode = ScfEnergy;
      energyIndex = 0;
    }
  }
  // make sure 'return_moments_flag' is set correctly
  switch (evec_generation_mode)
  {
    case Constant : break;
    case Random : break;
    case WangLandau_1d :
      return_moments_flag = true;
      generator_needs_moment = true;
      break;
    case ExhaustiveIsing : break;
    case WangLandau_2d :
      return_moments_flag = true;
      generator_needs_moment = true;
      break;
    default:
      std::cout << " ERROR: UNKNOWN EVEC GENERATION MODE\n";
      exit(1);
  }
}


void LSMS_Hamiltonian::write_header()
{
  if(rank != 0) return;
    std::cout << "LSMS_3" << std::endl;
//    std::cout << " SVN revision " << SVN_REV << std::endl << std::endl;
//#ifdef USE_PAPI
//    std::cout << " Using Papi counters" << std::endl << std::endl; 
//#endif
    std::cout << " Size of LSMS instances = " << size_lsms << " atoms\n";
    std::cout << " Number of LSMS instances = " << num_lsms << std::endl;
    std::cout << " LSMS Energy calculated using ";
    switch (energyCalculationMode)
    {
      case OneStepEnergy:
        std::cout << "oneStepEnergy [frozen potential band energy]" << std::endl;
        break;
      case MultiStepEnergy:
        std::cout << "multiStepEnergy [frozen potential band energy with converged Fermi energy]" << std::endl;
        break;
      case ScfEnergy:
        std::cout << "scfEnergy [self-consistent total energy]" << std::endl;
        break;
      default:
        std::cout << "UNKNOWN ENERGY CALCULATION METHOD" << std::endl;
        exit(1);
    }
    if( isOccupancySim ) 
      std::cout << " Exploring occupancy configurational space\n";
    if( isSpinSim ) 
      std::cout << " Exploring spin configurational space\n";
    if( isOccupancySim && isSpinSim ) {
      std::cout << " Error: Exploring both occupancy and spin configurational\n";
      std::cout << "  has not been tested. Please re-compile, test, then run\n";
      std::cout << "  desired system\n";
      exit(1);
    }
    if (restrict_steps)
      std::cout << " Number of gWL steps = " << num_steps << std::endl;
    if (restrict_time)
      std::cout << " Maximum walltime = " << max_time << "s\n";
    std::cout << " Processor alignment (process allocation quantization) = " << align << std::endl;
    switch (evec_generation_mode)
    {
      case Constant :
        std::cout << " Constant moments direction along "
                  << ev0[0] << " " << ev0[1] << " " << ev0[2] << std::endl;
        break;
      case Random :
        std::cout << " Random distribution of moments (no Wang-Landau)" << std::endl;
        break;
      case WangLandau_1d :
        std::cout << " Wang-Landau for one continuous variable (energy)" << std::endl;
        break;
      case ExhaustiveIsing :
        std::cout << " Exhaustive Ising sampling" << std::endl;
        break;
      case WangLandau_2d :
        std::cout << " Wang-Landau for two continuous variable (energy, ";
        switch (second_dimension)
        {
          case MagneticMoment  : std::cout << "magnitude of magnetization)"; break;
          case MagneticMomentX : std::cout<<"x component of magnetization)"; break;
          case MagneticMomentY : std::cout<<"y component of magnetization)"; break;
          case MagneticMomentZ : std::cout<<"z component of magnetization)"; break;
        }
        std::cout << std::endl;
        break;
      default:
        std::cout << " ERROR: UNKNOWN EVEC GENERATION MODE\n";
        exit(1);
    }
    if (step_out_flag)
      std::cout << " Step output written to: " << step_out_name << std::endl;
    std::cout << std::endl;

    if (step_out_flag && (evec_generation_mode == WangLandau_1d))
    {
      // step_out_flag = false;
      snprintf(wl_step_out_name, 127, "wl1d_%s", step_out_name);
      wl_stepf = wl_step_out_name;
    }

    if (step_out_flag)
    {
      step_out_file.open(step_out_name);
      step_out_file << "#";
      //for(int i=0; i<argc; i++) step_out_file << " " << argv[i];
      step_out_file << std::endl << size_lsms << std::endl;
    }
}


void LSMS_Hamiltonian::init()
{
    // build the communicators
    int color = MPI_UNDEFINED;
    int s = align;
    comm_size = (size-align) / num_lsms;
    for (int i=0; i<num_lsms; i++)
    {
      if ((world_rank >= s) && (world_rank < s+comm_size))
      {
        my_group = i;
        color = i;
      }
      lsms_rank0[i] = s;
      s += comm_size;
    }
    if (world_rank == 0) color = num_lsms;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &local_comm);
    snprintf(prefix, 38, "Group %4d: ", my_group);
}


void LSMS_Hamiltonian::send_alloy_description(LSMS& lsms_calc)
{
   // send alloy description back to master node
   if( rank == 0 ) {
     // browngrg 10/14/2015 Not sure this is debugged. site_alloyclass gets malloced in lsms_calc, and freed here?
     int nclasses;
     int *site_alloyclass;
     lsms_calc.getAlloyInfo(alloyDesc, &(site_alloyclass));
     // ^-- allocates and fills alloy class definition for each site
  
     nclasses = alloyDesc.size();
     MPI_Send(&nclasses, 1, MPI_INTEGER, 0, 1012, MPI_COMM_WORLD);
  
     std::vector<int> ncomps(nclasses);
       for(int i = 0; i < nclasses; i++)
       ncomps[i] = alloyDesc[i].size();
  
     MPI_Send(&(ncomps[0]), nclasses, MPI_INTEGER, 0, 1012, MPI_COMM_WORLD);
     for(int i = 0; i < nclasses; i++) 
       MPI_Send(alloyDesc[i].data(), ncomps[i]*sizeof(AtomType), MPI_BYTE, 0, 1012, MPI_COMM_WORLD);
  
     MPI_Send(&(site_alloyclass[0]), size_lsms, MPI_INTEGER, 0, 1012, MPI_COMM_WORLD);
     free(site_alloyclass);
     // dprintf("Walker %d: Sent master alloy description.\n",my_group);
   }
}


int LSMS_Hamiltonian::get_op(LSMS_Config& sigma, std::vector<double>& recv_buffer)
{
   int op;
   if (rank == 0)
   {
     // recieve command in op code
     // data represents either site spins or site occupancies 
     MPI_Recv(&(recv_buffer[0]), 4*size_lsms, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
     // use op code to distinguish whether occupancy or spin variables recieved
     op = status.MPI_TAG;
     // dprintf("Walker %d: Recieved op code %d.\n", my_group, op);
     if( op == 51 || op == 61 ) {
       for(int i = 0; i < 4*size_lsms; i++)
         sigma.evec[i] = recv_buffer[i];
     }
     else if( op == 52 || op == 62 ) {
       for(int i = 0; i < size_lsms; i++)
         sigma.occ[i] = int(0.5 + recv_buffer[i]);
     }
   }
   MPI_Bcast(&op, 1, MPI_INT, 0, local_comm);
   return op;
}


void LSMS_Hamiltonian::calc_energy(LSMS& lsms_calc, LSMS_Config& sigma, MoveChoice_t local_MoveChoice)
{
    // if move type is Spin, set evec
    //  otherwise set occupancies
    if( local_MoveChoice == SpinMove ) 
    {
      if(lsms_calc.vSpinShiftFlag())
        lsms_calc.setEvecAndSpinPotentialShift4( &(sigma.evec[0]) );
      else
        lsms_calc.setEvec( &(sigma.evec[0]) );
    } 
    else 
    {
      lsms_calc.setOccupancies( &(sigma.occ[0]) );
    }
    // For SCF calculations the Total energy is relevant for MC
    // For Frozen potential calculations its the band energy
    if (energyCalculationMode == OneStepEnergy)
      sigma.energy = lsms_calc.oneStepEnergy(&sigma.band_energy);
    else if (energyCalculationMode == MultiStepEnergy)
      sigma.band_energy = sigma.energy = lsms_calc.multiStepEnergy();
    else if (energyCalculationMode == ScfEnergy)
      sigma.energy = lsms_calc.scfEnergy(&sigma.band_energy);
    else
    {
      std::cout << "ERROR: Unknown energy calculation mode for lsms_calc in wl-lsms main!\n";
      MPI_Abort(MPI_COMM_WORLD, 5);
    }
    r_values[0] = sigma.energy;
    r_values[1] = sigma.band_energy;
    // store the type of move performed as an extra variable in return values
    if( local_MoveChoice == SpinMove )
      r_values[2] = 0x00;
    else if( local_MoveChoice == OccupancyMove )
      r_values[2] = 0x01;

    if (return_moments_flag)
      lsms_calc.getMag( &(r_values[R_VALUE_OFFSET]) );
    if (rank == 0)
    {
      if (return_moments_flag)
        MPI_Send( &(r_values[0]), R_VALUE_OFFSET + 3*size_lsms, MPI_DOUBLE, 0, 1005, MPI_COMM_WORLD);
      else 
        MPI_Send( &(r_values[0]), R_VALUE_OFFSET, MPI_DOUBLE, 0, 1005, MPI_COMM_WORLD);
    }
}


void LSMS_Hamiltonian::set_config(LSMS& lsms_calc, LSMS_Config& sigma, MoveChoice_t local_MoveChoice)
{
    // set configration without calculation
    if( MoveChoice == SpinMove ) {
      if(lsms_calc.vSpinShiftFlag())
        lsms_calc.setEvecAndSpinPotentialShift4( &(sigma.evec[0]) );
      else
        lsms_calc.setEvec( &(sigma.evec[0]) );
    } 
    else {
      if(rank==0) {
        /* dprintf("Walker %d: Setting occupancies: ",my_group);
        for(int i = 0; i < size_lsms; i++)
          dprintf("%d ",occ[i]);
        dprintf("\n"); */
      }

      lsms_calc.setOccupancies( &(sigma.occ[0]) );
    }
}


void LSMS_Hamiltonian::send_num_sites(LSMS& lsms_calc)
{
    // This is never received by the master
    i_values[0] = lsms_calc.numSpins();
    MPI_Send( &(i_values[0]), 10, MPI_INT, 0, 1010, MPI_COMM_WORLD);
}


void LSMS_Hamiltonian::send_shifter_info(LSMS& lsms_calc)
{
    r_values[0] = -1.0; if(lsms_calc.vSpinShiftFlag()==true) r_values[0]=1.0;
    r_values[1]=lsms_calc.potentialMinShift();
    r_values[2]=lsms_calc.potentialMaxShift();
}


/* recognized opcodes:
  5: calculate energy
    51: calculate energy for spin change
    52: calculate energy for occupancy change

  6: set configuration only
    61: set spin variables
    62: set occupancy variables

  recognized energy calculation modes:
  OneStepEnergy : calclulate frozen potential band energy in one step (don't converge Ef)
  use only if the Fermi energy will not change due to MC steps!
  The only method available in LSMS_1.9
  MultiStepEnergy : calculate frozen potential band energy after converging Fermi energy
  This should be the new default method. If the Fermi energy doesn't change
  multiStepEnergy only performs one step and should be equivalent to oneStepEnergy
  The tolerance for Ef convergence can be set with LSMS::setEfTol(Real).
  The default tolerance is set in the LSMS::LSMS constructor (currently 1.0e-6).
  The maximum number of steps is read from the LSMS input file 'nscf' parameter.
  ScfEnergy : this will calculate the selfconsistent total energy.
  The maximum number of steps is read from the LSMS input file 'nscf' parameter.

  10: get number of sites
  11: get potential shifter information
  12: get alloy description
*/
void LSMS_Hamiltonian::do_slave()
{
   LSMS_Config sigma;

   // recieve either 'evec' or 'occupancy' variable set
   recv_buffer.resize(4*size_lsms);
   sigma.occ.resize(size_lsms);
   sigma.evec.resize(4*size_lsms);
   r_values.resize( R_VALUE_OFFSET + 3*(size_lsms+1) );
   i_values.resize(10);

   MPI_Comm_rank(local_comm, &rank);
   snprintf(prefix, 38, "%d_", my_group);
   // to use the ramdisk on jaguarpf:
   // snprintf(prefix, 38, "/tmp/ompi/%d_", my_group);
   // This is from lsmsClass.hpp
   LSMS lsms_calc(local_comm, i_lsms_name, prefix, my_group);
   snprintf(prefix, 38, "Group %4d: ", my_group);

   // wait for commands from master
   bool finished = false;
   while (!finished)
   {
      int op = get_op(sigma,recv_buffer);
      int calc = op/10;
      int imove = op - 10*calc;
      switch(imove)
      {
         default:
         case 1: MoveChoice = SpinMove; break;
         case 2: MoveChoice = OccupancyMove; break;
      }
      switch(calc)
      {
      case 5: 
         calc_energy(lsms_calc,sigma,MoveChoice);
         break;
      case 6: 
         set_config(lsms_calc,sigma,MoveChoice);
         break;
      case 10:
         send_num_sites(lsms_calc);
         break;
      case 11:
         send_shifter_info(lsms_calc);
         break;
      case 12:
         send_alloy_description(lsms_calc);
         break;
      default:
         // printf("world rank %d: recieved exit\n",world_rank);
         if(rank==0) printf("Walker %d: Exiting simulation.\n", my_group);
         if(rank==0) printf("Wang-Landau walker %d performed %ld energy contour integrations.\n",my_group,lsms_calc.energyLoopCount);
         // send this to the root for statistics
         long int sb[2];
         sb[0]=my_group;
         sb[1]=lsms_calc.energyLoopCount;
         MPI_Gather(sb,2,MPI_LONG,NULL,2,MPI_LONG,0,MPI_COMM_WORLD);
         finished = true;
      }
   }
}
  
void LSMS_Hamiltonian::recieve_alloy_description(EvecGenerator *generator)
{
    int nclasses;
    MPI_Recv(&nclasses, 1, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    std::vector<int> ncomps(nclasses);
    std::vector<int> site_alloyclass(size_lsms);
    MPI_Recv(&(ncomps[0]), nclasses, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    alloyDesc.resize(nclasses);
    for(int i = 0; i < nclasses; i++) {
      alloyDesc[i].resize(ncomps[i]);
      MPI_Recv(alloyDesc[i].data(), ncomps[i]*sizeof(AtomType), MPI_BYTE, 
          lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    }
    MPI_Recv(&(site_alloyclass[0]), size_lsms, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    generator->setAlloyClasses(alloyDesc, &(site_alloyclass[0]));
    printf("Master: Recieved alloy description.\n");
}


void LSMS_Hamiltonian::send_evec(LSMS_Config& sigma, int itarget, bool vSpinShiftFlag)
{
   if(vSpinShiftFlag)
   {
     for (int j=0; j<size_lsms; j++)
     {
       sigma.evecsAndSpinShifts[4*j+0] = sigma.evec[3*j+0];
       sigma.evecsAndSpinShifts[4*j+1] = sigma.evec[3*j+1];
       sigma.evecsAndSpinShifts[4*j+2] = sigma.evec[3*j+2];
       sigma.evecsAndSpinShifts[4*j+3] = sigma.vSpinShifts[j];
     }
     MPI_Send(&(sigma.evecsAndSpinShifts[0]), 4*size_lsms, MPI_DOUBLE, itarget /* lsms_rank0[i]*/, 51, MPI_COMM_WORLD);
   }
   else
     MPI_Send(&(sigma.evec[0]), 3*size_lsms, MPI_DOUBLE, itarget /*lsms_rank0[i]*/, 51, MPI_COMM_WORLD);
}


void LSMS_Hamiltonian::send_occ(LSMS_Config& sigma, int itarget)
{
// Does the buffer have to persist after the call?
// std::vector<double> send_buffer(size_lsms);
// for(int j = 0; j < size_lsms; j++)
//   send_buffer[j] = double(sigma.occ[j]);
// MPI_Send( &(send_buffer[0]), size_lsms, MPI_DOUBLE, itarget, 52, MPI_COMM_WORLD);
   MPI_Send( &(sigma.occ[0]), size_lsms, MPI_DOUBLE, itarget, 52, MPI_COMM_WORLD);
}


void LSMS_Hamiltonian::do_master()
    {
      std::vector<LSMS_Config> sigma;

      int running;
      std::vector<double> send_buffer;  // moved to send_occ?
//    std::vector<double> r_values;
//    std::vector<int> i_values(10);
      i_values.resize(10);
//    int *init_steps;
      int total_init_steps;
      bool accepted;

        
      char *wl_inf = NULL;
      char *wl_outf = NULL;
      if (gWL_in_name) wl_inf = gWL_in_name;
      if (gWL_out_name) wl_outf = gWL_out_name;
        
      PotentialShifter potentialShifter;
      EvecGenerator *generator;

      sigma.resize(num_lsms);
      for(int i=0; i<num_lsms; i++)
      {
        sigma[i].evec.resize( 3*size_lsms  );
        sigma[i].occ.resize(size_lsms);
        sigma[i].init_steps = initial_steps;
      }
      if(potentialShifter.vSpinShiftFlag)
      {
        for(int i=0; i<num_lsms; i++)
        {
          sigma[i].vSpinShifts.resize(size_lsms);
          for (int j = 0; j < size_lsms; j++) 
            sigma[i].vSpinShifts[j] = 0.0;
          sigma[i].evecsAndSpinShifts.resize(4*size_lsms);
        }
      }
      send_buffer.resize( 4*size_lsms );

      total_init_steps = num_lsms * initial_steps;
        
      r_values.resize( R_VALUE_OFFSET + 3*(size_lsms+1) );

      //SMURF SMURF SMURF SMURF SMURF SMURF
      //browngrg 10/14/2015 Had to modify generators to deal with moving evecs,occs to LSMS_Config

      // Initialize the correct evec generator
      switch (evec_generation_mode)
      {
        case Random : generator = new RandomEvecGenerator(size_lsms);
          break;
        case Constant: generator = new ConstantEvecGenerator(size_lsms, ev0, num_lsms);
          break;
      case WangLandau_1d : 
          //generator = new WL1dEvecGenerator<std::mt19937>(size_lsms, num_lsms, evecs, potentialShifter, wl_inf, wl_outf, wl_stepf, occs);
          generator = new WL1dEvecGenerator<std::mt19937>(size_lsms, num_lsms, NULL, potentialShifter, wl_inf, wl_outf, wl_stepf, NULL);
          break;
        case ExhaustiveIsing : 
          //generator = new ExhaustiveIsing1dEvecGenerator(size_lsms, num_lsms, evecs, wl_inf, wl_outf);
          generator = new ExhaustiveIsing1dEvecGenerator(size_lsms, num_lsms, NULL, wl_inf, wl_outf);
          break;
        case WangLandau_2d : 
          //generator = new WL2dEvecGenerator<std::mt19937>(size_lsms, num_lsms, evecs, wl_inf, wl_outf, wl_stepf);
          generator = new WL2dEvecGenerator<std::mt19937>(size_lsms, num_lsms, NULL, wl_inf, wl_outf, wl_stepf);
          break;
        default : std::cerr << "The code should never arrive here: UNKNOWN EVEC GENERATION MODE\n";
          exit(1);
      }

      //SMURF SMURF SMURF SMURF SMURF SMURF

      // get alloy description from one of the LSMS instances
      // 'nclasses' is number of alloy mixing classes
      //  components within the same mixing class may substitute for each other
      //  'ncomps[ac]' is number of components for mixing class 'ac'
      //  'site_alloyclass[i]' is the mixing class for site position 'i'
      //  'alloyDesc[ac][i]' is atom description of ith component of mixing class 'ac'

      // signal op code 12 = send master (me) description of alloy
      if( isOccupancySim ) {
        MPI_Send( &(send_buffer[0]), 1, MPI_DOUBLE, lsms_rank0[0], 12, MPI_COMM_WORLD);
        printf("Master: Requesting alloy description.\n");
        recieve_alloy_description(generator);
      }

      // Generate the initial spin configuration 
      if (potentialShifter.vSpinShiftFlag)
        for(int i=0; i<num_lsms; i++)
          generator -> initializeEvecAndPotentialShift(i, &(sigma[i].evec[0]), &(sigma[i].vSpinShifts[0]));
      else 
        for(int i=0; i<num_lsms; i++)
          generator -> initializeEvec(i, &(sigma[i].evec[0]) );

      // generate initial occupancies
      if( isOccupancySim ) {
        printf("Master: Generating initial occupancies.\n");
        for(int i = 0; i < num_lsms; i++)
          generator->initializeOccupancies(i,&(sigma[i].occ[0]));
      }

      // if both spin and occupancy variables are needed
      // then we need to set both sets of variables before beginning calculations
      if( isSpinSim && isOccupancySim )
      for(int i = 0; i < num_lsms; i++) {
        send_evec(sigma[i],lsms_rank0[i],potentialShifter.vSpinShiftFlag);
        send_occ(sigma[i],lsms_rank0[i]);
      }
    
      std::cout << "This is the master node\n";

      // issue initial commands to all LSMS instances
      running = 0;
      bool more_work = true;
      if (total_init_steps > 0)
      {

       
        for(int i=0; i<num_lsms; i++)
        {
          std::cout << "starting initial calculation in group " << i << std::endl;

          // determine whether to send spin or occupancy change
          MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);

          if( MoveChoice == SpinMove ) {
            // dprintf("Master: Sending trial spins to Walker %d.\n", i);
            send_evec(sigma[i],lsms_rank0[i],potentialShifter.vSpinShiftFlag);
          }
          else if( MoveChoice == OccupancyMove ) {
            // dprintf("Master: Sending trial occupancies to Walker %d.\n", i);
            send_occ(sigma[i],lsms_rank0[i]);
          }
  
          num_steps--;
          running++;
          stepCount++;
          walkerSteps[i+1]++;
          if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
        }
        // first deal with the initial steps:
        while (running > 0)
        {
          if (return_moments_flag)
            MPI_Recv( &(r_values[0]), R_VALUE_OFFSET+3*size_lsms, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          else
            MPI_Recv( &(r_values[0]), R_VALUE_OFFSET, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          running--;
         
          // std::cout << "received energy E_tot =" << r_values[0] << std::endl;
          // std::cout << "    band energy E_band =" << r_values[1] << std::endl;
          printf("received energy E_tot = %25.15f\n",r_values[0]);
          printf("    band energy E_band= %25.15f\n",r_values[1]);
 
          if (total_init_steps > 0)
          {
            int r_group = (status.MPI_SOURCE - align) / comm_size;
            std::cout << "starting additional calculation in group " << r_group << std::endl;

            MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);
            if (sigma[r_group].init_steps > 0)
            {
              if( MoveChoice == SpinMove ) 
                generator -> generateUnsampledEvec(r_group, &(sigma[r_group].evec[0]), r_values[energyIndex]);
              else if( MoveChoice == OccupancyMove ) 
                generator -> generateUnsampledOcc(r_group, &(sigma[r_group].occ[0]));
              //more_work = !(generator -> generateUnsampledEvec(r_group, evecs[r_group], r_values[energyIndex]));
              sigma[r_group].init_steps;
              total_init_steps--;
            }

            // send configuration
            if( MoveChoice == SpinMove )  { 
              // dprintf("Master: Sending trial spins to Walker %d.\n", r_group);
              send_evec(sigma[r_group],lsms_rank0[r_group],potentialShifter.vSpinShiftFlag);
            }
            else if( MoveChoice == OccupancyMove ) {
              // dprintf("Master: Sending trial occupancies to Walker %d.\n", r_group);
              send_occ(sigma[r_group],lsms_rank0[r_group]);
            }
                  
            num_steps--;
            running++;
            stepCount++;
            walkerSteps[r_group+1]++;
            if (restrict_steps && num_steps <= 0) more_work = false;
            if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
            walltime = MPI_Wtime() - walltime_0;
            if (restrict_time && walltime >= max_time) more_work = false;
            if (restrict_time) std::cout << "      " << max_time - walltime << " seconds remaining\n";
          }
              
        }
      }
      more_work = true;
      running = 0;
      for (int i=0; i<num_lsms; i++)
      {
        std::cout << "starting main calculation in group " << i << std::endl;
        // select whether to perform spin or occupancy move
        MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);
        if( MoveChoice == SpinMove ) 
        {
          // dprintf("Master: Sending trial spins to Walker %d.\n", i);
          send_evec(sigma[i],lsms_rank0[i],potentialShifter.vSpinShiftFlag);
        }
        else if( MoveChoice == OccupancyMove ) 
        {
          send_occ(sigma[i],lsms_rank0[i]);
        }
        num_steps--;
        running++;
        stepCount++;
        walkerSteps[i+1]++;
        if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
      }
        
      generator -> startSampling(isSpinSim, isOccupancySim);

      // define variables to keep track whether spin (or occupancy)
      // configuration needs to be reverted to previous values
      std::vector<bool> acceptedSpinMove(num_lsms);
      std::vector<bool> acceptedOccMove(num_lsms);

      for(int i = 0; i < num_lsms; i++)
        acceptedSpinMove[i] = acceptedOccMove[i] = true;

      // wait for results and issue new commands or wind down
      while (running > 0)
      {
        accepted = false;

        MPI_Recv( &(r_values[0]), R_VALUE_OFFSET+3*size_lsms, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        running--;
        printf("received energy E_tot = %25.15f\n",r_values[0]);
        printf("    band energy E_band= %25.15f\n",r_values[1]);
        // printf("from status.MPI_SOURCE=%d\n",status.MPI_SOURCE);
        energy_accumulator += r_values[0];
        energies_accumulated++;

        // determine whether returning from a spin or occupancy move
        int r_group = (status.MPI_SOURCE - align) / comm_size;
        if( r_values[2] == 0x00 ) {
          prevMoveChoice = SpinMove;
          // dprintf("Master: Recieved energy from Walker %d after spin trial.\n",r_group);
        }
        else if( r_values[2] = 0x01 ) {
          prevMoveChoice = OccupancyMove;
          // dprintf("Master: Recieved energy from Walker %d after occupancy trial.\n",r_group);
        }
 
        if (more_work)
        {
          int r_group = (status.MPI_SOURCE - align) / comm_size;
          std::cout << "starting additional calculation in group " << r_group << std::endl;
              
          if (generator_needs_moment)
          {
            double m0, m1, m2;
            m0 = 0.0; m1 = 0.0; m2 = 0.0;
            for(int i=0; i<3*size_lsms; i+=3)
            {
              m0 += r_values[R_VALUE_OFFSET+i];
              m1 += r_values[R_VALUE_OFFSET+i+1];
              m2 += r_values[R_VALUE_OFFSET+i+2];
            }
            switch (second_dimension)
            {
              case  MagneticMoment : magnetization = std::sqrt(m0*m0 + m1*m1 + m2*m2); break;
              case  MagneticMomentX : magnetization = m0; break;
              case  MagneticMomentY : magnetization = m1; break;
              case  MagneticMomentZ : magnetization = m2; break;
            }

            // todo: separately keep track of whether last evec or occ move accepted
            // Determine if the configuration is accepted
            accepted = generator -> determineAcceptance(r_group, r_values[energyIndex], magnetization);
          } 
          else {
            // Determine if the configuration is accepted
            accepted = generator -> determineAcceptance(r_group, r_values[energyIndex]);
          }

          if(accepted)
            printf("Master: Accepted Walker %d trial move.\n",r_group);
          else
            printf("Master: Rejected Walker %d trial move.\n",r_group);

           // need to separate track of whether the last spin or occupancy change was accepted
           // this ensures the Wang-Landau generator can properly restore the old configuration
          if( prevMoveChoice == SpinMove ) 
            acceptedSpinMove[r_group] = accepted;
          else if( prevMoveChoice == OccupancyMove )
            acceptedOccMove[r_group] = accepted;

          // dprintf("Master: Updating histogram after Walker %d move.\n",r_group);
          if (potentialShifter.vSpinShiftFlag)
          {
            if (generator -> updateHistogram(r_group, &(sigma[r_group].evec[0]), &(sigma[r_group].vSpinShifts[0]), accepted))
              more_work = false;
          } else {
            if (generator -> updateHistogram(r_group, &(sigma[r_group].evec[0]), accepted))
              more_work = false;
          }
          // dprintf("Master: Updating log after Walker %d move.\n",r_group);
          generator->updateLog(r_group, &(sigma[r_group].evec[0]), &(sigma[r_group].occ[0]), 
              r_values[energyIndex], accepted, prevMoveChoice, 
              isSpinSim, isOccupancySim);

          // Prepare a new configuration
          MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);

          if( MoveChoice == SpinMove ) {
            // dprintf("Master: Generating trial spins for Walker %d.\n",r_group);
            generator -> generateEvec(r_group, &(sigma[r_group].evec[0]), acceptedSpinMove[r_group]);
            if (potentialShifter.vSpinShiftFlag)
              generator -> generatePotentialShift(r_group, &(sigma[r_group].vSpinShifts[0]), acceptedSpinMove[r_group]);
          }
          else if( MoveChoice == OccupancyMove ) {
            // dprintf("Master: Generating trial occupancies for Walker %d.\n",r_group);
            generator->generateOccupancies(r_group, &(sigma[r_group].occ[0]), acceptedOccMove[r_group]);
          }

          if( MoveChoice == SpinMove ) 
          { 
            send_evec(sigma[r_group],lsms_rank0[r_group],potentialShifter.vSpinShiftFlag);
          }
          else if( MoveChoice == OccupancyMove ) 
          {
            send_occ(sigma[r_group],lsms_rank0[r_group]);
          }
  
          num_steps--;
          running++;
          stepCount++;
          walkerSteps[r_group+1]++;
          if (restrict_steps && num_steps <= 0) more_work = false;
          if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
          walltime = MPI_Wtime() - walltime_0;
          if (restrict_time && walltime >= max_time) more_work = false;
          if (restrict_time) std::cout << "      " << max_time - walltime << " seconds remaining\n";
        }
        else
        {
          // send an exit message to this instance of LSMS
          int r_group = (status.MPI_SOURCE - align) / comm_size;

          MPI_Send( &(sigma[r_group].evec[0]), 3*size_lsms, MPI_DOUBLE, lsms_rank0[r_group], 2, MPI_COMM_WORLD);
        }

        if (step_out_flag && accepted)
        {
          step_out_file << "# iteration " << energies_accumulated << std::endl;
          step_out_file.precision(15);
          step_out_file << energies_accumulated << std::endl;
          step_out_file << r_values[0] << "  " << r_values[1] << std::endl;
          for (int j=0; j<3*size_lsms; j+=3)
            step_out_file << r_values[j+R_VALUE_OFFSET] << "  " << r_values[j+R_VALUE_OFFSET+1]
                          << "  " << r_values[j+R_VALUE_OFFSET+2] << std::endl;
        }
        // write restart file every restartWriteFrequency seconds
        if (walltime > nextWriteTime)
        {
          generator -> writeState("WLrestart.jsn");
          nextWriteTime += restartWriteFrequency;
        }

      }
      generator -> writeState("WLrestart.jsn");

      // gather statistics of energy loop counts
     long int sb[2];
     sb[0]=-1;
     sb[1]=0;
     // std::vector<long int> rb(2*(num_lsms+1));
     std::vector<long int> rb(2*size);
     MPI_Gather(sb,2,MPI_LONG,&rb[0],2,MPI_LONG,0,MPI_COMM_WORLD);
     walkerSteps[0]=stepCount;  // report the total step count for the WL master
     std::ofstream statOut; statOut.open("energyLoopCount.statistics");
       statOut<<rb[2*0]<<"   "<<walkerSteps[0]<<"   "<<rb[2*0+1]<<std::endl;
     for(int i=0; i<num_lsms; i++)
       statOut<<rb[2*lsms_rank0[i]]<<"   "<<walkerSteps[i+1]<<"   "<<rb[2*lsms_rank0[i]+1]<<std::endl;
     statOut.close();
/*
  if(evec_generation_mode==WangLandau_1d)
  (static_cast<WL1dEvecGenerator<std::mt19937> *>(generator))->writeState("WLrestart.state");
  if(evec_generation_mode==ExhaustiveIsing)
  (static_cast<ExhaustiveIsing1dEvecGenerator *>(generator))->writeState("WLrestart.state");
*/

    }

void LSMS_Hamiltonian::old_dosim()
{
  //this->parse_command_line(argc,argv);
  if( evec_random_flag ) evec_generation_mode = Random;
  restrict_steps = (num_steps>0);
  restrict_time  = (walltime>0);
  if( strcmp(step_out_name,"")!=0 )
  {
     step_out_flag = true;
     return_moments_flag = true;
  }
  this->set_modes();

  // initialize MPI:
  lsms_rank0.resize(num_lsms+1) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  world_rank = rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Construct new MPI info and compare to old
  // Not compatible with Master/Slave
  //int NWalk = num_lsms;
  //MPI_Struct mp_world = MPI_Struct::world();
  //MPI_Gang mp_gang;
  //std::cout << __FILE__ << ":" << __LINE__ << " NWalk=" << NWalk << std::endl;
  //mp_gang.init(NWalk,mp_world);
  //if( mp_gang.igang!=my_group ) std::cout << "**MPI(" << mp_world.iproc << "): mp_gang.igang=" << mp_gang.igang << " my_group=" << my_group << std::endl;

  walltime_0 = MPI_Wtime();

  // PAPI START

  write_header();

  if (generator_needs_moment) return_moments_flag = true;

  if (num_lsms <= 1)
  {
     std::cout << "num_lsms==1 known to not work" << std::endl;
     exit(1);
  }
  walkerSteps.resize(num_lsms+1); for(int i=0; i<num_lsms+1; i++) walkerSteps[i]=0;
  this->init();
  // now we get ready to do some calculations...
  if (my_group >= 0)
  {
     do_slave();
  }
  else if (world_rank == 0)
  {
     do_master();
     if (step_out_flag)
     {
       step_out_file << "# end\n-1\n"
                     << energy_accumulator / double(energies_accumulated) << std::endl;
       step_out_file.close();
     }
     std::cout << "Finished all scheduled calculations. Freeing resources.\n";
     std::cout<<"Energy mean = "<<energy_accumulator/double(energies_accumulated)<<"Ry\n";
  }
  // make sure averyone arrives here:
  MPI_Bcast (stupid, 37, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (world_rank == 0)
     MPI_Comm_free(&local_comm);
  else if (my_group >= 0)
     MPI_Comm_free(&local_comm);
  // PAPI STOP
  if (world_rank == 0)
  {
    double walltime = MPI_Wtime() - walltime_0;
    std::cout << " WL-LSMS finished in " << walltime << " seconds.\n";
    std::cout << " Monte-Carlo steps / walltime = "
              << double(stepCount) / walltime << "/sec\n";
  }
}


#if 0

// #define USE_PAPI 1
#ifdef USE_PAPI
#include <papi.h>
#endif

////////////////////////////////// PAPI START /////////////////////////////////////////////
#ifdef USE_PAPI
#define NUM_PAPI_EVENTS 4
  int hw_counters = PAPI_num_counters();
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS]; // = {PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_FP_OPS,PAPI_VEC_INS};
  char *papi_event_name[] = {"PAPI_TOT_INS","PAPI_FP_OPS",
                             "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:OP_TYPE",
                             "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:OP_TYPE"};
  // "RETIRED_INSTRUCTIONS",
  // "RETIRED_MMX_AND_FP_INSTRUCTIONS:PACKED_SSE_AND_SSE2",
  // "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:1",
  // "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:1"
  // get events from names:
  for(int i=0; i<NUM_PAPI_EVENTS; i++)
  {
    if(PAPI_event_name_to_code(papi_event_name[i],&papi_events[i]) != PAPI_OK)
    {
      // printline("Error in obtaining PAPI event code for: "+ttos(papi_event_name[i]),
      //           std::cerr,parameters.myrankWorld);
      // printline("Skipping all following events",
      //           std::cerr,parameters.myrankWorld);
      if(hw_counters>i) hw_counters=i;
    }
  }
  long long papi_values[NUM_PAPI_EVENTS+4];
  // printline("PAPI: "+ttos(hw_counters)+" counters available",std::cout,parameters.myrankWorld);
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  long long papi_real_cyc_0 = PAPI_get_real_cyc();
  long long papi_real_usec_0 = PAPI_get_real_usec();
  long long papi_virt_cyc_0 = PAPI_get_virt_cyc();
  long long papi_virt_usec_0 = PAPI_get_virt_usec();
  PAPI_start_counters(papi_events,hw_counters);
#endif

////////////////////////////////////////////// PAPI STOP ///////////////////////////////////////////////////
#ifdef USE_PAPI
  PAPI_stop_counters(papi_values,hw_counters);
  papi_values[hw_counters  ] = PAPI_get_real_cyc()-papi_real_cyc_0;
  papi_values[hw_counters+1] = PAPI_get_real_usec()-papi_real_usec_0;
  papi_values[hw_counters+2] = PAPI_get_virt_cyc()-papi_virt_cyc_0;
  papi_values[hw_counters+3] = PAPI_get_virt_usec()-papi_virt_usec_0;
  long long accumulated_counters[NUM_PAPI_EVENTS+4];
  MPI_Reduce(papi_values,accumulated_counters,hw_counters+4,
             MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  if(world_rank==0)
  {
    for(int i=0; i<hw_counters; i++)
    {
      std::cout<<"Accumulated: "<<(papi_event_name[i])<<" = "<<(accumulated_counters[i])<<"\n";
    }
    std::cout<<"PAPI accumulated real cycles : "<<(accumulated_counters[hw_counters])<<"\n";
    std::cout<<"PAPI accumulated user cycles : "<<(accumulated_counters[hw_counters+2])<<"\n";
    double gflops_papi = ((double)accumulated_counters[1])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_double = ((double)accumulated_counters[2])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_single = ((double)accumulated_counters[3])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gips = ((double)accumulated_counters[0])/(1000.0*(double)papi_values[hw_counters+1]);
    std::cout<<"PAPI_FP_OPS real GFLOP/s : "<<(gflops_papi)<<"\n";
    std::cout<<"PAPI hw double real GFLOP/s : "<<(gflops_hw_double)<<"\n";
    std::cout<<"PAPI hw single real GFLOP/s : "<<(gflops_hw_single)<<"\n";
    std::cout<<"PAPI real GINST/s : "<<(gips)<<"\n";
  }
#endif


#endif

#endif
