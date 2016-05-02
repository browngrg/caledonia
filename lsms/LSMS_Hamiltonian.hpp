#ifndef CALEDONIA_LSMS_HAMILTONIAN_HPP
#define CALEDONIA_LSMS_HAMILTONIAN_HPP


#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <iostream>
#include <fstream>

#include "SystemParameters.hpp"

//#include "EvecGenerator.h"
// Header for the evec generator class to be used in lsms.cc
// also example implementations:
// ConstantEvecGenerator
// RandomEvecGenerator

#ifndef EVEC_GENERATOR_H
#define EVEC_GENERATOR_H

#include <stdio.h>
#include <vector>
#include "random_evec.h"

typedef enum {SpinMove, OccupancyMove} MoveChoice_t;

/*
YingWai's note  (Dec 5, 13):
1. The split of the original generateEvec into "determineAcceptance, 
   updateHistogram and generateEvec" is only implemented in 
   WangLandau.h and WangLandau2d.h
2. generateEvec is renamed to updateHistogram in ExhautiveIsing.h.
   Nothing has been changed there otherwise.
   Main program (wl_lsms.cpp) was changed in a way where ExhaustiveIsing 
   was NOT considered. Use with caution!
*/

class EvecGenerator
{
 public:

  virtual double getDOSRatio(int instance, double energy)
  { return 0.0; }

// YingWai: do we really need the function overloading for determineAcceptance?   (Dec 4, 13)
  virtual bool determineAcceptance(int instance, double energy)
  { return false; }

  virtual bool determineAcceptance(int instance, double energy, double magnetization)
  {
    return determineAcceptance(instance, energy);
  }

  virtual bool updateHistogram(int instance, double *evecs, bool accepted)
  { return false; }
/*
  virtual bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted)
  {
    return updateHistogram(instance, evecs, accepted);
  }
*/
  virtual bool updateHistogram(int instance, double *evecs, double *potentialShifts, bool accepted)
  {
    return updateHistogram(instance, evecs, accepted);
  }

  virtual void writeDOS(const char* fname) { std::cout << __FILE__ << ":" << __LINE__ << " BASE::writeDOS" << std::endl; return; }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, int check)
  { return false; }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, int check)
  { 
    return updateHistogramFromRE(instance, evecs, energy, check);
  }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double *potentialShifts, int check)
  {
    return updateHistogramFromRE(instance, evecs, energy, check);
  }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, double *potentialShifts, int check)
  {
    return updateHistogramFromRE(instance, evecs, energy, check);
  }
/*
// YingWai: keep these function overloadings for updateHistogram for the moment,
//          but might want to get rid of them later when code is stable  (Dec 4, 13)

  virtual bool updateHistogram(int instance, double *evecs, double energy)
  { return false; }
  
  virtual bool updateHistogram(int instance, double *evecs, double energy, bool *accepted)
  {
    *accepted = true;
    return updateHistogram(instance, evecs, energy);
  }

  virtual bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted)
  {
    *accepted = true;
    return updateHistogram(instance, evecs, energy);
  }
*/

  virtual void generatePotentialShift(int instance, double *potentialShifts, bool accepted) {};

  virtual void generateEvec(int instance, double *evecs, bool accepted) {};

  virtual void initializeEvecAndPotentialShift(int instance, double *evecs, double *potentialShifts) {};

  virtual void initializeEvec(int instance, double *evecs)
  { 
    //generateEvec(instance, evecs, 0.0); 
    bool accepted {false};
    generateEvec(instance, evecs, accepted); 
  }

  virtual void generateUnsampledEvec(int instance, double *evecs, double energy) 
  { 
    //generateEvec(instance, evecs, energy);
    bool accepted {false};
    generateEvec(instance, evecs, accepted); 
  }

  virtual void startSampling(bool isspin = true, bool isocc = false) {;}

  virtual void writeState(const char *name) {;}

  void setVerbosity(int v) { verbosity = v; }

  int verbosity;

  // dummy functions to be used for occupancy variables
  // these routines were added to the Evec class because they share Monte Carlo logic
  virtual MoveChoice_t selectMoveType(bool isSpinSim, bool isOccSim)  { return SpinMove; }
  virtual void setAlloyClasses(const AlloyMixingDesc&, int* siteclass) {}
  virtual void generateOccupancies(int instance, int *occ, bool acceptedOcc) {}
  virtual void generateUnsampledOcc(int inst, int *occ) {}
  virtual void initializeOccupancies(int inst, int *occ) {}
  virtual void updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
      bool isspin, bool isocc) {}
};


class ConstantEvecGenerator : public EvecGenerator
{
 public:

  ConstantEvecGenerator(size_t num_spins, double evec[3], int nw = 0)
  {
    printf("Constructor of ConstantEvecGenerator called.\n");
    n_spins = num_spins;
    setEvec(evec);
    n_walker = nw;
    if (nw > 0)
    {
      walker_step.resize(nw);
      for (int i=0; i<nw; i++) walker_step[i] = 0;
    }
  }

  void setEvec(double evec[3])
  {
    ev[0] = evec[0];
    ev[1] = evec[1];
    ev[2] = evec[2];
  }

  void generatePotentialShift(int instance, double *potentialShifts, bool accepted)
  {
    for (size_t i=0; i<n_spins; i++)
      potentialShifts[i] = 0.0;
  }

  void generateEvec(int instance, double *evecs, bool accepted)
  {
    for (size_t j=0; j<3*n_spins; j+=3)
    {
      evecs[j] = ev[0];
      evecs[j+1] = ev[1];
      evecs[j+2] = ev[2];
    }
    if (instance < n_walker) walker_step[instance]++;
    //return false;
  }

  void initializeEvecAndPotentialShift(int instance, double *evecs, double *potentialShifts)
  {
    initializeEvec(instance, evecs);
    generatePotentialShift(instance, potentialShifts, true);
  }

  void initializeEvec(int instance, double *evecs)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
    {
      evecs[j] = ev[0];
      evecs[j+1] = ev[1];
      evecs[j+2] = ev[2];
    }
  }

 private:

  size_t n_spins;
  double ev[3];
  int n_walker;
  std::vector<int> walker_step;

};


class RandomEvecGenerator : public EvecGenerator
{
 public:

  RandomEvecGenerator(size_t num_spins)
  { n_spins = num_spins; }

  void generateEvec(int inst, double *evecs, double energy)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
      random_evec(&evecs[j]);
    //return false;
  }

  void initializeEvec(int instance, double *evecs)
  { generateEvec(instance, evecs, 0.0); }

 private:
  size_t n_spins;

};

#endif




//#include "WangLandau.h"
#ifndef LSMS_WANG_LANDAU_H
#define LSMS_WANG_LANDAU_H

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// we use BOOST for the random number generator
// #include <boost/random.hpp>
#include <random>
#include "../../mjson/json.h"
//#include "EvecGenerator.h"
#include "Graph1dMoments.hpp"
#include "../Potential/PotentialShifter.hpp"

void inline performGlobalUpdate(Graph1dMoments<double,double> &g, double kappa, double lambda, double omega)
{
  for(int i=0; i<g.getN(); i++)
    if(g[i]>omega) g[i]+=kappa*std::exp(-lambda/(g[i]-omega));
}

class StatesWriter
{
public:
  StatesWriter(const char *filename=NULL)
  {
    if(filename==NULL) writeFlag=false;
    else {writeFlag=true; of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(17);}
  }
  ~StatesWriter() {if(writeFlag) of.close();}
  void writeHeader(double lnF, int numWalk, int numSpin, double **spins, 
      bool isSpinSim = true, bool isOccSim = false, int **occup = NULL)
  {

    if( spins==NULL ) return;   // browngrg 10/14/2015 to allow evec->LSMS_Config
    if( !writeFlag ) return;

    if( isSpinSim ) {
      of << lnF << " " << numWalk << " " << numSpin << std::endl;
      for(int i = 0; i < numWalk; i++) {
        of << i;
        for(int j = 0; j < 3*numSpin; j++) 
          of << " " << spins[i][j];
        of << std::endl;
      }
    }

    if( isOccSim ) return;
  }
  void writeChange(int iWalk, int numRet, int ispin, double *ev, double E)
  {
    if(writeFlag)
      of<<iWalk<<" "<<numRet<<" "<<ispin<<" "<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<" "<<E<<std::endl;
  }
  void writeChangeOcc(int iWalk, int numRet, int i, int j, int occ_i, int occ_j, double E)
  {
    if(!writeFlag) return;
      of<<iWalk<<" "<<numRet<<" ("<<i<<"<-->"<<j<<") "<<occ_i<<" "<< occ_j <<" "<<E<<std::endl;
  }
  void newFile(const char *filename=NULL)
  {
    if(writeFlag) of.close();
    writeFlag=false;
    if(filename!=NULL){
      writeFlag=true;
      of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);
    }
  }
private:
  bool writeFlag;
  std::ofstream of;
};

template<class RNG = std::mt19937>
// template<class RNG = boost::mt19937>
class WL1dEvecGenerator : public EvecGenerator
{
 public:
  WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                    PotentialShifter &potentialShifter,
                    const char *init_file_name=NULL, const char *out_file_name=NULL,
                    const char *states_filename=NULL, int ** occ_p = NULL);

  bool determineAcceptance(int instance, double energy);
  bool determineAcceptance(int instance, double energy, double magnetization);

  bool updateHistogram(int instance, double *evecs, bool accepted);
/*
  bool updateHistogram(int instance, double *evecs, double energy, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy) {bool h; return updateHistogram(instance, evecs, energy, &h);}
  bool updateHistogram(int instance, double *evecs) {std::cerr<<"Need energy for WL1dEvecGenerator\n"; exit(1);}
*/

  void generateEvec(int instance, double *evecs, bool accepted);
  //void generateEvec(int instance, double *evecs, double energy);
  void generatePotentialShift(int instance, double *potentialShifts, bool accepted);

  void initializeEvecAndPotentialShift(int inst, double *evecs, double *potentialShifts);
  void initializeEvec(int instance, double *evecs);

  void generateUnsampledEvec(int instance, double *evecs, double energy)
  {
    initializeEvec(instance, evecs); 
    //return false;
  }

  void startSampling(bool isSpinSim = true, bool isOccSim = false)
  { sw.writeHeader(gamma, n_walkers, n_spins, evecs_pointer, isSpinSim, isOccSim, occupancy_ptr); }

  void writeState(const char *name);
  void writeDOS(const char *name);

  // Wang-Landau for occupancy variables
  MoveChoice_t selectMoveType(bool isSpinSim, bool isOccSim);
  void setAlloyClasses(const AlloyMixingDesc&, int* siteclass);
  void generateOccupancies(int instance, int *occ, bool acceptedOcc);
  void generateUnsampledOcc(int inst, int *occ);
  void initializeOccupancies(int inst, int *occ);
  void updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
      bool isspin, bool isocc);

 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;
  int n_initialized_from_file_occ;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
  // boost::uniform_real<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd01; //(0.0,1.0);

  /*
  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  */

  double hMinimum;
  Graph1dMoments<double,double> dos, histo;
  Kernel1d<double,double> dosKernel, histoKernel, nullKernel;
  KernelType kernelType;

  unsigned long accept, reject, acceptSinceLastChange;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  bool fixEnergyWindow;
  
  // instance specific:
  std::vector<long> ref0, ref1;
  std::vector<double> position, magnetizationAtPosition;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<int> lastChangePotentialShift;
  std::vector<int> lastAccepted;
  std::vector<int> lastAcceptedPotentialShiftIndex;
  std::vector<double> lastAcceptedEnergy;
  std::vector<double> lastAcceptedEvec;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]
  std::vector<double> oldPotentialShift;
  std::vector<double> lastAcceptedPotentialShift;

  std::vector<int> numRetentions;

  // changes to accomodate alloying (i.e. variable site occupancies)
  std::vector< std::pair<int,int> > lastSwapOcc;
  std::vector< std::pair<int,int> > lastAcceptedSwap;
  std::vector< std::pair<int,int> > lastAcceptedOcc;
  AlloyMixingDesc alloyDesc;
  std::vector<int> siteAlloyClass;
  std::vector<int> numSites_per_AlloyClass;
  int ** occupancy_ptr;

  char *statesFile;
  StatesWriter sw;

  int changeMode;
  bool histogramUpdateMode;
  int updatesPerBin;

  struct {double kappa, lambda, omega; int frequency, changes;} globalUpdate;

#ifdef ISING
    void inline random_evec_1(double ev[3])
  {
    ev[0]=ev[1]=0.0;
    ev[2]=1.0;
    if(rng()%2 == 0) ev[2]=-ev[2];
  }
#else
  void inline random_evec_1(double ev[3])
  {
    double x,y,z;
    do {
      x = rnd(rng);
      y = rnd(rng);
    } while(x*x+y*y>1);
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r;
    y *= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z;
    r=1.0/sqrt(x*x+y*y+z*z);
    ev[0]=x*r; ev[1]=y*r; ev[2]=z*r;
  }
#endif

#ifdef ISING    
  void inline random_evec(double ev[3])
  {
    ev[2]=-ev[2];
  }
#else
  void inline random_evec(double ev[3])
  {
    double x, y, z;
    do {
      x = rnd(rng); y = rnd(rng);
    } while(x*x+y*y>1); 
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r; y*= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z; 
    // Project out the parallel component;
    r = x*ev[0] + y*ev[1] + z*ev[2];
    x -= r*ev[0]; y -= r*ev[1]; z -= r*ev[2];
    r = x*x + y*y + z*z;
    double t = 1-0.3*rnd(rng);
    ev[0] *= t; ev[1] *= t; ev[2] *= t;
    r = sqrt((1-t*t)/r);
    ev[0] += x*r; ev[1] += y*r; ev[2] += z*r;
    r=1.0/sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    ev[0]*=r; ev[1]*=r; ev[2]*=r;
    
    /*  
    ev[2]=1.0-2.0*rnd(rng);
    // ev[2]=rnd11(rng);
    double phi=2.0*M_PI*rnd(rng);
    // double phi=rnd0pi(rng);
    double cos_theta=sqrt(1-ev[2]*ev[2]);
    ev[0]=cos_theta*cos(phi);
    ev[1]=cos_theta*sin(phi);
    */  
  }
#endif

  double minVxShift {0.0}, maxVxShift {0.0}, rangeVxShift {0.0};

  double inline randomPotentialShift()
  {
    return minVxShift + rangeVxShift * rnd01(rng);
  }
};


template<class RNG>
WL1dEvecGenerator<RNG>::WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                                          PotentialShifter &potentialShifter,
            const char *init_file_name, const char *out_file_name, 
                                          const char *states_filename, int **occ_p)
: sw(states_filename)
{

  std::cout << __FILE__ << ":" << __LINE__ << std::endl;

  if (potentialShifter.vSpinShiftFlag) {
    minVxShift = potentialShifter.minShift;
    maxVxShift = potentialShifter.maxShift;
    rangeVxShift = maxVxShift - minVxShift;
  }

  verbosity=0;

  changeMode = 4;

  globalUpdate.frequency=0;
  globalUpdate.changes=0;
  globalUpdate.kappa=1.0;
  globalUpdate.lambda=1.0;
  globalUpdate.omega=0.5;

  long nX=-1;
  histogramUpdateMode=false;
  updatesPerBin=100;
  double interval=0.01;
  double kernelWidth=0.1;
  double xMin=-std::numeric_limits<double>::max();// double xMin=-HUGE;
  double xMax=1.0;
  bool readPositions=false;
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  n_initialized_from_file_occ = 0;
  evecs_pointer = ev_p;
  occupancy_ptr = occ_p;
  ref0.resize(n_walkers);  for(int i=0; i<n_walkers; i++) ref0[i]=-1; // ref0[i]=HUGE;
  ref1.resize(n_walkers);
  position.resize(n_walkers);
  magnetizationAtPosition.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);

  oldPotentialShift.resize(n_walkers);

  lastChangePotentialShift.resize(n_walkers);
  lastAcceptedPotentialShiftIndex.resize(n_walkers);
  lastAcceptedPotentialShift.resize(n_walkers);

  lastSwapOcc.resize(n_walkers);
  lastAcceptedSwap.resize(n_walkers);
  lastAcceptedOcc.resize(n_walkers);
  for(int i=0; i<n_walkers; i++)
    lastSwapOcc[i].first = -1;

  fixEnergyWindow = false;
  
  // fprintf(stderr, "WARNING: Fixing random number seed!\n");
  // rng.seed(1);

  for(int i=0; i<n_walkers; i++)
  {
    lastAccepted[i]=-2;
    lastAcceptedPotentialShiftIndex[i] = -2;
    lastAcceptedPotentialShift[i] = 0.0;
  }
  for(int i=0; i<3*n_walkers; i++)  lastAcceptedEvec[i]=0.0;

  // browngrg 10/14/2015 
  // This is done in the do_master() loop
  //for(int i=0; i<n_walkers; i++) initializeEvec(i,evecs_pointer[i]);

  oldSpin.resize(3*n_walkers);

  statesFile=NULL;
  if(states_filename!=NULL)
    {
      statesFile=(char *)malloc(sizeof(char)*(1+strlen(states_filename)));
      strcpy(statesFile,states_filename);
    }

  numRetentions.resize(n_walkers);
  for(int i=0; i<n_walkers; i++) numRetentions[i]=0;

  /*
  nX = -1;
  xMin = -HUGE; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  */
  dos_out_name="";//"dos1d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  cycleCount=0;
  hMinimum= 1;     //   10
  acceptSinceLastChange=accept=reject=0;
  gammaFinal=1.e-6;
  flipPerUpdate=1; //  100
  updateCycle= 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.75;
  updatesPerBin =100;

  kernelType=None;

  // special processing flags:
  int clearHistogram = 0;
  int setFirstWalkerToFM = 0;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(interval);
  histo.setDeltaAndClear(interval);

  if(init_file_name!=NULL && init_file_name[0]!=0)
  {
    std::string label, value;

    std::cout<<"Reading Wang-Landau configuration from: "<<init_file_name<<std::endl;

    dos_out_name=std::string(init_file_name)+".out";
    if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
    std::ifstream inp(init_file_name);
    std::ostringstream buff;

    std::string line;
    while(std::getline(inp,line)) 
      buff << line << std::endl;

    inp.close();

    std::string fileString = buff.str();
    const char* fileChars  = fileString.c_str();
    json_t *json_root=NULL;

    json_parse_document(&json_root, (char *)fileChars);

    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    if(json_root == NULL || json_root->type != JSON_OBJECT)
    {
      std::ostringstream message;
      std::cerr << "In WL1dEvecGenerator(" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;
      if(label=="xMin") xMin=atof(it->child->text);
      else if(label=="xMax") xMax=atof(it->child->text);
      else if(label=="fixEnergyWindow")
      {
        if( atoi(it->child->text) != 0 ) fixEnergyWindow = true;
      } 
      else if(label=="interval") interval=atof(it->child->text);
      else if(label=="kernelWidth") kernelWidth=atof(it->child->text);
      else if(label=="kernelType")
      {
        std::string strValue(it->child->text);
        kernelType=getKernelType(strValue);
      }
      else if(label=="gamma") gamma=atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal=atof(it->child->text);
      else if(label=="nX")
      {
        nX=atoi(it->child->text);
  dos.setRangeAndClear(xMin,xMax,nX);
  histo.setRangeAndClear(xMin,xMax,nX);
  /*
        if(dos!=NULL) free(dos);
        if(histo!=NULL) free(histo);
        dos=(double *)calloc(nX,sizeof(double));
        histo=(int *)calloc(nX,sizeof(int));
  */
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="cycleCount") cycleCount=atoi(it->child->text);
      else if(label=="changeMode") changeMode=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
      else if(label=="histogramMinimum") hMinimum=atof(it->child->text);
      else if(label=="updatesPerBin") updatesPerBin=atoi(it->child->text);
      else if(label=="globalUpdate.frequency") globalUpdate.frequency=atoi(it->child->text);
      else if(label=="globalUpdate.changes") globalUpdate.changes=atoi(it->child->text);
      else if(label=="globalUpdate.kappa") globalUpdate.kappa=atof(it->child->text);
      else if(label=="globalUpdate.lambda") globalUpdate.lambda=atof(it->child->text);
      else if(label=="globalUpdate.omega") globalUpdate.omega=atof(it->child->text);
      else if(label=="seed") rng.seed(atoi(it->child->text));
      else if(label=="accept") accept=atol(it->child->text);
      else if(label=="acceptSinceLastChange") acceptSinceLastChange=atol(it->child->text);
      else if(label=="reject") reject=atol(it->child->text);
//*
      else if(label=="rngState")
      {
        std::string strValue(it->child->text);
        std::stringstream strStream(strValue, std::stringstream::in);
        strStream>>rng;
      }
//*/
      else if(label=="dos")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          dos[j++]=atof(i->text);
        }
        if(j!=dos.getN()) {std::cout<<"ERROR #(dos) "<<j<<" != nX "<<dos.getN()<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          histo[j++]=atof(i->text);
        }
        if(j!=histo.getN()) {std::cout<<"ERROR #(histo) "<<j<<" != nX "<<histo.getN()<<std::endl; exit(1);}
      }
      else if(label=="moments")
      {
        json_t *a = it->child->child;
        int k = atoi(a->text);
        dos.setNumberOfMoments(k);
        a=a->next;
        int jj=0;
        for(json_t *j=a->child; j!=NULL; j=j->next)
        {
          dos.setNumberOfSamplesAtIdx(jj,atoi(j->text));
          jj++;
        }
        json_t *i=a->next;
        for(int kk=0; kk<k; kk++)
        {
          int jj=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
          {
            dos.setMomentAtIdx(jj,kk,atof(j->text));
            jj++;
          }
          i=i->next;
        }
      }
      else if(label=="moments.k")
      {
        dos.setNumberOfMoments(atoi(it->child->text));
      }
      // additional read for occupancy scheme
      // mirrors 'evecs' read
      else if( label == "occupancies" ) {
        json_t *a = it->child;
        n_initialized_from_file_occ=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file_occ<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              occupancy_ptr[n_initialized_from_file_occ][k++]=atoi(j->text);
            }
            // initialize lastSwapOcc to -1 to indicate no previous swap
            lastSwapOcc[n_initialized_from_file_occ].first = -1;
           n_initialized_from_file_occ++;
           if(k != n_spins) {std::cout<<"ERROR #(occs) "<<k<<" != n_spins "<<n_spins<<std::endl; exit(1);}
          }
        } 
      }
      else if( label == "lastSwapOcc" ) {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              if( k == 0 ) lastSwapOcc[n_initialized].first  = atoi(j->text);
              if( k == 1 ) lastSwapOcc[n_initialized].second = atoi(j->text);
              k++;
            }
            n_initialized++;
            if(k!=2) {std::cout<<"ERROR #(lastSwapOcc) "<<k<<" != 2"<<std::endl; exit(1);}
          }
        } 
      }
      //
      else if(label=="evecs")
      {
// browngrg 10/14/2015
if(evecs_pointer==NULL)
{
   std::cout << __FILE__ << ":" << __LINE__ << " Can't read evecs from json file" << std::endl;
   continue;
}
        json_t *a = it->child;
        n_initialized_from_file=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              evecs_pointer[n_initialized_from_file][k++]=atof(j->text);
            }
// initialize oldSpin and lastChange to point to site 0
            lastChange[n_initialized_from_file]=0;
            oldSpin[  3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][0];
            oldSpin[1+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][1];
            oldSpin[2+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][2];
//
            n_initialized_from_file++;
            if(k!=3*n_spins) {std::cout<<"ERROR #(evecs) "<<k<<" != 3*n_spins "<<3*n_spins<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="oldSpin")
      {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              oldSpin[n_initialized*3+k++]=atof(j->text);
            }
            n_initialized++;
            if(k!=3) {std::cout<<"ERROR #(oldSpin) "<<k<<" != 3"<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="lastChange")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers) lastChange[j++]=atoi(i->text);
          n_initialized++;
        }
      }
      else if(label=="position")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            position[j++]=atof(i->text);
            n_initialized++;
          }
        }
        readPositions=true;
      }
      else if(label=="magnetizationAtPosition")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            magnetizationAtPosition[j++]=atof(i->text);
            n_initialized++;
          }
        }
//        readMagnetizationAtPositions=true;
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else if(label=="clearHistogram") clearHistogram=atoi(it->child->text);
      else if(label=="setFirstWalkerToFM") setFirstWalkerToFM=atoi(it->child->text);
      else std::cout << __FILE__ << ":" << __LINE__ << " WARNING: unknown label: " << label << std::endl;
    }

    // set xMax and xMin or interval depending on nX:
    if(nX>0)
    {
      interval=(xMax-xMin)/double(nX);
      dos.setRange(xMin,xMax); histo.setRange(xMin,xMax);
    } else {
      dos.setDeltaAndClear(interval); histo.setDeltaAndClear(interval);
    }

    json_free_value(&json_root);
  }

  dosKernel.setWidthAndClear(interval,kernelWidth);
  histoKernel.setWidthAndClear(interval,kernelWidth);
  nullKernel.setWidthAndClear(interval,kernelWidth);
  initKernel(kernelType,dosKernel);
  dosKernel.scale(gamma/dosKernel(0.0));
  initKernel(kernelType,histoKernel);
  histoKernel.scale(1.0/histoKernel(0.0));
  initKernel(kernelType,nullKernel);
  nullKernel.scale(0.0);
//  histoKernel.scale(interval);

  if(readPositions)
    for(int i=0; i<n_walkers; i++)
    {
      ref0[i]=dos.idx(position[i]);
      lastAccepted[i] = -1;
      lastAcceptedEnergy[i] = position[i];
    }
      
  if(clearHistogram!=0)
  {
    std::cout<<"Clearing Wang-Landau Histogram.\n";
    histo.clear();
  }

  if(setFirstWalkerToFM!=0)
  { if(setFirstWalkerToFM<0)
    {
      std::cout<<"Setting all walkers to FM.\n";
      setFirstWalkerToFM=n_walkers;
    } else {
      std::cout<<"Setting first "
               << std::min(setFirstWalkerToFM,n_walkers)
               << " walkers to FM.\n";
    }
// browngrg 10/14/2015
/*
    for(int n=0; n<std::min(setFirstWalkerToFM,n_walkers); n++)
      for(int i=0; i<n_spins; i++)
      {
        evecs_pointer[n][  3*i]=0.0;
        evecs_pointer[n][1+3*i]=0.0;
        evecs_pointer[n][2+3*i]=1.0;
      }
*/
// initialize oldSpin and lastChange to point to site 0
   lastChange[0]=0;
// browngrg 10/14/2015
/*
   oldSpin[0]=evecs_pointer[0][0];
   oldSpin[1]=evecs_pointer[0][1];
   oldSpin[2]=evecs_pointer[0][2];
*/
// initialize oldPotentialShift and lastChangePotentialShift to point to site 0
    lastChangePotentialShift[0] = 0;
    oldPotentialShift[0] = 0.0;
  }

  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout << __FILE__ << ":" << __LINE__ << " Wang-Landau output will be written to: \"" << dos_out_name << "\"" << std::endl;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvecAndPotentialShift(int inst, double *evecs, double *potentialShifts)
{
  for (size_t j=0; j<n_spins; j++)
    potentialShifts[j] = 0.0;

  initializeEvec(inst, evecs);
}

template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  bool firstFerromagnetic=true;
  if(inst>=n_initialized_from_file)
    {
      if(firstFerromagnetic)
        if(inst==0)
    for(size_t j=0; j<3*n_spins; j+=3)
      {evecs[j]=0.0; evecs[j+1]=0.0; evecs[j+2]=1.0;}
        else if(inst==1)
          for(size_t j=0; j<3*n_spins; j+=3)
            {evecs[j+2]= (j/3)%2 ? 1.0 : -1.0 ; evecs[j+1]=0.0; evecs[j]=0.0;}
        else
    for(size_t j=0; j<3*n_spins; j+=3)
      random_evec_1(&evecs[j]);
      else
        for(size_t j=0; j<3*n_spins; j+=3)
          random_evec_1(&evecs[j]);
    }

  out[inst]=false;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::writeState(const char* name)
{
// browngrg 10/14/2015
return;
  // if(!syncronizeGraphs(dos,histo))
  if(dos.getN()!=histo.getN())
  {
    std::cout<<"Histogramm size dosn't match DOS! Clearing histogramm!\n";
    histo.setRangeAndClear(dos.getMinX(),dos.getMaxX(),dos.getN());
  }
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(15);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    if( fixEnergyWindow )
      ofile<<"\"fixEnergyWindow\" : 1,\n";
    ofile<<"\"nX\" : " << dos.getN() << ",\n";
    if(kernelType!=None)
      ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma << ",\n";
    ofile<<"\"accept\" : " << accept <<",\n";
    ofile<<"\"acceptSinceLastChange\" : " << acceptSinceLastChange <<",\n";
    ofile<<"\"reject\" : " << reject <<",\n";
    ofile<<"\"gammaFinal\" : " << gammaFinal << ",\n";

    ofile<<"\"changeMode\" : "<< changeMode <<",\n";
    ofile<<"\"flatnessCriterion\" : "<< flatnessCriterion <<",\n";
    ofile<<"\"histogramMinimum\" : "<< hMinimum <<",\n";
    ofile<<"\"updatesPerBin\" : "<< updatesPerBin <<",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << ",\n";
    ofile<<"\"updateCycle\" : " << updateCycle << ",\n";

    ofile<<"\"globalUpdate.frequency\" : "<<  globalUpdate.frequency << ",\n";
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"globalUpdate.kappa\" : "<<  globalUpdate.kappa << ",\n";
    ofile<<"\"globalUpdate.lambda\" : "<<  globalUpdate.lambda << ",\n";
    ofile<<"\"globalUpdate.omega\" : "<<  globalUpdate.omega << ",\n";

    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"rngState\" : \""<<rng<<"\",\n";
    // additional output for occupancies
    ofile << "\"occupancies\" : [\n";
    for(int i = 0; i < n_walkers; i++) {
      ofile << "[ ";
      int j = 0, stride = 20;
      while( j < n_spins ) {
        ofile << occupancy_ptr[i][j];
        if( j+1 != n_spins ) ofile << ", ";
        if( j+1 % stride == 0 ) ofile << "\n";
        j++;
      }
      ofile << " ]";
      if( i+1 != n_walkers ) ofile << ", ";
      ofile << "\n";
    }
    ofile << "],\n";
    ofile << "\"lastSwapOcc\" : [\n";
    for(int i = 0; i < n_walkers; i++) {
      ofile << "[ ";
      ofile << lastSwapOcc[i].first;
      ofile << ", ";
      ofile << lastSwapOcc[i].second;
      ofile << " ]";
      if( i+1 != n_walkers ) ofile << ", ";
      ofile << "\n";
    }
    ofile << "],\n";
    //
    ofile<<"\"evecs\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[\n";
      for(int j=0; j<3*n_spins; j+=3)
        ofile<<evecs_pointer[i][j]<<", "<<evecs_pointer[i][j+1]
             <<", "<<evecs_pointer[i][j+2]<<((j==3*n_spins-3)?"\n":",\n");
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"oldSpin\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[ "<<oldSpin[3*i]<<", "<<oldSpin[3*i+1]
             <<", "<<oldSpin[3*i+2];
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"lastChange\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< lastChange[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"position\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< position[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"magnetizationAtPosition\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< magnetizationAtPosition[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"stepsSinceLastHistogramUpdate\" : " << stepsSinceLastHistogramUpdate << ",\n";
    ofile<<"\"numberOfUpdatesSinceLastBoost\" : " << numberOfUpdatesSinceLastBoost << ",\n";
    ofile<<"\"modificationFactorChanges\" : " << modificationFactorChanges << ",\n";
    ofile<<"\"cycleCount\" : " << cycleCount << "\n";
    ofile<<"}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
} 

template<class RNG>
void WL1dEvecGenerator<RNG>::writeDOS(const char* name)
{
  std::ofstream ofile(name);
  // write dos;
  // we are using JSON as our file format
  
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(15);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << "," <<std::endl;
    ofile<<"\"xMax\" : " << dos.getMaxX() << "," <<std::endl;
    ofile<<"\"nX\" : " << dos.getN() << "," <<std::endl;
    ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma <<"," <<std::endl;
    ofile<<"\"gammaFinal\" : " << gammaFinal << "," <<std::endl;
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << "," <<std::endl;
    ofile<<"\"updateCycle\" : " << updateCycle << "," <<std::endl;
    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"]\n}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy)
{
  return determineAcceptance(instance, energy, 0.0);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy, double magnetization)
{

  bool accept_step;
  // +++++++++ Added by Odbadrakh and Don on Aug 11, 2010
  // +++ If a flip of spin results in energy in the same bin as previous accepted 
  // +++ state, it is accepted all the time. We change it by accepting if the
  // +++ new state has lower density of states according to a linear spline of DOS
  
  double dos_differ;
  double energy_differ;
  double to_go_or_not;

  // +++++++++ end of modification. Follow the new variables to see the actual changes 
  // energy between xMin and xMax ?
  int grid = dos.idx(energy);

  // dos.test();

  stepsSinceLastHistogramUpdate++; // counter[0]++

  // record initial energy for step out
  if (lastAccepted[instance] == -2)
  {
    position[instance] = lastAcceptedEnergy[instance] = energy;
    lastAccepted[instance] = -1;
  }
  
  if (grid < 0 || grid > dos.getN()-1)
  {

    // guard for fixed window
    if( fixEnergyWindow ) {

      // default action is reject
      accept_step = false;

      // not in window yet
      if( ref0[instance] < 0 || ref0[instance] > dos.getN()-1 ) { 
        out[instance] = true;
        if( grid < 0 && energy - position[instance] > 0.0) accept_step = true;
        if( grid > dos.getN()-1 && energy - position[instance] < 0.0) accept_step = true;
      }

      // already inside window
      else
        out[instance] = false;

      return accept_step;
    }

    dos.extendTo(energy - dosKernel.getWidth()); histo.extendTo(energy - histoKernel.getWidth());
    dos.extendTo(energy + dosKernel.getWidth()); histo.extendTo(energy + histoKernel.getWidth());
    grid = dos.idx(energy);
    // we need to adjust the grid point of the walker positions (ref0) for all walkers
    for (long ii=0; ii<n_walkers; ii++)
      if (ref0[ii] >= 0) ref0[ii] = dos.idx(position[ii]);
  }
  numRetentions[instance]++;
  out[instance] = false;
  ref1[instance] = grid;
  if (ref0[instance] < 0 || ref0[instance] > dos.getN()-1 )
  {
      accept_step = true;
      ref0[instance] = ref1[instance];
      position[instance] = energy;
      magnetizationAtPosition[instance] = magnetization;
  }
  else 
  {
    if (abs(ref1[instance] - ref0[instance]) < 1 ) 
    {
// Actual change made by Odbadrakh, Aug 30, 2010
      dos_differ = dos[ref0[instance]-1] - dos[ref0[instance]];
      energy_differ = energy - lastAcceptedEnergy[instance];
      to_go_or_not = dos_differ * energy_differ;
      if (to_go_or_not >= 0.0)        // accepts all downhill changes
      {
        accept_step = true;
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance] = magnetization;
      } 
      else       // uphill moves
      {
        if(rnd(rng) < exp(to_go_or_not/dos.getDelta()))
        { //std::cout<<" dos "<<instance<<"  weight "<<dos.getDelta()<<std::endl;
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance] = magnetization;
        }
        else
        {
          accept_step = false;
        }
      }
    }
    else
    { 
      if (dos[ref1[instance]] <= dos[ref0[instance]])
      {
        accept_step = true;
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance] = magnetization;
      }
      else
      {
        if(rnd(rng) < exp(dos[ref0[instance]] - dos[ref1[instance]]))
        {
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance] = magnetization;
        }
        else 
        {
          accept_step = false;
        }
      }
      
    }

 }

// End of change made on Aug 30, 2010
  if (verbosity > 2)
    std::cout << "WangLandau 1d EvecGenerator step "
              << modificationFactorChanges << ":" << numberOfUpdatesSinceLastBoost << ":"
              << stepsSinceLastHistogramUpdate << " nX=" << dos.getN()
        << " [" << dos.getMinX() << ", " << dos.getMaxX() << "] "
              << (accept_step ? "accepted" : "rejected")
              << " Energy = " << energy << ", Instance " << instance << std::endl;

  return accept_step;

}


/*
template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, double energy, bool *accepted)
{
  return updateHistogram(instance, evecs, energy, 0.0, accepted);
}
*/

template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, bool accepted)
{

// Update histogram and DOS
  if(stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      // addKernel(dos,dosKernel,energy);
      if (!histogramUpdateMode)             // standard Wang-Landau
      {
        addKernel(dos, dosKernel, position[instance]);
        dos.addMomentsAtIdx(dos.idx(position[instance]), magnetizationAtPosition[instance]);
        // if(accepted) addKernel(histo,histoKernel,position[instance]);
        addKernel(histo, histoKernel, position[instance]);
        // addKernel(dos, dosKernel, energy);
        // addKernel(histo, histoKernel, energy);
      }
      else
      {
        // addKernel(dos, nullKernel, position[instance]);
        if(accepted) addKernel(histo, histoKernel, position[instance]);
      }
    } 
    else
    {
      if(!fixEnergyWindow) {
        std::cerr << "ATTENTION: We should never reach this place in WL1dEvecGenerator!\n";
        exit(1);
      }
    }
  }

// 1/t algorithm, change gamma at every step
// Reference: Phys. Rev. E 75, 046701 (2007).
  if(changeMode & (8+16+32))
  {
    long n = accept + reject;
    double dn;

    if ((changeMode &  (8+16)) == 8) n = accept;
    else if ((changeMode &  (8+16)) == 16) n = reject;
    else if ((changeMode &  (8+16)) == 8+16) n = accept + reject;

    if (changeMode & 32) dn = double(n) / double(n_walkers);
    else dn = double(n);

    gamma = 1.0 / dn;
    if (gamma > 1.0) gamma = 1.0;

    initKernel(kernelType, dosKernel);
    dosKernel.scale(gamma);
  }

// 1. write configuration
// 2. check histogram flatness
// 3. perform global update
  if (cycleCount >= updateCycle)
  {
    cycleCount = 0;
    // syncronizeGraphs(dos, histo);
    if (dos.getN() != histo.getN())
    {
      std::cout << "Histogramm size dosn't match DOS! Clearing histogramm!\n";
      histo.setRangeAndClear(dos.getMinX(), dos.getMaxX(), dos.getN());
    }

    writeState("WL1d.state");

    if (!histogramUpdateMode)        // Normal Wang-Landau
    {
      // calculate minimum nonzero histogram
      // we look only at the histogram inside the energy interval that was actually sampled if we use kernel updates
      double hMin, hMax, hMean;

      if (kernelType == None)
      {
        hMean = histo.getMeanY();
        histo.getMinMaxY(hMin,hMax);
      }
      else 
      {
        hMean = histo.getMeanYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                         dos.getMaxX() - dosKernel.getWidth());
        histo.getMinMaxYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                   dos.getMaxX() - dosKernel.getWidth(),
                                   hMin, hMax);
      }

      // double currentFlatnessCriterion = double(hMin-hMax) / hMean;
      double currentFlatnessCriterion = double(hMin) / hMean;
      
      std::cout << "# acceptence ratio = " << double(accept) / double(accept+reject) << "\n";
      std::cout << "# current flatness = " << currentFlatnessCriterion << (changeMode & 4 ? " *":"") << "\n";
      std::cout << "# current histogram minimum = " << hMin << (changeMode & 2 ? " *":"") << "\n";
      std::cout << "# average accepted steps/bin since last gamma change = "
                << double(acceptSinceLastChange) / double(histo.getN())
                << (changeMode & 1 ? " *":"") << "\n";

      if (changeMode != 0 && changeMode < 8)
      {
        // perform global update
        if (globalUpdate.frequency > 0 && (globalUpdate.frequency*histo.getN()) < acceptSinceLastChange) 
        {
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDOS(fn);
          globalUpdate.changes++;
          std::cout<<"# global update "<<globalUpdate.changes<<std::endl;
          performGlobalUpdate(dos,globalUpdate.kappa,globalUpdate.lambda,globalUpdate.omega);
          histo.clear();
          acceptSinceLastChange=0;
        }
        else if (((acceptSinceLastChange>=histo.getN()*updatesPerBin || !(changeMode &1)) &&
                  (hMin >= hMinimum || !(changeMode & 2)) &&
                  (currentFlatnessCriterion>flatnessCriterion || !(changeMode & 4)) ))
        {
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " is finished.\n";
          modificationFactorChanges++; // counter[2]

          // write dos;
          // we are using JSON as our file format
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDOS(fn);
          std::cout << __FILE__ << ":" << __LINE__ << " fn=" << fn << " dos_out_name=" << dos_out_name << std::endl;
          // dos_out_name local?
          writeDOS(dos_out_name.data());

          // clear the Histogram
          histo.clear();
          acceptSinceLastChange = 0;
 
          // change gamma
          dosKernel.scale(0.5);
          gamma = 0.5 * gamma;
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
          if(statesFile != NULL)
    {
            // char fn[256];
            snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
            sw.newFile(fn);
            sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
          }
        }
      }
    } 
    else            // histogram update mode != 0
    {
      if (acceptSinceLastChange >= histo.getN() * updatesPerBin)
      {
        acceptSinceLastChange = 0;
        for (int ix=0; ix<dos.getN(); ix++)
          dos[ix] += gamma * histo[ix];
        histo.clear();
        gamma = 0.5 * gamma;

        std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
        if (statesFile != NULL)
        {
          char fn[256];
          snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
          sw.newFile(fn);
          sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
        }
      }
    }
  }

  if (accepted) 
  {
    // Suffian: Moved to routine updateLog(...)
    /* sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
    lastAccepted[instance] = lastChange[instance];
    lastAcceptedEnergy[instance] = position[instance];
    //lastAcceptedEnergy[instance] = energy;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]]; */
    accept++;
    acceptSinceLastChange++;
  }
  else reject++;

  if (gamma < gammaFinal) return true;
  else return false;

}

template<class RNG>
void WL1dEvecGenerator<RNG>::updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
    bool isSpinSim, bool isOccupancySim) {

  // use static variables to keep track of the last move type
  static bool firstrun = true;
  static std::vector<MoveChoice_t> lastAcceptedMove;
  static FILE* filep;

  // initialize static variables
  if( firstrun ) {
    lastAcceptedMove.resize(n_walkers);

    if( isOccupancySim ) {
      filep = fopen("occ_vs_energies","w");
      fprintf(filep,"# site occupancies, energy\n");
    }
    firstrun = false;
  }

  // Suffian: my preference is to print all configurations and energy
  if( isOccupancySim ) {
    for(int i = 0; i < n_spins; i++)
      fprintf(filep, "%1d", occ[i]);
    fprintf(filep, "   %20.15f", energy);
    fprintf(filep, "   %5d", instance);
    if( accepted ) fprintf(filep, "   accepted");
    fprintf(filep, "\n");
    fflush(filep);
  }

  // return unless state changed
  if( !accepted ) return;

  // print output corresponding to previous state
  if( lastAcceptedMove[instance] == SpinMove )
    sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
  else
    sw.writeChangeOcc(instance, numRetentions[instance], 
        lastAcceptedSwap[instance].first,
        lastAcceptedSwap[instance].second,
        lastAcceptedOcc[instance].first, 
        lastAcceptedOcc[instance].second, 
        lastAcceptedEnergy[instance]); 

  // remember the new configuration
  lastAcceptedEnergy[instance] = position[instance];
  if( MoveChoice == SpinMove ) {
    lastAccepted[instance] = lastChange[instance];
    //lastAcceptedEnergy[instance] = energy;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]]; 
  }
  else {
    lastAcceptedSwap[instance].first  = lastSwapOcc[instance].first;
    lastAcceptedSwap[instance].second = lastSwapOcc[instance].second;
    lastAcceptedOcc[instance].first  = occ[ lastSwapOcc[instance].first ];
    lastAcceptedOcc[instance].second = occ[ lastSwapOcc[instance].second ];
  }

  lastAcceptedMove[instance] = MoveChoice;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generatePotentialShift(int instance, double *potentialShifts, bool accepted)
{

  if (accepted)
  {
    lastChangePotentialShift[instance] = int(rnd(rng)*n_spins);
    oldPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
    potentialShifts[lastChangePotentialShift[instance]] = randomPotentialShift();
  } else {
    potentialShifts[lastChangePotentialShift[instance]] = oldPotentialShift[instance];
    lastChangePotentialShift[instance] = int(rnd(rng)*n_spins);
    oldPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
    potentialShifts[lastChangePotentialShift[instance]] = randomPotentialShift();
  }

}

template<class RNG>
void WL1dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, bool accepted)
{
  if (accepted)
  {
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
    numRetentions[instance] = 0;
  }
  else 
  {
    evecs[  3*lastChange[instance]] = oldSpin[  3*instance];
    evecs[1+3*lastChange[instance]] = oldSpin[1+3*instance];
    evecs[2+3*lastChange[instance]] = oldSpin[2+3*instance];
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
  }

}

// additional routines for Wang-Landau for alloying
// swapping atoms preserves concentration

template<class RNG>
MoveChoice_t WL1dEvecGenerator<RNG>::selectMoveType(bool isSpinSim, bool isOccSim) {

  // choose between spin and occupancy equally
  // is this a good choice??
  if( isSpinSim && isOccSim ) {
    int c = int(rnd(rng)*2);
    if( c == 0 ) 
      return SpinMove;
    else
      return OccupancyMove;
  }

  if( isSpinSim )
    return SpinMove;
  if( isOccSim )
    return OccupancyMove;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::setAlloyClasses(const AlloyMixingDesc& alloyDesc, int *siteclass) {
  this->alloyDesc = alloyDesc;

  siteAlloyClass.resize(n_spins);
  for(int i = 0; i < n_spins; i++)
    siteAlloyClass[i] = siteclass[i];

  // remember number of sites per alloy class 
  int nclasses = alloyDesc.size();
  numSites_per_AlloyClass.resize(nclasses);
  int *nsites = numSites_per_AlloyClass.data();
  for(int i = 0; i < nclasses; i++)
    nsites[i] = 0;
  for(int i = 0; i < n_spins; i++)
    nsites[ siteAlloyClass[i] ]++;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generateOccupancies(int instance, int *occ, bool acceptedOcc) {

  int temp, atom_i, atom_j;

  // if previous move exists and not accepted 
  // revert occupancies
  if( !acceptedOcc && lastSwapOcc[instance].first >= 0 ) {
    atom_i = lastSwapOcc[instance].first;
    atom_j = lastSwapOcc[instance].second;

    temp = occ[atom_i];
    occ[atom_i] = occ[atom_j];
    occ[atom_j] = temp;

    /*
    printf("Master: Reverted occupancies:");
    for(int i = 0; i < n_spins; i++)
      printf("%d ",occ[i]);
    printf("\n"); */
  }
  else
    numRetentions[instance] = 0;

  // pick an alloy class
  int ac = int(rnd(rng)*alloyDesc.size());

  // pick two distinct atoms within mixing class
  // here we assume there are few alloy classes
  do{ atom_i = int(rnd(rng)*n_spins); } 
    while( siteAlloyClass[atom_i] != ac );
  do{ atom_j = int(rnd(rng)*n_spins); } 
    while( siteAlloyClass[atom_j] != ac || occ[atom_j] == occ[atom_i] );

  // swap atoms
  temp = occ[atom_i];
  occ[atom_i] = occ[atom_j];
  occ[atom_j] = temp;

  // remember which atoms were swapped
  lastSwapOcc[instance].first  = atom_i;
  lastSwapOcc[instance].second = atom_j;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generateUnsampledOcc(int inst, int *occ) {

  // determine number atoms of each component type given concentration
  std::vector< std::vector<int> > atomsLeft;
  atomsLeft.resize(alloyDesc.size());
  for(int i = 0; i < atomsLeft.size(); i++) {
    atomsLeft[i].resize(alloyDesc[i].size());
    for(int j = 0; j < atomsLeft[i].size(); j++) {

      // how many sites for this alloy class and component?
      double nfrac = alloyDesc[i][j].conc * numSites_per_AlloyClass[i];
      double error = fmod(nfrac, 1.0);

      // ensure concentrations commensurate with number of sites
      const double tol = 1.e-06;
      if( tol < error && error < 1.0-tol ) {
        printf("error: alloy concentrations are incommensurate with supercell\n");
        printf("alloy class = %d, component = %d, concentration = %f\n", 
            i+1, j+1, alloyDesc[i][j].conc);
        printf("desired sites = %f\n", nfrac);
        exit(1);
      }

      atomsLeft[i][j] = int(nfrac + 0.5);
    }
  }

  // distribute atoms across sites
  for(int i = 0; i < n_spins; i++) {

     // pick an atom from the right class
     int type, ac = siteAlloyClass[i];
     do { type = int( rnd(rng) * alloyDesc[ac].size() ); }
       while( atomsLeft[ac][type] == 0 );

     // assign to site
     occ[i] = type; atomsLeft[ac][type]--;
  }
}
 
template<class RNG>
void WL1dEvecGenerator<RNG>::initializeOccupancies(int inst, int *occ) {

  // Warning: do nothing if walker occupancies read from file
  if( inst < n_initialized_from_file_occ )
    return;

  generateUnsampledOcc(inst, occ);
}


#endif

//#include "WangLandau2d.h"

// -*- mode: c++ -*-
#ifndef LSMS_WANG_LANDAU_2d_H
#define LSMS_WANG_LANDAU_2d_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
// we use BOOST for the random number generator
// #include <boost/random.hpp>
#include <random>
#include "../../mjson/json.h"
#include "EvecGenerator.h"
#include "Graph2d.hpp"

class StatesWriter2d
{
public:
  StatesWriter2d(const char *filename=NULL)
  {
    if(filename==NULL) writeFlag=false;
    else {writeFlag=true; of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);}
  }
  ~StatesWriter2d() {if(writeFlag) of.close();}
  void writeHeader(double lnF, int numWalk, int numSpin, double **spins)
  {
// browngrg 20/14/2015
if(spins==NULL) return;
    if(writeFlag)
      {
	of<<lnF<<" "<<numWalk<<" "<<numSpin<<std::endl;
	for(int i=0; i<numWalk; i++)
	  {
	    of<<i;
	    for(int j=0; j<3*numSpin; j++) of<<" "<<spins[i][j];
	    of<<std::endl;
	  }
      }
  }
  void writeChange(int iWalk, int numRet, int ispin, double *ev, double E)
  {
    if(writeFlag)
      of<<iWalk<<" "<<numRet<<" "<<ispin<<" "<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<" "<<E<<std::endl;
  }
  void newFile(const char *filename=NULL)
  {
    if(writeFlag) of.close();
    writeFlag=false;
    if(filename!=NULL){
      writeFlag=true;
      of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);
    }
  }
private:
  bool writeFlag;
  std::ofstream of;
};

// template<class RNG = boost::mt19937>
template<class RNG = std::mt19937>
class WL2dEvecGenerator : public EvecGenerator
{
 public:
  WL2dEvecGenerator(int num_spins, int num_instances, double **ev_p,
                    const char *init_file_name = NULL,
                    const char *out_file_name = NULL,
                    const char *states_filename = NULL);

  bool determineAcceptance(int instance, double energy);
  bool determineAcceptance(int instance, double energy, double magnetization);

  bool updateHistogram(int instance, double *evecs, bool accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted);
/*
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization)
  {
    bool h;
    return updateHistogram(instance, evecs, energy, magnetization, &h);
  }
  bool updateHistogram(int instance, double *evecs)
  {
    std::cerr << "Need energy for WL2dEvecGenerator\n"; 
    exit(1);
  }
*/

  void generateEvec(int instance, double *evecs, bool accepted);
  //void generateEvec(int instance, double *evecs, double energy);

  void initializeEvec(int instance, double *evecs);

  void generateUnsampledEvec(int instance, double *evecs, double energy)
  {
    initializeEvec(instance, evecs);
    //return false;
  }

  void startSampling(void)
  { sw.writeHeader(gamma, n_walkers, n_spins, evecs_pointer); }

  void writeState(const char *name);

 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
  // boost::uniform_real<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);

  /*
  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  */

  double hMinimum;
  Graph2d<double,double> dos, histo;
  Kernel2d<double,double> dosKernel, histoKernel;
  KernelType kernelType;

  unsigned long accept, reject;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  // instance specific:
  std::vector<double> positionX, positionY;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<int> lastAccepted;
  std::vector<double> lastAcceptedEnergy;
  std::vector<double> lastAcceptedMagnetization;
  std::vector<double> lastAcceptedEvec;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]

  std::vector<int> numRetentions;

  char *statesFile;
  StatesWriter2d sw;


  void inline random_evec_1(double ev[3])
  {
    double x,y,z;
    do {
      x = rnd(rng);
      y = rnd(rng);
    } while(x*x+y*y>1);
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r;
    y *= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z;
    r=1.0/sqrt(x*x+y*y+z*z);
    ev[0]=x*r; ev[1]=y*r; ev[2]=z*r;
  }

  void inline random_evec(double ev[3])
  {
    double x, y, z;
    do {
      x = rnd(rng); y = rnd(rng);
    } while(x*x+y*y>1); 
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r; y*= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z; 
    // Project out the parallel component;
    r = x*ev[0] + y*ev[1] + z*ev[2];
    x -= r*ev[0]; y -= r*ev[1]; z -= r*ev[2];
    r = x*x + y*y + z*z;
    double t = 1-0.3*rnd(rng);
    ev[0] *= t; ev[1] *= t; ev[2] *= t;
    r = sqrt((1-t*t)/r);
    ev[0] += x*r; ev[1] += y*r; ev[2] += z*r;
    r=1.0/sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    ev[0]*=r; ev[1]*=r; ev[2]*=r;
    
    /*  
    ev[2]=1.0-2.0*rnd(rng);
    // ev[2]=rnd11(rng);
    double phi=2.0*M_PI*rnd(rng);
    // double phi=rnd0pi(rng);
    double cos_theta=sqrt(1-ev[2]*ev[2]);
    ev[0]=cos_theta*cos(phi);
    ev[1]=cos_theta*sin(phi);
    */  
  }

};

template<class RNG>
WL2dEvecGenerator<RNG>::WL2dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
					  const char *init_file_name, const char *out_file_name, const char *states_filename)
  : sw(states_filename)
{
  long nX=-1;
  double intervalEnergy=0.1;
  double intervalMagnetization=0.1;
  double kernelWidthEnergy=1.0;
  double kernelWidthMagnetization=1.0
;
  double xMin=-std::numeric_limits<double>::max();
  double xMax=1.0;
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  evecs_pointer = ev_p;
  positionX.resize(n_walkers);
  positionY.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedMagnetization.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);
  for(int i=0; i<n_walkers; i++)
  {
    lastAccepted[i]=-2;
  }
  for(int i=0; i<3*n_walkers; i++)  lastAcceptedEvec[i]=0.0;

  oldSpin.resize(3*n_walkers);

  statesFile=NULL;
  if(states_filename!=NULL)
    {
      statesFile=(char *)malloc(sizeof(char)*(1+strlen(states_filename)));
      strcpy(statesFile,states_filename);
    }

  numRetentions.resize(n_walkers);
  for(int i=0; i<n_walkers; i++) numRetentions[i]=0;

  /*
  nX = -1;
  xMin = -HUGE; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  */
  dos_out_name="dos2d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  cycleCount=0;
  hMinimum= 1;     //   10
  accept=reject=0;
  gammaFinal=1.e-6;
  flipPerUpdate=4; //  100
  updateCycle= 1000; // 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.75;

  kernelType=Epanechnikov;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(intervalEnergy,intervalMagnetization);
  histo.setDeltaAndClear(intervalEnergy,intervalMagnetization);
  if(init_file_name!=NULL && init_file_name[0]!=0)
  {
    std::string label, value;

    std::cout<<"Reading Wang-Landau configuration from: "<<init_file_name<<std::endl;

    dos_out_name=std::string(init_file_name)+".out";
    if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
    std::ifstream inp(init_file_name);
    std::ostringstream buff;

    std::string line;
    while(std::getline(inp,line)) 
      buff << line << std::endl;

    inp.close();

    std::string fileString = buff.str();
    const char* fileChars  = fileString.c_str();
    json_t *json_root=NULL;

    json_parse_document(&json_root, (char *)fileChars);

    if(json_root == NULL || json_root->type != JSON_OBJECT)
    {
      std::ostringstream message;
      std::cerr << "In WL2dEvecGenerator(" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
  
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;
      if(label=="xMin") xMin=atof(it->child->text);
      else if(label=="xMax") xMax=atof(it->child->text);
      else if(label=="intervalEnergy") intervalEnergy=atof(it->child->text);
      else if(label=="intervalMagnetization") intervalMagnetization=atof(it->child->text);
      else if(label=="kernelWidthEnergy") kernelWidthEnergy=atof(it->child->text);
      else if(label=="kernelWidthMagnetization") kernelWidthMagnetization=atof(it->child->text);
      else if(label=="kernelType")
      {
        std::string strValue(it->child->text);
        kernelType=getKernelType(strValue);
      }
      else if(label=="gamma") gamma=atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal=atof(it->child->text);
      else if(label=="nX")
      {
        nX=atoi(it->child->text);
	dos.setXRangeAndClear(xMin,xMax,nX);
	histo.setXRangeAndClear(xMin,xMax,nX);
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="cycleCount") cycleCount=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
      else if(label=="seed") rng.seed(atoi(it->child->text));
      else if(label=="accept") accept=atol(it->child->text);
      else if(label=="reject") reject=atol(it->child->text);
//*
      else if(label=="rngState")
      {
        std::string strValue(it->child->text);
        std::stringstream strStream(strValue, std::stringstream::in);
        strStream>>rng;
      }
//*/
      else if(label=="dos")
      {
        json_t *a = it->child;
        int ix=0;
        int iy;
        long Ny;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          Ny=atoi(i->text); i=i->next;
          double minY=atof(i->text); i=i->next;
          dos.setYminAndClear(ix,minY,Ny);
          iy=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
            dos[ix][iy++]=atof(j->text);
          ix++;
        }
        if(ix!=dos.getNx()) {std::cout<<"ERROR #(dos) "<<ix<<" != nX "<<dos.getNx()<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int ix=0;
        int iy;
        long Ny;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          Ny=atoi(i->text); i=i->next;
          double minY=atof(i->text); i=i->next;
          histo.setYminAndClear(ix,minY,Ny);
          iy=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
            histo[ix][iy++]=atof(j->text);
          ix++;
        }
        if(ix!=histo.getNx()) {std::cout<<"ERROR #(histo) "<<ix<<" != nX "<<histo.getNx()<<std::endl; exit(1);}
      }
      else if(label=="evecs")
      {
// browngrg 10/14/2015
if( evecs_pointer==NULL )
{
   std::cout << __FILE__ << ":" << __LINE__ << " cannot read evecs" << std::endl;
   continue;
}
        json_t *a = it->child;
        n_initialized_from_file=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              evecs_pointer[n_initialized_from_file][k++]=atof(j->text);
            }
            n_initialized_from_file++;
            if(k!=3*n_spins) {std::cout<<"ERROR #(evecs) "<<k<<" != 3*n_spins "<<3*n_spins<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="oldSpin")
      {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              oldSpin[n_initialized*3+k++]=atof(j->text);
            }
            n_initialized++;
            if(k!=3) {std::cout<<"ERROR #(oldSpin) "<<k<<" != 3"<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="lastChange")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers) lastChange[j++]=atoi(i->text);
          n_initialized++;
        }
      }
      else if(label=="position")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            positionX[j]=atof(i->text);
            i=i->next;
            positionY[j]=atof(i->text);
            n_initialized++;
          }
        }
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else std::cout << __FILE__ << ":" << __LINE__ << " WARNING: unknown label: " << label << std::endl;
    }

    // set xMax and xMin or interval depending on nX:
    if(nX>0)
    {
      intervalEnergy=(xMax-xMin)/double(nX);
      dos.setXRange(xMin,xMax); histo.setXRange(xMin,xMax);
    } else {
      dos.setDeltaAndClear(intervalEnergy,intervalMagnetization);
      histo.setDeltaAndClear(intervalEnergy,intervalMagnetization);
    }

    json_free_value(&json_root);
  }
  dosKernel.setWidthAndClear(intervalEnergy, intervalMagnetization,
                             kernelWidthEnergy, kernelWidthMagnetization);
  histoKernel.setWidthAndClear(intervalEnergy, intervalMagnetization,
                             kernelWidthEnergy, kernelWidthMagnetization);
  initKernel(kernelType,dosKernel);
  dosKernel.scale(gamma/dosKernel(0.0,0.0));
  initKernel(kernelType,histoKernel);
  histoKernel.scale(1.0/histoKernel(0.0,0.0));
//  histoKernel.scale(interval);

  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout<<"Wang-Landau output will be written to: "<<dos_out_name<<std::endl;
}

template<class RNG>
void WL2dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  if(inst>=n_initialized_from_file)
    {
      if(inst==0)
	for(size_t j=0; j<3*n_spins; j+=3)
	  {evecs[j]=0.0; evecs[j+1]=0.0; evecs[j+2]=1.0;}
      else if(inst==1)
        for(size_t j=0; j<3*n_spins; j+=3)
          {evecs[j+2]= (j/3)%2 ? 1.0 : -1.0 ; evecs[j+1]=0.0; evecs[j]=0.0;}
      else
	for(size_t j=0; j<3*n_spins; j+=3)
	  random_evec_1(&evecs[j]);
    }

  out[inst]=false;
}

template<class RNG>
void WL2dEvecGenerator<RNG>::writeState(const char* name)
{
// browngrg 10/14/2015
return;
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(8);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    ofile<<"\"nX\" : " << dos.getNx() << ",\n";
    ofile<<"\"intervalMagnetization\" : " << dos.getDeltaY() << ",\n";
    ofile<<"\"kernelWidthEnergy\" : " << dosKernel.getWidthX() << ",\n";
    ofile<<"\"kernelWidthMagnetization\" : " << dosKernel.getWidthY() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma << ",\n";
    ofile<<"\"accept\" : " << accept <<",\n";
    ofile<<"\"reject\" : " << reject <<",\n";
    ofile<<"\"gammaFinal\" : " << gammaFinal << ",\n";
    ofile<<"\"flatnessCriterion\" : "<< flatnessCriterion <<",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << ",\n";
    ofile<<"\"updateCycle\" : " << updateCycle << ",\n";
    ofile<<"\"dos\" : ["<<std::endl;
    for(int ix=0; ix<dos.getNx(); ix++)
    {
      ofile<<dos.getNy(ix)<<", "<<dos.getMinY(ix)<<", ["<<std::endl;
      for(int iy=0; iy<dos.getNy(ix); iy++)
        ofile<<dos[ix][iy]<<((iy==dos.getNy(ix)-1)?"\n":",\n");
      ofile<<((ix==dos.getNx()-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int ix=0; ix<dos.getNx(); ix++) 
    { 
      ofile<<histo.getNy(ix)<<", "<<dos.getMinY(ix)<<", ["<<std::endl;
      for(int iy=0; iy<histo.getNy(ix); iy++)
        ofile<<histo[ix][iy]<<((iy==histo.getNy(ix)-1)?"\n":",\n");
      ofile<<((ix==histo.getNx()-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"rngState\" : \""<<rng<<"\",\n";
    ofile<<"\"evecs\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[\n";
      for(int j=0; j<3*n_spins; j+=3)
        ofile<<evecs_pointer[i][j]<<", "<<evecs_pointer[i][j+1]
             <<", "<<evecs_pointer[i][j+2]<<((j==3*n_spins-3)?"\n":",\n");
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"oldSpin\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[ "<<oldSpin[3*i]<<", "<<oldSpin[3*i+1]
             <<", "<<oldSpin[3*i+2];
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"lastChange\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< lastChange[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"position\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< positionX[i]<<", "<<positionY[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"stepsSinceLastHistogramUpdate\" : " << stepsSinceLastHistogramUpdate << ",\n";
    ofile<<"\"numberOfUpdatesSinceLastBoost\" : " << numberOfUpdatesSinceLastBoost << ",\n";
    ofile<<"\"modificationFactorChanges\" : " << modificationFactorChanges << ",\n";
    ofile<<"\"cycleCount\" : " << cycleCount << "\n";
    ofile<<"}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
} 


template<class RNG>
bool WL2dEvecGenerator<RNG>::determineAcceptance(int instance, double energy)
{
  return determineAcceptance(instance, energy, 0.0);
}


template<class RNG>
bool WL2dEvecGenerator<RNG>::determineAcceptance(int instance, double energy, double magnetization)
{

  bool accept_step;
 
  // energy between xMin and xMax ?
  long gridX = dos.idxX(energy);
  // we have to postpone geting the y position until we have checked the validity of gridX
  //  int gridY = dos.idxY(gridX, magnetization);

  stepsSinceLastHistogramUpdate++; // counter[0]++

  // record initial position for step out
  if (lastAccepted[instance] == -2)
  {
    positionX[instance] = lastAcceptedEnergy[instance] = energy;
    positionY[instance] = lastAcceptedMagnetization[instance] = magnetization;
    lastAccepted[instance] = -1;
  }

  // out of energy range
  if (gridX < 0 || gridX > dos.getNx()-1)
  {
    dos.extendToX(energy - dosKernel.getWidthX()); histo.extendToX(energy - histoKernel.getWidthX());
    dos.extendToX(energy + dosKernel.getWidthX()); histo.extendToX(energy + histoKernel.getWidthX());
    gridX = dos.idxX(energy);  // gridY=dos.idxY(gridX, magnetization);
  }
  long gridY = dos.idxY(gridX, magnetization);
  // out of magnetization range
  if (gridY < 0 || gridY > dos.getNy(gridX)-1)
  {
    dos.extendToY(gridX, magnetization - dosKernel.getWidthY());
    dos.extendToY(gridX, magnetization + dosKernel.getWidthY());
    histo.extendToY(gridX, magnetization - histoKernel.getWidthY());
    histo.extendToY(gridX, magnetization + histoKernel.getWidthY());
    gridY = dos.idxY(gridX, magnetization);
  }
 
  numRetentions[instance]++;
  out[instance] = false;
  // ref1X[instance] = gridX; ref1Y[instance] = gridY;
  long refX = dos.idxX(positionX[instance]);
  long refY = dos.idxY(refX,positionY[instance]);

  if (lastAccepted[instance] < 0)
  {
    accept_step = true;
    positionX[instance] = energy;
    positionY[instance] = magnetization;
  }
  else
  {
    if (dos[gridX][gridY] <= dos[refX][refY])
    {
      accept_step = true;
      positionX[instance] = energy;
      positionY[instance] = magnetization;
    } 
    else 
    {
      if(rnd(rng) < exp(dos[refX][refY] - dos[gridX][gridY]))
      {
        accept_step = true;
        positionX[instance] = energy;
	positionY[instance] = magnetization;
      } 
      else 
      {
        accept_step = false;
      }
    } 

  }

  if (verbosity > 0)
    std::cout << "WangLandau 2d EvecGenerator step "
              << modificationFactorChanges << ":" << numberOfUpdatesSinceLastBoost << ":"
              << stepsSinceLastHistogramUpdate << " nX=" << dos.getNx()
              << " [" << dos.getMinX() << ", " << dos.getMaxX() << "] "
              << (accept_step ? "accepted" : "rejected")
              << " Energy = " << energy << ", Magnetization = " << magnetization
	      << ", Instance " << instance << std::endl;

  return accept_step;

}


template<class RNG>
bool WL2dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, bool accepted)
{
  std::cerr << "Need energy and magnetization for updateHistogram in WL2dEvecGenerator.\n";
  std::cerr << "Either wl_lsms or updateHistogram needs to be amended!\n"; 
  exit(1);
}


template<class RNG>
bool WL2dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted)
{

// Update histogram and DOS
  if (stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      //addKernel(dos, dosKernel, positionX[instance], positionY[instance]);
      //addKernel(histo, histoKernel, positionX[instance], positionY[instance]);
      addKernel(dos, dosKernel, energy, magnetization);
      addKernel(histo, histoKernel, energy, magnetization);
    } 
    else
    {
      std::cerr << "ATTENTION: We should never reach this place in WL2dEvecGenerator!\n";
      exit(1);
    }
  }

// 1. write configuration
// 2. check histogram flatness
  if (cycleCount >= updateCycle)
  {
    cycleCount = 0;
    // syncronizeGraphs(dos,histo);
    writeState("WL2d.state");
    // calculate minimum nonzero histogram
    /*
    int hMin = histo[0];
    int hMax = histo[0];
    double hMean = 0.0;
    for(int i=1; i<nX; i++)
    {
      hMean += double(histo[i]);
      if(histo[i]>=hMax) hMax=histo[i];
      if(histo[i]<hMin) hMin=histo[i];
    }
    hMean /= double(nX);
    double currentFlatnessCriterion = double(hMax-hMin)/hMean;
    
    std::cout <<"# accepence ratio = "<<double(accept)/double(accept+reject)<<"\n";
    std::cout <<"# current flatness = "<<currentFlatnessCriterion<<"\n";
    */

    //if(hMin >= hMinimum && currentFlatnessCriterion<=flatnessCriterion)
    // we look only at the histogram inside the energy interval that was actually sampled
    double hMin = histo.getMinValWithBorders(dosKernel.getWidthX(), dosKernel.getWidthY());
    std::cout << "# acceptence ratio = " << double(accept) / double(accept+reject) << "\n";
    std::cout << "# current histogram minimum = " << hMin << "\n";
    if (hMin >= hMinimum)
    {
      std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " is finished.\n";
      modificationFactorChanges++; // counter[2]

      // write dos;
      // we are using JSON as our file format
      writeState(dos_out_name.data());
 
      // else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";

      // clear the Histogram
      histo.clear();

      // change gamma
      dosKernel.scale(0.5);
      gamma = 0.5 * gamma;
      std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
      if (statesFile != NULL)
      {
        char fn[256];
        snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
        sw.newFile(fn);
        sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
      }
    }
  }


  if (accepted)
  {
    sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
    lastAccepted[instance] = lastChange[instance];
    lastAcceptedEnergy[instance] = energy;
    lastAcceptedMagnetization[instance] = magnetization;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]];
    accept++;
  }
  else reject++;

  if (gamma<gammaFinal) return true;
  else return false;

}

template<class RNG>
void WL2dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, bool accepted)
{
  if (accepted)
  {
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
    numRetentions[instance] = 0;
  }
  else
  {
    evecs[  3*lastChange[instance]] = oldSpin[  3*instance];
    evecs[1+3*lastChange[instance]] = oldSpin[1+3*instance];
    evecs[2+3*lastChange[instance]] = oldSpin[2+3*instance];
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
  }

}


#endif





















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
#include "ExhaustiveIsing.h"

#include "Heisenberg_Hamiltonian.hpp"


























// Occupation simulations in magnetic alloys demand evec
struct LSMS_Config
{
public:
   double energy;
   double band_energy;
   double wl2nd;
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
  int nwalk;                 // number of parallel LSMS instances
  int natom;                // number of atoms in a lsms instance
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

   std::vector<double> send_buffer;
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
     // size has to be align + natom*nwalk
     align = 1;
     nwalk = 1;
     natom = -1;
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
   void do_without_master(); 

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
    if (!strcmp("-nwalk", argv[i])) nwalk = atoi(argv[++i]);
    if (!strcmp("-natom", argv[i])) natom = atoi(argv[++i]);
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
  options.add_option("nwalk",      "Number of lsms calculators",     ' ', &(nwalk) );
  options.add_option("natom",     "Number of sites in calculation", ' ', &(natom) );
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
    std::cout << " Size of LSMS instances = " << natom << " atoms\n";
    std::cout << " Number of LSMS instances = " << nwalk << std::endl;
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

    if (step_out_flag)
    {
      step_out_file.open(step_out_name);
      step_out_file << "#";
      //for(int i=0; i<argc; i++) step_out_file << " " << argv[i];
      step_out_file << std::endl << natom << std::endl;
    }
}


void LSMS_Hamiltonian::init()
{
    // build the communicators
    int color = MPI_UNDEFINED;
    int s = align;
    comm_size = (size-align) / nwalk;
    for (int i=0; i<nwalk; i++)
    {
      if ((world_rank >= s) && (world_rank < s+comm_size))
      {
        my_group = i;
        color = i;
      }
      lsms_rank0[i] = s;
      s += comm_size;
    }
    if (world_rank == 0) color = nwalk;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &local_comm);
    snprintf(prefix, 38, "Group %4d: ", my_group);
}


int LSMS_Hamiltonian::get_op(LSMS_Config& sigma, std::vector<double>& recv_buffer)
{
   int op;
   if (rank == 0)
   {
     // recieve command in op code
     // data represents either site spins or site occupancies 
     MPI_Recv(&(recv_buffer[0]), 4*natom, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
     // use op code to distinguish whether occupancy or spin variables recieved
     op = status.MPI_TAG;
     // dprintf("Walker %d: Recieved op code %d.\n", my_group, op);
     if( op == 51 || op == 61 ) {
       for(int i = 0; i < 4*natom; i++)
         sigma.evec[i] = recv_buffer[i];
     }
     else if( op == 52 || op == 62 ) {
       for(int i = 0; i < natom; i++)
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
        MPI_Send( &(r_values[0]), R_VALUE_OFFSET + 3*natom, MPI_DOUBLE, 0, 1005, MPI_COMM_WORLD);
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
        for(int i = 0; i < natom; i++)
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
   recv_buffer.resize(4*natom);
   sigma.occ.resize(natom);
   sigma.evec.resize(4*natom);
   r_values.resize( R_VALUE_OFFSET + 3*(natom+1) );
   i_values.resize(10);

   MPI_Comm_rank(local_comm, &rank);
   snprintf(prefix, 38, "%d_", my_group);
   // to use the ramdisk on jaguarpf:
   // snprintf(prefix, 38, "/tmp/ompi/%d_", my_group);
   // This is from lsmsClass.hpp
   LSMS lsms_calc(local_comm, i_lsms_name, prefix, my_group);
   snprintf(prefix, 38, "Group %4d: ", my_group);

   std::string wl_inf;
   if (gWL_in_name) 
   {
      wl_inf = gWL_in_name;
   }
   std::string wl_outf;
   if (gWL_out_name)
   {
      char buffer[128];
      snprintf(buffer,128,"%s_%d",gWL_out_name,my_group);
      wl_outf = buffer;
   }
   std::string wl_stepf;
   if (step_out_flag && (evec_generation_mode == WangLandau_1d))
   {
     // step_out_flag = false;
     char buffer[128];
     snprintf(buffer,127,"wl1d_%s_%d",step_out_name,my_group);
     wl_stepf = buffer;
   }
   std::cout << __FILE__ << ":" << __LINE__ << " wang-landau slave output=\"" << wl_outf << "\"" << std::endl;

   PotentialShifter potentialShifter;
   EvecGenerator *generator;

   // vSpinShiftFlag is hard-wired inside Potential Shifter
   if(potentialShifter.vSpinShiftFlag)
   {
      sigma.vSpinShifts.resize(natom);
      for (int j = 0; j < natom; j++) 
        sigma.vSpinShifts[j] = 0.0;
      sigma.evecsAndSpinShifts.resize(4*natom);
   }

   // Factory generator here
   switch (evec_generation_mode)
   {
   case Random : 
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      generator = new RandomEvecGenerator(natom);
      break;
   case Constant:
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      generator = new ConstantEvecGenerator(natom, ev0, nwalk);
      break;
   case WangLandau_1d : 
      //generator = new WL1dEvecGenerator<std::mt19937>(natom, nwalk, evecs, potentialShifter, wl_inf, wl_outf.c_str(), wl_stepf, occs);
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      generator = new WL1dEvecGenerator<std::mt19937>(natom, nwalk, NULL, potentialShifter, wl_inf.c_str(), wl_outf.c_str(), wl_stepf.c_str(), NULL);
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      break;
   case ExhaustiveIsing : 
      //generator = new ExhaustiveIsing1dEvecGenerator(natom, nwalk, evecs, wl_inf, wl_outf);
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      generator = new ExhaustiveIsing1dEvecGenerator(natom, nwalk, NULL, wl_inf.c_str(), wl_outf.c_str());
      break;
   case WangLandau_2d : 
      //generator = new WL2dEvecGenerator<std::mt19937>(natom, nwalk, evecs, wl_inf, wl_outf, wl_stepf);
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      generator = new WL2dEvecGenerator<std::mt19937>(natom, nwalk, NULL, wl_inf.c_str(), wl_outf.c_str(), wl_stepf.c_str());
      break;
   default : std::cerr << "The code should never arrive here: UNKNOWN EVEC GENERATION MODE\n";
      exit(1);
   }

   // Generate the initial spin configuration 
   if (potentialShifter.vSpinShiftFlag)
      generator -> initializeEvecAndPotentialShift(0, &(sigma.evec[0]), &(sigma.vSpinShifts[0]));
   else 
      generator -> initializeEvec(0, &(sigma.evec[0]) );
   if( isOccupancySim ) 
   {
      generator->initializeOccupancies(0,&(sigma.occ[0]));
   }

   MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);

   generator -> startSampling(isSpinSim, isOccupancySim);

   // wait for commands from master
   bool finished = false;
   while (!finished)
   {
      int iwalk = 0;
      bool accepted = true;
      int op = get_op(sigma,recv_buffer);
      int calc = op/10;
      int imove = op - 10*calc;
      switch(imove)
      {
         default:
         case 1: MoveChoice = SpinMove; break;
         case 2: MoveChoice = OccupancyMove; break;
      }
      switch( calc)
      {
      case 5:
         calc_energy(lsms_calc,sigma,MoveChoice);
         if (potentialShifter.vSpinShiftFlag)
         {
            generator -> updateHistogram(iwalk, &(sigma.evec[0]), &(sigma.vSpinShifts[0]), accepted);
         } 
         else 
         {
            generator -> updateHistogram(iwalk, &(sigma.evec[0]), accepted);
         }
         // dprintf("Master: Updating log after Walker %d move.\n",iwalk);
         generator->updateLog(iwalk, &(sigma.evec[0]), &(sigma.occ[0]), 
             r_values[energyIndex], accepted, prevMoveChoice, isSpinSim, isOccupancySim);
         std::cout << "slave " << my_group << " " << __FILE__ << ":" << __LINE__ << std::endl;
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
   generator->writeDOS(wl_outf.c_str());
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
  
     MPI_Send(&(site_alloyclass[0]), natom, MPI_INTEGER, 0, 1012, MPI_COMM_WORLD);
     free(site_alloyclass);
     // dprintf("Walker %d: Sent master alloy description.\n",my_group);
   }
}


void LSMS_Hamiltonian::recieve_alloy_description(EvecGenerator *generator)
{
    int nclasses;
    MPI_Recv(&nclasses, 1, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    std::vector<int> ncomps(nclasses);
    std::vector<int> site_alloyclass(natom);
    MPI_Recv(&(ncomps[0]), nclasses, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    alloyDesc.resize(nclasses);
    for(int i = 0; i < nclasses; i++) {
      alloyDesc[i].resize(ncomps[i]);
      MPI_Recv(alloyDesc[i].data(), ncomps[i]*sizeof(AtomType), MPI_BYTE, 
          lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    }
    MPI_Recv(&(site_alloyclass[0]), natom, MPI_INTEGER, lsms_rank0[0], 1012, MPI_COMM_WORLD, &status);
    generator->setAlloyClasses(alloyDesc, &(site_alloyclass[0]));
    printf("Master: Recieved alloy description.\n");
}


void LSMS_Hamiltonian::send_evec(LSMS_Config& sigma, int itarget, bool vSpinShiftFlag)
{
   if(vSpinShiftFlag)
   {
     for (int j=0; j<natom; j++)
     {
       sigma.evecsAndSpinShifts[4*j+0] = sigma.evec[3*j+0];
       sigma.evecsAndSpinShifts[4*j+1] = sigma.evec[3*j+1];
       sigma.evecsAndSpinShifts[4*j+2] = sigma.evec[3*j+2];
       sigma.evecsAndSpinShifts[4*j+3] = sigma.vSpinShifts[j];
     }
     MPI_Send(&(sigma.evecsAndSpinShifts[0]), 4*natom, MPI_DOUBLE, itarget /* lsms_rank0[i]*/, 51, MPI_COMM_WORLD);
   }
   else
     MPI_Send(&(sigma.evec[0]), 3*natom, MPI_DOUBLE, itarget /*lsms_rank0[i]*/, 51, MPI_COMM_WORLD);
}


void LSMS_Hamiltonian::send_occ(LSMS_Config& sigma, int itarget)
{
// Does the buffer have to persist after the call?
// std::vector<double> send_buffer(natom);
// for(int j = 0; j < natom; j++)
//   send_buffer[j] = double(sigma.occ[j]);
// MPI_Send( &(send_buffer[0]), natom, MPI_DOUBLE, itarget, 52, MPI_COMM_WORLD);
   MPI_Send( &(sigma.occ[0]), natom, MPI_DOUBLE, itarget, 52, MPI_COMM_WORLD);
}

void LSMS_Hamiltonian::do_without_master()
{
   LSMS_Config sigma;

   int running;
   int total_init_steps;
   bool accepted;

   MPI_Comm_rank(local_comm, &rank);
   snprintf(prefix, 38, "%d_", my_group);
   // to use the ramdisk on jaguarpf:
   // snprintf(prefix, 38, "/tmp/ompi/%d_", my_group);
   // This is from lsmsClass.hpp
   LSMS lsms_calc(local_comm, i_lsms_name, prefix, my_group);
   snprintf(prefix, 38, "Group %4d: ", my_group);

   // recieve either 'evec' or 'occupancy' variable set
   recv_buffer.resize(4*natom);
   sigma.occ.resize(natom);
   sigma.evec.resize(4*natom);
   r_values.resize( R_VALUE_OFFSET + 3*(natom+1) );
   i_values.resize(10);

}

void LSMS_Hamiltonian::do_master()
    {
      std::vector<LSMS_Config> sigma;

      int running;
      int total_init_steps;
      bool accepted;

      char *wl_inf = NULL;
      if (gWL_in_name) wl_inf = gWL_in_name;
      char *wl_outf = NULL;
      if (gWL_out_name)
      {
         wl_outf = gWL_out_name;
      }
      char *wl_stepf = NULL;
      if (step_out_flag && (evec_generation_mode == WangLandau_1d))
      {
        // step_out_flag = false;
        snprintf(wl_step_out_name, 127, "wl1d_%s", step_out_name);
        wl_stepf = wl_step_out_name;
      }
      std::cout << __FILE__ << ":" << __LINE__ << " wang-landau output=\"" << wl_outf << "\"" << std::endl;
        
      PotentialShifter potentialShifter;
      EvecGenerator *generator;

      sigma.resize(nwalk);
      for(int i=0; i<nwalk; i++)
      {
        sigma[i].evec.resize( 3*natom  );
        sigma[i].occ.resize(natom);
        sigma[i].init_steps = initial_steps;
      }
      if(potentialShifter.vSpinShiftFlag)
      {
        for(int i=0; i<nwalk; i++)
        {
          sigma[i].vSpinShifts.resize(natom);
          for (int j = 0; j < natom; j++) 
            sigma[i].vSpinShifts[j] = 0.0;
          sigma[i].evecsAndSpinShifts.resize(4*natom);
        }
      }
      send_buffer.resize( 4*natom );
      i_values.resize(10);

      total_init_steps = nwalk * initial_steps;
        
      r_values.resize( R_VALUE_OFFSET + 3*(natom+1) );

      //SMURF SMURF SMURF SMURF SMURF SMURF
      //browngrg 10/14/2015 Had to modify generators to deal with moving evecs,occs to LSMS_Config


      // Initialize the correct evec generator
      switch (evec_generation_mode)
      {
      case Random : 
         std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         generator = new RandomEvecGenerator(natom);
         break;
      case Constant:
         generator = new ConstantEvecGenerator(natom, ev0, nwalk);
         break;
      case WangLandau_1d : 
         //generator = new WL1dEvecGenerator<std::mt19937>(natom, nwalk, evecs, potentialShifter, wl_inf, wl_outf, wl_stepf, occs);
         generator = new WL1dEvecGenerator<std::mt19937>(natom, nwalk, NULL, potentialShifter, wl_inf, wl_outf, wl_stepf, NULL);
         break;
      case ExhaustiveIsing : 
         //generator = new ExhaustiveIsing1dEvecGenerator(natom, nwalk, evecs, wl_inf, wl_outf);
         generator = new ExhaustiveIsing1dEvecGenerator(natom, nwalk, NULL, wl_inf, wl_outf);
         break;
      case WangLandau_2d : 
         //generator = new WL2dEvecGenerator<std::mt19937>(natom, nwalk, evecs, wl_inf, wl_outf, wl_stepf);
         generator = new WL2dEvecGenerator<std::mt19937>(natom, nwalk, NULL, wl_inf, wl_outf, wl_stepf);
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
        for(int i=0; i<nwalk; i++)
          generator -> initializeEvecAndPotentialShift(i, &(sigma[i].evec[0]), &(sigma[i].vSpinShifts[0]));
      else 
        for(int i=0; i<nwalk; i++)
          generator -> initializeEvec(i, &(sigma[i].evec[0]) );

      // generate initial occupancies
      if( isOccupancySim ) {
        printf("Master: Generating initial occupancies.\n");
        for(int i = 0; i < nwalk; i++)
          generator->initializeOccupancies(i,&(sigma[i].occ[0]));
      }

      // if both spin and occupancy variables are needed
      // then we need to set both sets of variables before beginning calculations
      if( isSpinSim && isOccupancySim )
      for(int i = 0; i < nwalk; i++) {
        send_evec(sigma[i],lsms_rank0[i],potentialShifter.vSpinShiftFlag);
        send_occ(sigma[i],lsms_rank0[i]);
      }
    
      std::cout << __FILE__ << ":" << __LINE__ << " This is the master node\n";

      // issue initial commands to all LSMS instances
      running = 0;
      bool more_work = true;
      if (total_init_steps > 0)
      {

       
        for(int i=0; i<nwalk; i++)
        {
          std::cout << __FILE__ << ":" << __LINE__ << "starting initial calculation in group " << i << std::endl;

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
          // get results from last calc_energy push
          if (return_moments_flag)
            MPI_Recv( &(r_values[0]), R_VALUE_OFFSET+3*natom, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          else
            MPI_Recv( &(r_values[0]), R_VALUE_OFFSET, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          running--;
         
          // std::cout << "received energy E_tot =" << r_values[0] << std::endl;
          // std::cout << "    band energy E_band =" << r_values[1] << std::endl;
          printf("received energy E_tot = %25.15f\n",r_values[0]);
          printf("    band energy E_band= %25.15f\n",r_values[1]);
 
          if (total_init_steps > 0)
          {
            int iwalk = (status.MPI_SOURCE - align) / comm_size;
            std::cout << __FILE__ << ":" << __LINE__ << "starting additional calculation in group " << iwalk << std::endl;

            MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);
            if (sigma[iwalk].init_steps > 0)
            {
              if( MoveChoice == SpinMove ) 
                generator -> generateUnsampledEvec(iwalk, &(sigma[iwalk].evec[0]), r_values[energyIndex]);
              else if( MoveChoice == OccupancyMove ) 
                generator -> generateUnsampledOcc(iwalk, &(sigma[iwalk].occ[0]));
              //more_work = !(generator -> generateUnsampledEvec(iwalk, evecs[iwalk], r_values[energyIndex]));
              sigma[iwalk].init_steps;
              total_init_steps--;
            }

            // send configuration
            if( MoveChoice == SpinMove )  { 
              // dprintf("Master: Sending trial spins to Walker %d.\n", iwalk);
              send_evec(sigma[iwalk],lsms_rank0[iwalk],potentialShifter.vSpinShiftFlag);
            }
            else if( MoveChoice == OccupancyMove ) {
              // dprintf("Master: Sending trial occupancies to Walker %d.\n", iwalk);
              send_occ(sigma[iwalk],lsms_rank0[iwalk]);
            }
                  
            num_steps--;
            running++;
            stepCount++;
            walkerSteps[iwalk+1]++;
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
      for (int i=0; i<nwalk; i++)
      {
        std::cout << __FILE__ << ":" << __LINE__ << " starting main calculation in group " << i << std::endl;
        // select whether to perform spin or occupancy move
        MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
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
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      }
        
      generator -> startSampling(isSpinSim, isOccupancySim);

      // define variables to keep track whether spin (or occupancy)
      // configuration needs to be reverted to previous values
      std::vector<bool> acceptedSpinMove(nwalk);
      std::vector<bool> acceptedOccMove(nwalk);

      for(int i = 0; i < nwalk; i++)
        acceptedSpinMove[i] = acceptedOccMove[i] = true;

      // wait for results and issue new commands or wind down
      while (running > 0)
      {
        accepted = false;

      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        MPI_Recv( &(r_values[0]), R_VALUE_OFFSET+3*natom, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        running--;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        // determine whether returning from a spin or occupancy move
        if( r_values[2] == 0x00 ) { prevMoveChoice = SpinMove; } else if( r_values[2] = 0x01 ) { prevMoveChoice = OccupancyMove; }
        int iwalk = (status.MPI_SOURCE - align) / comm_size;
        printf("received energy E_tot = %25.15f\n",r_values[0]);
        printf("    band energy E_band= %25.15f\n",r_values[1]);
        // printf("from status.MPI_SOURCE=%d\n",status.MPI_SOURCE);
        energy_accumulator += r_values[0];
        energies_accumulated++;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;

        if (more_work)
        {
          std::cout << "starting additional calculation in group " << iwalk << std::endl;
              
          if (generator_needs_moment)
          {
            double m0, m1, m2;
            m0 = 0.0; m1 = 0.0; m2 = 0.0;
            for(int i=0; i<3*natom; i+=3)
            {
              m0 += r_values[R_VALUE_OFFSET+i];
              m1 += r_values[R_VALUE_OFFSET+i+1];
              m2 += r_values[R_VALUE_OFFSET+i+2];
            }
            switch (second_dimension)
            {
              case  MagneticMoment :  magnetization = std::sqrt(m0*m0 + m1*m1 + m2*m2); break;
              case  MagneticMomentX : magnetization = m0; break;
              case  MagneticMomentY : magnetization = m1; break;
              case  MagneticMomentZ : magnetization = m2; break;
            }

            // todo: separately keep track of whether last evec or occ move accepted
            // Determine if the configuration is accepted
            accepted = generator -> determineAcceptance(iwalk, r_values[energyIndex], magnetization);
          } 
          else {
            // Determine if the configuration is accepted
            accepted = generator -> determineAcceptance(iwalk, r_values[energyIndex]);
          }

          if(accepted)
            printf("Master: Accepted Walker %d trial move.\n",iwalk);
          else
            printf("Master: Rejected Walker %d trial move.\n",iwalk);

           // need to separate track of whether the last spin or occupancy change was accepted
           // this ensures the Wang-Landau generator can properly restore the old configuration
          if( prevMoveChoice == SpinMove ) 
            acceptedSpinMove[iwalk] = accepted;
          else if( prevMoveChoice == OccupancyMove )
            acceptedOccMove[iwalk] = accepted;

          // dprintf("Master: Updating histogram after Walker %d move.\n",iwalk);
          if (potentialShifter.vSpinShiftFlag)
          {
            if (generator -> updateHistogram(iwalk, &(sigma[iwalk].evec[0]), &(sigma[iwalk].vSpinShifts[0]), accepted))
              more_work = false;
          } 
          else 
          {
            if (generator -> updateHistogram(iwalk, &(sigma[iwalk].evec[0]), accepted))
              more_work = false;
          }
          // dprintf("Master: Updating log after Walker %d move.\n",iwalk);
          generator->updateLog(iwalk, &(sigma[iwalk].evec[0]), &(sigma[iwalk].occ[0]), 
              r_values[energyIndex], accepted, prevMoveChoice, 
              isSpinSim, isOccupancySim);

          // Prepare a new configuration
          MoveChoice = generator->selectMoveType(isSpinSim, isOccupancySim);

          if( MoveChoice == SpinMove ) {
            // dprintf("Master: Generating trial spins for Walker %d.\n",iwalk);
            generator -> generateEvec(iwalk, &(sigma[iwalk].evec[0]), acceptedSpinMove[iwalk]);
            if (potentialShifter.vSpinShiftFlag)
              generator -> generatePotentialShift(iwalk, &(sigma[iwalk].vSpinShifts[0]), acceptedSpinMove[iwalk]);
          }
          else if( MoveChoice == OccupancyMove ) {
            // dprintf("Master: Generating trial occupancies for Walker %d.\n",iwalk);
            generator->generateOccupancies(iwalk, &(sigma[iwalk].occ[0]), acceptedOccMove[iwalk]);
          }

          if( MoveChoice == SpinMove ) 
          { 
            send_evec(sigma[iwalk],lsms_rank0[iwalk],potentialShifter.vSpinShiftFlag);
          }
          else if( MoveChoice == OccupancyMove ) 
          {
            send_occ(sigma[iwalk],lsms_rank0[iwalk]);
          }
  
          num_steps--;
          running++;
          stepCount++;
          walkerSteps[iwalk+1]++;
          if (restrict_steps && num_steps <= 0) more_work = false;
          if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
          walltime = MPI_Wtime() - walltime_0;
          if (restrict_time && walltime >= max_time) more_work = false;
          if (restrict_time) std::cout << "      " << max_time - walltime << " seconds remaining\n";
        }
        else
        {
          // send an exit message to this instance of LSMS
          int iwalk = (status.MPI_SOURCE - align) / comm_size;

          MPI_Send( &(sigma[iwalk].evec[0]), 3*natom, MPI_DOUBLE, lsms_rank0[iwalk], 2, MPI_COMM_WORLD);
        }

        if (step_out_flag && accepted)
        {
          step_out_file << "# iteration " << energies_accumulated << std::endl;
          step_out_file.precision(15);
          step_out_file << energies_accumulated << std::endl;
          step_out_file << r_values[0] << "  " << r_values[1] << std::endl;
          for (int j=0; j<3*natom; j+=3)
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
     // std::vector<long int> rb(2*(nwalk+1));
     std::vector<long int> rb(2*size);
     MPI_Gather(sb,2,MPI_LONG,&rb[0],2,MPI_LONG,0,MPI_COMM_WORLD);
     walkerSteps[0]=stepCount;  // report the total step count for the WL master
     std::ofstream statOut; statOut.open("energyLoopCount.statistics");
       statOut<<rb[2*0]<<"   "<<walkerSteps[0]<<"   "<<rb[2*0+1]<<std::endl;
     for(int i=0; i<nwalk; i++)
       statOut<<rb[2*lsms_rank0[i]]<<"   "<<walkerSteps[i+1]<<"   "<<rb[2*lsms_rank0[i]+1]<<std::endl;
     statOut.close();
/*
  if(evec_generation_mode==WangLandau_1d)
  (static_cast<WL1dEvecGenerator<std::mt19937> *>(generator))->writeState("WLrestart.state");
  if(evec_generation_mode==ExhaustiveIsing)
  (static_cast<ExhaustiveIsing1dEvecGenerator *>(generator))->writeState("WLrestart.state");
*/

   generator->writeDOS(wl_outf);

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
  lsms_rank0.resize(nwalk+1) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  world_rank = rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Construct new MPI info and compare to old
  // Not compatible with Master/Slave
  //int NWalk = nwalk;
  //MPI_Struct mp_world = MPI_Struct::world();
  //MPI_Gang mp_gang;
  //std::cout << __FILE__ << ":" << __LINE__ << " NWalk=" << NWalk << std::endl;
  //mp_gang.init(NWalk,mp_world);
  //if( mp_gang.igang!=my_group ) std::cout << "**MPI(" << mp_world.iproc << "): mp_gang.igang=" << mp_gang.igang << " my_group=" << my_group << std::endl;

  walltime_0 = MPI_Wtime();

  // PAPI START

  write_header();


  if (generator_needs_moment) return_moments_flag = true;

  if (nwalk <= 1)
  {
     std::cout << "nwalk==1 known to not work" << std::endl;
     exit(1);
  }
  walkerSteps.resize(nwalk+1); for(int i=0; i<nwalk+1; i++) walkerSteps[i]=0;
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
