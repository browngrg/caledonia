Build with: make EHModel2.exe
Run with EHModel2.exe -i EHModel-ProgramOptions.txt
Get a brief description of options with EHMdole2.exe -h  or EHModel2.exe --help
You can use the ProgramOptionsOut.txt from one as input for the next run.

Here is an explanation of the parameters in the input file

help= 0

   Setting this to 1 generates the same output as the -h or --help flags

verbose= 0

   Setting this to 1 generates more output. Its really just debugging 
   information so keep this set to 0

input= EHModel-ProgramOptions.txt

   This specifies the input file, like the -i or --input flags

output= ProgramOptionsOut.txt

   This is the filename used to output the complete list of options,
   along with useful run information

H= 0

   This is the applied magnetic field

ehfile= EHModel-Ni2MnGa.txt

   This is the file of exchange interactions

ehfilet= EHModelPairs

   This describes the exchange interaction file. The choices are
   EHModelPairs and EHModelKiJij. The second is for anistropy problems.

Elo= -1400

   The global lowest energy in the calculation. Individual windows will
   have different values

Ehi= 0

   The global highest energy in he calculation.

Ebin= 4

   This is the width of the bin in energy units. It should be the
   smaller than the average Monte Carlo step. If Ebin is bigger than
   the MC the Transition Matrix estimate of the entropy is no longer
   quantitative (the shape is right, but there is an unknown 
   multiplicative factor). If Ebin gets too small causes problems, too.

numwalk= 16

   Number of walkers per MPI process. If nproc>1, then each processor
   get assigned to one window. Then the number of walkers in each
   window is the product of numwalk times the number of processors
   assigned to each window. If the run involves only one process,
   particularly when not in MPI, then this is the total number
   of walkers. Then numwalk needs to be a multiple of the
   number of windows.

overlap= 0.5

   This is the fractional overlap between adjacent windows. A value
   of 0.5 means that the windows overlap by half, and that each 
   energy is covered by two windows (except the very highest and
   lowest energies). For replica exchange to work well this needs
   to be a large fraction, between 0.5 and 0.8

numbin= 200

   Number of energy bins in each window. If this number is too small
   you lose resolution. If the number is too big, the bins are small
   and the jumps in the transition matrix are too big. You want the 
   bins to be about the change in energy for each Monte Carlo move.

Q= 0.05

   This is the convergence criterion, essentilly the flatness
   parameter. Here, however, Q is calculated as the variation
   of the visit parameter from the normalized mean. Values of
   Q very close to zero mean a very flat visit histogram, and
   values far from zero indicate lots of fluctuations or large
   systematic variation. Values between 0.2 and 0.01 are reasonable.
   The calculation resets the visit histogram periodically,
   otherwise the systematic corrections at the beginning of a
   stage of convergence take a long time to become smaller than
   the total number of counts. The number of steps between
   resets only grows. At any stage it will be at least as
   many steps as it took to find a flat visit histogram at
   the previous level. After that, it checks for flatness after
   every nchunk steps.

nchunk= 1000000

   This is the number of Monte Carlo steps to do between output and
   checks for convergence. The number of steps between zeroing the
   Wang-Landau visit histogram grows by nchunk, too. 

maxupdate= 15

   Maximum number of changes in the Wang-Landau parameter before
   the program stops. Some windows will get ahead of others. The
   program will not exit until all windows have reached 
   iupdate=maxupdate. The other windows will push ahead because
   replica exchange needs them to continue.

update= 0

   Set this to 0 if you want the input DOS from estdos (see
   below) to stay the same throughout the entire run. Set this
   parameter to 1 if you want the Transition Matrix results to
   update this as the run continues.

bywindow= 0

   Set this to 1 if you want more output: DOS broken down by
   window, Walker configurations, and the Transition Matrix
   data that can be reloaded by subsequent runs.

useITTM= 0

   Set this to 1 if you want to use Transition Matrix 
   data from previous runs.

estdos= SeedDOS.csv

   This is a file with the initial value of Sfixed. It needs
   to cover at least the range Elo to Ehi. If it is bigger,
   then only the needed values are retained. Linear interpolation
   is used to convert to the energies in the calculation. The
   file format is a simple x,y csv file with energy in the first
   column and lng in the second column. Extra columns are
   ignored.

Greg Brown (browngrg@comcast.net,gbrown@fsu.edu)
March 9, 2015

