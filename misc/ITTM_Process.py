#!/usr/bin/python3

from subprocess import call,DEVNULL
from os import listdir
import re

all_files = listdir(".")
level1_search = re.compile("WangLandau-01-\d*\.csv",re.IGNORECASE)
level1_files = [ f for f in all_files if level1_search.search(f) ]
print("Level1 files = {0}  = {1}".format(len(level1_files),level1_files))
levels_search = re.compile("WangLandau-\d*-00\.csv")
levels_files = [ f for f in all_files if levels_search.search(f) ]
print("Levels files = {0}  = {1}".format(len(levels_files),levels_files))

NLevel = len(levels_files)+1
NWin = len(level1_files)
Rebin = False
Ebin = 4
Emin = [ -15000 + 2500*iwin + 1 for iwin in range(NWin) ]

for iwin in range(NWin):
   print("win={0} start={1}".format(iwin,Emin[iwin]))

for ilevel in range(1,NLevel):
   if Rebin:
      for iwin in range(NWin):
         cin  = open("WangLandauWalkerC-{0:02d}-{1:02d}.csv".format(ilevel,iwin),"r")
         cout = open("WangLandauWalkerC2-{0:02d}-{1:02d}.csv".format(ilevel,iwin),"w")
         call(["../../src/ITTM_CRebin.exe"],stdin=cin,stdout=cout)
         cin.close()
         cin.close()
         cin = open("WangLandauWalkerC2-{0:02d}-{1:02d}.csv".format(ilevel,iwin),"r")
         cout = open("ITTMDOS-{0:02d}-{1:02d}.csv".format(ilevel,iwin),"w")
         call(["../../src/ITTM_DOS.exe","{}".format(Ebin),"{}".format(Emin[iwin])],stdin=cin,stdout=cout)
         cin.close()
         cin.close()
      files = [ "ITTMDOS-{0:02d}-{1:02d}.csv".format(ilevel,iwin) for iwin in range(NWin) ]
      call(["../../src/DOSStitch3.exe","2"]+files,stdout=DEVNULL)
      call(["mv","DOSStitchResult.csv","DOS-{0:02d}-ITTM.csv".format(ilevel)])
   else:
      files = [ "WangLandau-{0:02d}-{1:02d}.csv".format(ilevel,iwin) for iwin in range(NWin) ]
      call(["../../src/DOSStitch3.exe","7"]+files,stdout=DEVNULL)
      call(["mv","DOSStitchResult.csv","DOS-{0:02d}-ITTM.csv".format(ilevel)])
   files = [ "WangLandau-{0:02d}-{1:02d}.csv".format(ilevel,iwin) for iwin in range(NWin) ]
   call(["../../src/DOSStitch3.exe","2"]+files,stdout=DEVNULL)
   call(["mv","DOSStitchResult.csv","DOS-{0:02d}-WL.csv".format(ilevel)])


