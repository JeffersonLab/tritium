#!/bin/tcsh
source /u/apps/root/6.10.02//setroot_CUE.csh
setenv ANALYZER /work/halla/triton/MARATHON_Replay/analyzer1.6/
setenv DB_DIR /work/halla/triton/MARATHON_Replay/DB/

setenv LD_LIBRARY_PATH ${ANALYZER}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${ANALYZER}:${PATH}
