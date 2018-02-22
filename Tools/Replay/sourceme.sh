#!/bin/bash
source /apps/root/6.10.02/setroot_CUE.bash
export ANALYZER=/work/halla/triton/MARATHON_Replay/analyzer1.6

export DB_DIR=/work/halla/triton/MARATHON_Replay/DB

export LD_LIBRARY_PATH=${ANALYZER}::$ROOTSYS/lib:${LD_LIBRARY_PATH}
export PATH=${ROOTSYS}/bin:${ANALYZER}:${PATH}
