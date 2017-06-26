#!/bin/bash

LIST="22013 22011 22008 22007 22006"
N=2000000

for i in $LIST
do
  echo "$i"
  : > right_gmp_$i.root.DpCut.ncut
  
    if [[ -e "right_gmp_$i.root.DpCut.cut" ]]; then
      while read l
      do
        if [[ "$l" == "fEvtHdr.fRun==0" ]];then
          echo -e "0" >> right_gmp_$i.root.DpCut.ncut 
        else
          echo -e "$N" >> right_gmp_$i.root.DpCut.ncut
        fi
      done < right_gmp_$i.root.DpCut.cut
    else  
        echo -e "fEvtHdr.fRun==0" >> right_gmp_$i.root.DpCut.cut
        echo -e "0" >> right_gmp_$i.root.DpCut.ncut    
    fi
    echo -e  "" >> right_gmp_$i.root.DpCut.cut 
    echo -e  "" >> right_gmp_$i.root.DpCut.ncut  
done
