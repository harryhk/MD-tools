#!/bin/bash
LO=8 
LOM=8
LNL=7  
LC=6 
LH1=7
LH2=8    
LP=15 
LOS=8
LP2=8  
LP3=9 
LC3=9    
LC2=8    
OW=8
HW=1




for i in `sed '7,60 !d' dopc-fixed.itp  | awk '{print $2}'`; do echo ${i} | awk '{printf "%5s", $1}' ; done
echo
for i in `sed '7,60 !d' dopc-fixed.itp  | awk '{print $2}'`; do echo ${!i} | awk '{printf "%5d", $1}';  done
echo
