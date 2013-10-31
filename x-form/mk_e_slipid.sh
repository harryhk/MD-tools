#!/bin/bash

# since it is all-atom model. Atom types can be determined from first letter. 

C=6 
H=1
O=8
P=15
N=7

itpfile=$1

if [ -z $1  ]; then 
	exit 1
fi 


for i in `cat $itpfile`; do echo -n "${i}   "  ; done
echo 
for i in `cat $itpfile`; do  j=${i:0:1}; echo -n "${!j}   "  ; done
echo 



#for i in `sed '36,173 !d' $itpfile  | awk '{print $2}'`; do echo ${i} | awk '{printf "%5s", $1}' ; done
#echo
#for i in `sed '36,173 !d' $itpfile  | awk '{print $2}'`; do  j=${i:0:1}; echo ${!j} | awk '{printf "%5d", $1}';  done
#echo
