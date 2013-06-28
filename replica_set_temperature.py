#!/bin/env python 

# calculate new replica exchange temperature according to the old temperature and exchange rate 


import sys, re
from lnx_util import parseInput, print_help
import subprocess


inputP= parseInput(sys.argv[1:])
paraOpt = '-flog -h'.split()

helpdoc=\
'''
Usage! ./prog.py
        -flog run1_0.log  ; log file from any replica
        -h   ; print help information 
'''

print_help(inputP, paraOpt, helpdoc)
finm = inputP['-flog']

rate = 0.25

exe = 'grep \'Repl  fake T\' run1_0.log '
p= subprocess.Popen( exe,  shell=True, stdout=subprocess.PIPE)
temperature = [ float(i) for i in   re.findall('\d*\.*\d*', p.stdout.readline() ) if len(i) >0 ]

exe = 'grep  -A 2 \'Repl  average probabilities:\' run1_0.log  | tail -n 1 '
p= subprocess.Popen( exe,  shell=True, stdout=subprocess.PIPE)
prob = [ float(i) for i in   re.findall('\d*\.*\d*', p.stdout.readline() ) if len(i) >0 ]

new_temp = [ 0 for i in range(len(temperature))  ]
new_temp[0]= 323.0

for i  in range(1, len(temperature)):
    new_temp[i] = new_temp[i-1] + ( temperature[i] - temperature[i-1] ) * prob[i-1] / rate

print temperature
print prob
print map(int, new_temp)
    
