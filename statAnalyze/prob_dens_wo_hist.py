#!/usr/bin/env python 

# generate probability density function without histogram. 
# ref . "From data to probability densities without histograms", 2008, CPC. 


# the programe takes input from stdin and write to stdout 

import numpy as np 
import sys 
import scipy 
import scipy.stats as stats
from  scipy.special import kolmogorov
import  matplotlib.pyplot as plt

M = 1000
Qcut = 0.8 
xmin = 0
xmax = 0 

def kvtest(cdf1, cdf2):
    n1 = len(cdf1)
    n2 = len(cdf2)
    if n1 != n2:
        print >> sys.stderr , "Wrong" 
        sys.exit(1)
    
    dn = max( np.abs( cdf1 - cdf2 )  )  
    dn = dn * np.sqrt( n1 ) 
    
    q = kolmogorov(dn)

    return q 


def probDensWoHist( data ) :
    '''
data is a list of raw data 
    '''
    data = np.array( data, dtype = float ) 

    data.sort()
    
    # plot 
    f1 = plt.figure()
    ax1 = f1.add_subplot(111)
    ax1.hist(data, 100)


    #n = len(data) 
    
    # accumulated cdf 
    data_nodup = []
    Fab = []
    i=0 
    ndata = len( data) 
    while i < ndata:  
        j = i+1
        while j < len(data):
            if data[j] == data[i]:
                j += 1
            else:
                break 
        
        data_nodup.append(data[i])
        Fab.append(  float(j)   / ndata  )
        i = j 


 
    data_nodup = np.array( data_nodup , dtype = float )
    Fab = np.array( Fab , dtype =float ) 
    
    data = data_nodup

    A = data[0]
    B = data[-1]
    L = (B-A) /2.0
    Linv = 1.0 / L 

    # ref function 
    N = len(data)
    F0 = np.array( [ (1.0 - 1.0/N ) / ( B - A ) * ( i- A )  + 1.0/N for i in data ] , dtype = float ) 

    # residue 
    R = Fab - F0

    
    R0 = np.ones_like( data )

    diList = {}
    
    dicos0 = 1/L * scipy.integrate.simps(  R , data   )
    R0 = R0 * dicos0  / 2.0

    diList[0] = dicos0

    for i in range(1,M+1):
        #print data 
        #print R *  scipy.sin( i * scipy.pi * ( data  - a )/ (b-a) )


        disin =  Linv * scipy.integrate.simps(  R *  scipy.sin(  i * scipy.pi *  data  * Linv ) , data   )
        dicos =  Linv * scipy.integrate.simps(  R *  scipy.cos(  i * scipy.pi *  data  * Linv ) , data   )
        R0 += disin * scipy.sin( Linv * i * scipy.pi * data   )
        R0 += dicos * scipy.cos( Linv * i * scipy.pi * data   )

        # convert R0 back to cdf and use Kolmogorov test for convergence. 
        tmp = R0 + F0
        

        diList[i]= (dicos, disin ) 

        #if i==3:
        #    plt.plot( data, tmp, data, Fab  )
        #    plt.show()
        #    sys.exit(1)

        tmpq =  kvtest(tmp, Fab ) 

        if tmpq > Qcut:
            break 
    
    if i == M:
        print >> sys.stderr, "M too small, or increase Qcut "
        sys.exit(1)

    tmpd = np.diff(tmp ) 
    tmpdat = np.diff(data) 
    
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)
    ax2.plot( data[1:], tmpd / tmpdat) 
    plt.show()
    sys.exit(1)

    # get probablity density 
    xmin = A 
    xmax = B 
    
    print i, diList , A, B

    x = np.linspace(xmin, xmax, 100)
    print x 
    probDens = np.zeros_like( x ) 
    for n in range(1, i+1):
        an, bn = diList[n]
        probDens = - an * Linv * n * scipy.pi * scipy.sin( n * scipy.pi * Linv * x ) \
                   + bn * Linv * n * scipy.pi * scipy.cos( n * scipy.pi * Linv * x )

    probDens += ( 1- 1.0/N ) / (B - A) 

    plt.plot( x, probDens )
    plt.show()

if __name__ == '__main__':
    

    data = [ i.strip() for i in sys.stdin.readlines()  ]

    probDensWoHist(data) 
