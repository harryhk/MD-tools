#!/usr/bin/evn python 

import sys
from common.lnx_util import parseInput, print_help 


if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -o -h '.split(' ')

    helpdoc = ''\
              'Usage prog.py  -f out.0  out.1 ...  ; diff outfiles from embrassing parallel\n'\
              '               -o out.xvg   ; concated output file \n'
    
    print_help(inputP, paraOpt, helpdoc)

    fout = open(inputP['-o'], 'w')

    result ={} # super big hash table to join result 
    
    flag=True

    for fn in inputP['-f']:
        
        fin = open(fn)
        # read files 
        header = fin.readline().strip().split()[1:]
        nheader = len( header) 

        for i in header:
            if not result.has_key(i):
                result[i]=[]
        
        while True:
            line = fin.readline().strip()
            if not line:
                break

            line = line.split() 
            nline = len( line)  
            
            if nline != nheader:
                print >> sys.stderr, "data format wrong!"
                sys.exit(1)

            if flag:
                
                for i in range( nline ):
                    result[ header[i] ].append(line[i] )


            else:
                for i in range(1, len(line) ) : # skip time column
                    result[ header[i] ].append( line[i] )


        flag = False

    # put together 
    index = sorted( result.keys() )[:-1]
    index.sort( key= lambda x : map( int, x.split('-') ) )

    # print header 
    fout.write('%15s' % '# time')
    for i in index:
        fout.write('%15s' % i)

    fout.write('\n')
    
    # write data
    n = len( result['time'] )
    print n
    for i in range(n):
        fout.write('%15s' % result['time'][i])
        for ii in index:
            fout.write('%15s' % result[ii][i] )
        fout.write('\n')

