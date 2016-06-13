import numpy as np
import gzip as g
import tempfile
import sys
import os



#INPUTS
# sys.argv[1] -- striped input isochrone file leaving only one line starting with # col_names
# sys.argv[2] -- map of columns names in output
# sys.argv[3] -- output file name

# clean input file first
# CLEANER PART OF THE SCRIPT

tmp = tempfile.NamedTemporaryFile()
os.system('sed -i "s/stage//" '+tmp.name)
bad=False

fout=open(tmp.name, 'w')
for line in open(sys.argv[1], 'r'):
    #if not comments, remove possible last stage information
    if line.split()[0] != '#':
        spline=line.split()[0:21]
    else: # do not limit if in header
        spline=line.split()
        
    if line.split()[-1] != 'LTP':
        if not bad:
            for i in spline:
                print>>fout,  i, 
            print >>fout,''
        elif bad and line.split()[0] == '#': # stop of deleting
            bad=False
            for i in spline:
                print >>fout, i, 
            print >>fout,''
    else: # start of deleting
        bad=True
fout.close()
# ==== END OF CLEANER PART ====================

fiso=open(tmp.name, 'r')
    
for line in fiso:
    if line.split()[1] == "Isochrone":
        METALICITY=line.split()[4]
        break
print METALICITY

data=np.genfromtxt(fiso,names=True)
print data.dtype.names

ages, index=np.unique(data['logageyr'], return_index=True)
file=open(sys.argv[3], "w")
columns = [line.split()[0] for line in open(sys.argv[2], "r")]
print >>file, ages.size
for i in range(ages.size):
    i_start = index[i]
    if i != len(index)-1:
        i_stop = index[i+1]
    else:
        i_stop = data['logageyr'].size
    #print >>file,">", i_stop-i_start, METALICITY, 10**ages[i], len(columns)
    print >>file,">", i_stop-10-i_start, METALICITY, 10**ages[i], len(columns) # TO fix strange AGB stars sequence which is most probable error in isochrones
    #for j in range(i_start, i_stop, 1):
    for j in range(i_start, i_stop-10, 1): # TO fix strange AGB stars sequence which is most probable error in isochrones
        for col in columns:
            if col == "M_ini" or col == "M_act":
                print >>file, np.log10(data[col][j]), 
            else:
                 print >>file, data[col][j], 
        print >>file, ''
tmp.close()
