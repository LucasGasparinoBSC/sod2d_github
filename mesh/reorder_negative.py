# Reorders elements with negative Jacobian
# Only for QUA04

import numpy as np
import array as ar

# GLL table with 1 pt.:

xgp = np.zeros((1,2))

xgp[0,0] = 0.00
xgp[0,1] = 0.00

# Derivative of shape function at GP:

s = xgp[0,0]
t = xgp[0,1]

dN = np.zeros((2,4))

dN[0,0] = -1.00+t
dN[1,0] = -1.00+s
dN[0,1] =  1.00-t
dN[1,1] = -1.00-s
dN[0,2] =  1.00+t
dN[1,2] =  1.00+s
dN[0,3] = -1.00-t
dN[1,3] =  1.00-s

# Read file *.geo.dat

oldfile="shock_tube.geo.dat"
newfile="reordered.geo.dat"

f=open(oldfile,'r')

line = True
inELEM = False
inCOOR = False

nelem = 120 # From mesh info
npoin = 147 # Idem
nnode = 4     # QUA04
connec = np.zeros((nelem,nnode)) # mesh table
coord = np.zeros((npoin,2))      # coordinates

print('Reading data...')
while line:
      line = f.readline()
      if 'END_SKEW' in line:
         line = False
      else:
         if 'ELEMENTS' in line:
            if 'END_ELEMENTS' in line:
               print('Finished element data!')
               inELEM = False
               line = True
            else:
               inELEM = True
               print('Reading element data...')
               line = True
         elif 'COORD' in line:
            if 'END_COORD' in line:
               print('Finished coordinate data!')
               inCOOR = False
               line = True
            else:
               inCOOR = True
               print('Reading coordinate data...')
               line = True
         else:
            if inELEM:
               stripline = line.strip().split()
               s = [int(i) for i in stripline]
               connec[s[0]-1,:] = s[1:]
               #print(connec[s[0]-1,:].astype(int))
               line = True
            elif inCOOR:
               stripline = line.strip().split()
               s = [float(i) for i in stripline]
               coord[int(s[0])-1,:] = s[1:]
               #print(coord[int(s[0]-1),:])
               line = True
            else:
               line = True

f.close()

# Start check

elcod = np.zeros((2,4))
Je = np.zeros((2,2))
for ielem in range(nelem):
    indList = connec[ielem,:].astype(int)-1
    elcod[0,:] = coord[indList,0]
    elcod[1,:] = coord[indList,1]
    Je = np.matmul(elcod,np.transpose(dN))
    detJe = Je[0,0]*Je[1,1]-Je[1,0]*Je[0,1]
