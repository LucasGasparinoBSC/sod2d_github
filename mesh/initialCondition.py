import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
velofile        = 'VELOC.alya'
pressfile        = 'PRESS.alya'
densifile        = 'DENSI.alya'

fCoord          =open(coordfile,'r')
fVel            =open(velofile,'w')
fPress           =open(pressfile,'w')
fDensi           =open(densifile,'w')

pi = 3.14159265
c = 0.0

print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])

    if x < c:
        vx = 0.0
        vy = 0.0
        pr = 1.0
        rho = 1.0
    elif x>= c:
        vx = 0.0
        vy = 0.0
        pr = 0.1
        rho = 0.125

    fVel.write('{} {} {}\n'.format(pid,vx,vy))
    fPress.write('{} {}\n'.format(pid,pr))
    fDensi.write('{} {}\n'.format(pid,rho))
        
fCoord.close()
fVel.close()
fPress.close()
fDensi.close()

print('---| End writing initial condition')
