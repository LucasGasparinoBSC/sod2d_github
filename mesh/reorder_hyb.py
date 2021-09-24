###################################################################
# Copyright 2005 - 2021 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



# Use after gmsh2alya for reordering HEX27 elements

import array as ar

oldfile="shock_tube.geo.dat"
newfile="replace.geo.dat"

f=open(oldfile,'r')
g=open(newfile,'w')

line = True
inELEM = False

while line:
      line = f.readline()
      if 'END_SKEW' in line:
         g.write(line)
         print('End Writting')
         print(line)
         line = False
      else:
         if 'ELEMENTS' in line:
            if 'END_ELEMENTS' in line:
               g.write(line)
               print('Finished element data!')
               inELEM = False
               line = True
            else:
               inELEM = True
               g.write(line)
               print('Reading element data...')
               line = True
         else:
            if inELEM:
               stripline = line.strip().split()
               s = [int(i) for i in stripline]
               nnodes = len(s)-1
               #print('nnodes = ', nnodes)

               if nnodes == 27:
                  #print('Reordering HEX27!')

                  a = [0] * 16
                  a[0] = s[12]
                  a[1] = s[14]
                  a[2] = s[10]
                  a[3] = s[17]
                  a[4] = s[19]
                  a[5] = s[20]
                  a[6] = s[18]
                  a[7] = s[11]
                  a[8] = s[13]
                  a[9] = s[15]
                  a[10] = s[16]
                  a[11] = s[23]
                  a[12] = s[24]
                  a[13] = s[22]
                  a[14] = s[25]
                  a[15] = s[21]

                  s[10] = a[0]
                  s[11] = a[1]
                  s[12] = a[2]
                  s[13] = a[3]
                  s[14] = a[4]
                  s[15] = a[5]
                  s[16] = a[6]
                  s[17] = a[7]
                  s[18] = a[8]
                  s[19] = a[9]
                  s[20] = a[10]
                  s[21] = a[11]
                  s[22] = a[12]
                  s[23] = a[13]
                  s[24] = a[14]
                  s[25] = a[15]

               elif nnodes == 18:
                  #print('Reordering PEN18!')

                  a = [0] * 7
                  a[0] = s[10]
                  a[1] = s[8]
                  a[2] = s[9]
                  a[3] = s[15]
                  a[4] = s[14]
                  a[5] = s[18]
                  a[6] = s[17]

                  s[8] = a[0]
                  s[9] = a[1]
                  s[10] = a[2]
                  s[14] = a[3]
                  s[15] = a[4]
                  s[17] = a[5]
                  s[18] = a[6]

               elif nnodes == 10:
                  #print('Reordering TET10!')

                  a = s[9]
                  b = s[10]
                  s[9] = b
                  s[10] = a

               elif nnodes == 14:
                  #print('Reordering PYR14!')

                  a = [0] * 5
                  a[0] = s[9]
                  a[1] = s[11]
                  a[2] = s[7]
                  a[3] = s[8]
                  a[4] = s[10]

                  s[7] = a[0]
                  s[8] = a[1]
                  s[9] = a[2]
                  s[10] = a[3]
                  s[11] = a[4]

               stripline = [str(i) for i in s]
               line = ' '.join(stripline)
               g.write(line+'\n')
               line = True
            else:
               g.write(line)
               line = True

f.close()
g.close()
