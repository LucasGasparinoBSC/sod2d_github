import sys
meshName=sys.argv[1]
geofile='{}.geo.dat'.format(meshName)
coord='{}.coord'.format(meshName)

f=open(geofile,'r')
g=open(coord,'w')


line = True
inCoord = False

boundMax=[float('-Inf'),float('-Inf'),float('-Inf')]
boundMin=[float('Inf' ),float('Inf' ),float('Inf' )]


while line:
	line = f.readline()
	
	if inCoord:
		if 'END' in line:
			if boundMin[2] == float('Inf' ):
				boundMin[2] = 0.0
			if boundMax[2] == float('-Inf'):
				boundMax[2] = 0.0
			print('--| Bounding Box: ({},{},{}),({},{},{})'.format(boundMin[0],boundMin[1],boundMin[2],boundMax[0],boundMax[1],boundMax[2]))
			print('--| End writing coordinates')
			inCoord = False
			line = False
		else:
			data = line.split()
			for i in range(len(data)-1):
				boundMax[i] = max(boundMax[i],float(data[i+1])) 
				boundMin[i] = min(boundMin[i],float(data[i+1])) 
			g.write(line)

	if line:
		if 'COORD' in line:
			inCoord = True
			print('--| Start writing coordinates')


f.close()
g.close()

