#!/usr/bin/env python

# Import modules
import numpy as np
from matplotlib.pyplot import *
from math import *

def aieq(xieq, zieq, xiplus1eq, ziplus1eq, Pxeq):
		return (xiplus1eq-Pxeq)+(ziplus1eq*(((xiplus1eq-Pxeq)-(xieq-Pxeq))/(zieq-ziplus1eq)))	
def thetaieq(xieq, zieq, Pxeq):
		return atan2(zieq,(xieq-Pxeq))
def thetaiplus1eq(xiplus1eq, ziplus1eq, Pxeq):
		return atan2(ziplus1eq,(xiplus1eq-Pxeq))
def phiieq(xieq, zieq, xiplus1eq, ziplus1eq, Pxeq):
		return atan2((ziplus1eq-zieq),((xiplus1eq-Pxeq)-(xieq-Pxeq)))

# Setup the subplot for the subsurface bodies
subplot(2,1,2)
maxdistance = 134	#kilometers
maxdepth = 24		#kilometers
Pxarray = np.linspace(-1, maxdistance, maxdistance)
Px = np.array(Pxarray).tolist()
G = 6.67
input = open('LS-18_Density1.txt','r')
readfile = input.read()
allbodies = readfile.split('Body ')[1:]
Zsummed = np.zeros(len(Px))
# Split the bodies at each body ^
for body in allbodies:
    lines = body.split('\n')
    # Splits the lines by each new line
    print (' ')
    print ('line=',lines)
    label = lines[0]
    print ('Body =', label)
    rho = float(lines[1])
    print ('Density =',rho)
    vertices = []
    alphabet = set('abcdefghijklmnopqrstuvwxyz')
    for line in lines[2:]:
        coordinates = line.split('\t')
        if len(coordinates)==2:
            v = (float(coordinates[0]),float(coordinates[1]))
            print (v)
            vertices.append(v)
    #Loops each vertice for all x values
    for i in range(len(vertices)):
        Z = []
        xi = vertices[i][0]
        zi = vertices[i][1]	
        if i < (len(vertices) - 1):
            xiplus1 = vertices[i+1][0]
            ziplus1 = vertices[i+1][1]
        elif i == (len(vertices) - 1):
            xiplus1 = vertices[0][0]
            ziplus1 = vertices[0][1]
        for Pxi in Px:
            if xiplus1 == ziplus1 == 0 or xi == zi:
                thetai = thetaieq(xi, zi, Pxi)
                thetaiplus1 = thetaiplus1eq(xiplus1, ziplus1, Pxi)
                if thetai == thetaiplus1:
                    Zeq = 0
                else:
                    Zeq = 0
            elif xi == xiplus1:
                thetai = thetaieq(xi, zi, Pxi)
                thetaiplus1 = thetaiplus1eq(xiplus1, ziplus1, Pxi)
                Zeq = xi*log(cos(thetai)/cos(thetaiplus1))
            elif zi == ziplus1:
                thetai = thetaieq(xi, zi, Pxi)
                thetaiplus1 = thetaiplus1eq(xiplus1, ziplus1, Pxi)
                Zeq = zi*(thetaiplus1 + thetai)
            elif xiplus1 == 0:
                ai = aieq(xi, zi, xiplus1, ziplus1, Pxi)
                thetai = thetaieq(xi, zi, Pxi)
                phii = phiieq(xi, zi, xiplus1, ziplus1, Pxi)
                Zeq = ai*sin(phii)*cos(phii)*(thetai-(pi/2)+tan(phii)*log(cos(thetai)*(tan(thetai)-tan(phii))))
            elif xi == 0:
                ai = aieq(xi, zi, xiplus1, ziplus1, Pxi)
                thetaiplus1 = thetaiplus1eq(xiplus1, ziplus1, Pxi)
                phii = phiieq(xi, zi, xiplus1, ziplus1, Pxi)
                Zeq = ai*sin(phii)*cos(phii)*(thetaiplus1-(pi/2)+tan(phii)*log(cos(thetaiplus1)*(tan(thetaiplus1)-tan(phii))))
            else:
                ai = aieq(xi, zi, xiplus1, ziplus1, Pxi)
                thetai = thetaieq(xi, zi, Pxi)
                thetaiplus1 = thetaiplus1eq(xiplus1, ziplus1, Pxi)
                phii = phiieq(xi, zi, xiplus1, ziplus1, Pxi)
                Zeq = ai*sin(phii)*cos(phii)*(thetai-thetaiplus1+tan(phii)*log((cos(thetai)*(tan(thetai)-tan(phii)))/(cos(thetaiplus1)*(tan(thetaiplus1)-tan(phii)))))
            Z.append(Zeq)
            Zarray = np.array(Z)
        print ('Zsum=',Zsummed)
        print ('rho=',rho)
        print ('Zarray=',Zarray)
        Zsummed = Zsummed + rho * Zarray
    print (rho)
    print ('Vertices =', vertices) 
    fill([i[0] for i in vertices], [i[1] for i in vertices], 'green')	
V = 2*G*Zsummed
depthparameters = [26,maxdistance,-maxdepth,0]
xlabel('Distance (km)')
ylabel('Depth (km)')
axis(depthparameters)

# Gravity Plot
subplot(2,1,1)
#observedgravity = np.genfromtxt('LS-08observedtest.txt') 
#x = observedgravity[:,0]
#y = observedgravity[:,1]
#plot(x, y, 'k.', label='Observed')
plot(Px, V, label = 'Calculated')
legend(loc='lower right', fontsize='medium')
ylabel('Gravity (mgals)')
xlim(0,100)
#ylim(-100,100)
show()

