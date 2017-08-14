#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

perR=[]
perL=[]
R=[]
L=[]

f=open('ok038dsp.dat')

for line in f:
    if "#" in line:
        continue
    else:
        line = ' '.join(line.split())
        line = line.split(' ')
    if 'R' in line[1]:
        R.append(float(line[6]))
        perR.append(float(line[5]))
    else:
        L.append(float(line[6]))
        perL.append(float(line[5]))

f.close()

perL=np.asarray(perL)
perR=np.asarray(perR)
L=np.asarray(L)*1000.
R=np.asarray(R)*1000.


print(np.mean(L[1./perL <=10.]))
print(np.mean(R[1./perR <=10.]))

pL=np.poly1d(np.polyfit(perL,L,1))
pR=np.poly1d(np.polyfit(perR,R,1))


import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


#pers =np.linspace(.2, 5., 100)

fig = plt.figure(1,figsize=(8,8))
plt.plot(1/perL, L,'.', markersize=18, label='Love')
#plt.plot(1/perL,pL(perL))
plt.plot(1/perR, R,'.',markersize=18, label='Rayleigh')
#plt.plot(pers,pR(pers))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase Velocity (m/s)')

plt.legend()
plt.savefig('PhaseVelocity.jpg',format='JPEG',dpi=400)
plt.savefig('PhaseVelocity.pdf',format='PDF',dpi=400)

        
