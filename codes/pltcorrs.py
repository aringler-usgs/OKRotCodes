#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def func(x, a, b, c):
    return a*np.exp(-b*x) + c










mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)

f = open('corrs','r')

OKR100=[]
OKR100D=[]
OKR110=[]
OKR110D=[]
OK38=[]
OK38D=[]


for line in f:
    if ('OKR1, 00' in line) and ('6OKR1' not in line):
        OKR100D.append(float(line.split(',')[4]))
        OKR100.append(float(line.split(',')[5]))
    elif ('OKR1, 10' in line) and ('6OKR1' not in line):
        OKR110D.append(float(line.split(',')[4]))
        OKR110.append(float(line.split(',')[5]))
    elif 'OK038' in line:
        OK38D.append(float(line.split(',')[4]))
        OK38.append(float(line.split(',')[5]))


OKR100=np.abs(np.asarray(OKR100))
OKR100D = np.abs(np.asarray(OKR100D))
OKR110D = np.abs(np.asarray(OKR110D))
OKR110 = np.abs(np.asarray(OKR110))
OK38= np.abs(np.asarray(OK38))
OK38D = np.abs(np.asarray(OK38D))


popt100D, pcov = curve_fit(func, OKR100D, OKR100)
popt110D, pcov = curve_fit(func, OKR110D, OKR110)
popt3800D, pcov = curve_fit(func, OK38D, OK38)
print(np.sort(OKR100D))
vals = np.linspace(0,20,200)
fig = plt.figure(1, figsize=(8,8))
plt.plot(OKR100D,OKR100,'.', label='OKR1 00',color='k', alpha=.3)
plt.plot(vals, func(vals, popt100D[0], popt100D[1], popt100D[2]), color='k', linewidth=3, label='OKR1 00 Decay: ' + str(round(popt100D[1],2)))
plt.plot(OKR110D,OKR110,'.', label='OKR1 10', color=(230./255.,159./255.,0./255.), alpha=.3)
plt.plot(vals, func(vals, popt110D[0], popt110D[1], popt110D[2]), color=(230./255.,159./255.,0./255.), linewidth=3, label='OKR1 10 Decay: ' + str(round(popt110D[1],2)))
plt.plot(OK38D, OK38, '.', label='OK038 00', color=(86./255.,180./255.,233./255.), alpha=.3)
plt.plot(vals, func(vals, popt3800D[0], popt3800D[1], popt3800D[2]), color=(86./255.,180./255.,233./255.), linewidth=3, label='OK038 00 Decay: ' + str(round(popt3800D[1],2)))
plt.xlabel('Relative Event Distance (km)')
plt.ylabel('Absolute Correlation')
plt.xlim((0.,20.))

plt.legend()
plt.savefig('CorrsvsDist.jpg', format='jpeg', dpi=400)
plt.savefig('CorrsvsDist.pdf', format='pdf', dpi=400)
#plt.show()
