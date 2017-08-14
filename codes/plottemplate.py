#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import glob
from obspy.core import UTCDateTime

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

temp='2'

sens = ['10_OKR1','00_OKR1', '00_OK038']
#sens=['00_OKR1','00_OK038']

goods=[]

for idx, sen in enumerate(sens):
    files  = glob.glob('TemplateResults_' + sen + '_' + temp + '_*')
    files =np.sort(files)
    #print(files)
    vals =[]
    for curfile in files:
        f = open(curfile,'r')
        data = False
        for line in f:
            if 'time' in line:
                timeStr = UTCDateTime(line.split(':')[1])
            if 'Magnitude:' in line:
                mag = line.split(':')
            if data:
                vals.append(float(line))
            if 'data' in line:
                data = True
        f.close()   
    vals = np.asarray(vals)
    vals=np.abs(vals)
    goods.append(vals)
    print(len(vals))
    
neweve1 = goods[0][goods[0] >.7]
neweve2 = goods[1][goods[1] >.7]
neweve3 = goods[2][goods[2] >.7]

print('High corr events: ' + str(len(neweve1)))
print('High corr events: ' + str(len(neweve2)))
print('High corr events: ' + str(len(neweve3)))


import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)
    
fig = plt.figure(1,figsize=(8,8))    
plt.subplots_adjust(hspace=0.001)
t2 = np.arange(len(goods[2]))/(2.*24*60*60.)

plt.subplot(2,1,1)
plt.plot(t2, goods[2], c='C1')
plt.plot(t2[goods[2] >.7], goods[2][goods[2]>.7], '.', c='C1', markersize=14)
t0 = np.arange(len(goods[0]))/(2.*24*60*60.)
plt.yticks([0., .2, .4, .6, .8])
plt.xticks([])
plt.plot(t0,goods[0], c='C2')
plt.plot(t0[goods[0] >.7], goods[0][goods[0]>.7], '.', c='C2', markersize=14)
plt.text(10.,.9,'(A)', fontsize=22)

plt.xlim((min(t0), max(t0)))
plt.subplot(2,1,2)
plt.plot(t2, goods[2], c='C1')
plt.plot(t2[goods[2] >.7], goods[2][goods[2]>.7], '.',c='C1', markersize=14)
t1 = np.arange(len(goods[1]))/(2.*24*60*60.)
plt.plot(t1,goods[1], c='C0')
plt.plot(t1[goods[1] >.7], goods[1][goods[1]>.7], '.',c='C0', markersize=14)
plt.yticks([0., .2, .4, .6, .8])
plt.xlim((min(t2), max(t2)))
plt.xlabel('Time (days)')
plt.text(10.,.9,'(B)', fontsize=22)

fig.text(0.05, 0.5,'Absolute Correlation', ha='center', va='center', rotation='vertical')   
plt.savefig('TemplateFig_' + temp + '.jpg',format='JPEG', dpi=400)
plt.savefig('TemplateFig_' + temp + '.pdf',format='PDF', dpi=400)

