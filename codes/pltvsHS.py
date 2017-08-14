#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import glob
from obspy.core import UTCDateTime
import sys

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

template='2'

sens = ['10_OKR1','00_OKR1', '00_OK038']
gooddata =[[],[],[]]

for day in range(118,128):
    temps=[]
    for idx, sen in enumerate(sens):
        f = open('TemplateResults_' + sen + '_' + template + '_' + str(day),'r')
        data = False
        vals=[]
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
        #vals = np.asarray(vals)
        #vals=np.abs(vals)
        temps.append(vals)
    minval = min([len(temps[0]), len(temps[1]), len(temps[2])])-1
    print(minval)
    for idx, temp in enumerate(temps):
        print(len(temp))
        gooddata[idx] += temp[:minval]
print(len(gooddata))

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

g0 = np.abs(np.asarray(gooddata[0]))
g1 = np.abs(np.asarray(gooddata[1]))
g2 = np.abs(np.asarray(gooddata[2]))

p02 = np.polyfit(g0,g2,1)
p12 = np.polyfit(g1,g2,1)
x = np.linspace(0,1,100)
    
plt.figure(1,figsize=(8,8))    
plt.subplots_adjust(hspace=0.001)
ax = plt.subplot(111)
ax.plot(g0,g2,'.', c='C1',label='00 OKR1', markersize=10, alpha=.5)
ax.plot(x, np.poly1d(p02)(x),c='C1', linewidth=4, label='00 OKR1 Slope: ' + str(round(p02[0],1)))

ax.plot(g1,g2,'.', c='C2', label='10 OKR1', markersize=10, alpha=.5)
ax.plot(x, np.poly1d(p12)(x),'--', c='C2', linewidth=3, label= '10 OKR1 Slope: ' + str(round(p02[0],1)))
plt.plot([.7, .7], [-1., 1.], 'k', linewidth=2)
plt.plot([-1., 1.], [.7, .7], 'k', linewidth=2)
plt.ylabel('Absolute Correlation (Translational)')
plt.xlabel('Absolute Correlation (Rotational)')
plt.xlim([0.0, 1.])
plt.ylim([0.0, 1.])

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * .1,
                 box.width, box.height * 0.9])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=False, shadow=False, ncol=2)

plt.savefig('TemplateCorrPlot_' + template + '.jpg',format='JPEG', dpi=400)
#plt.savefig('TemplateCorrPlot_' + template + '.pdf',format='PDF', dpi=400)
#plt.show()

sys.exit()


#for idx, sen in enumerate(sens):
    #files  = glob.glob('TemplateResults_' + sen + '_' + temp + '_*')
    #files =np.sort(files)
    ##print(files)
    #vals =[]
    #for curfile in files:
        #f = open(curfile,'r')
        #data = False
        #for line in f:
            #if 'time' in line:
                #timeStr = UTCDateTime(line.split(':')[1])
            #if 'Magnitude:' in line:
                #mag = line.split(':')
            #if data:
                #vals.append(float(line))
            #if 'data' in line:
                #data = True
        #f.close()   
    #vals = np.asarray(vals)
    #vals=np.abs(vals)
    #goods.append(vals)
    #print(len(vals))
    
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
