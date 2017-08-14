#!/usr/bin/env python
import matplotlib.pyplot as plt
import glob
import numpy as np
from scipy.optimize import fmin

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)



def getVars(var):
    res=[]
    files= glob.glob('Results_E*')
    for curfile in files:
        f = open(curfile,'r')
        if sum(1 for line in f) < 98:
            f.close()
            continue
        else:
            f.close()
            f= open(curfile,'r')
            for line in f:
                if var in line:
                    if var == 'Magnitude':
                        res.append(float(line.split(" ")[1]))
                    else:
                        res.append(float(line.split(":")[1]))
            f.close()
    res=np.asarray(res)
    return res

def getEvent():
    res=[]
    files= glob.glob('Results_E*')
    for curfile in files:
        f = open(curfile,'r')
        if sum(1 for line in f) < 98:
            f.close()
            continue
        else:
            f.close()
            res.append(int(curfile.split('_')[-1]))
    res=np.asarray(res)
    return res

##########################################################################################
#Plot of moment and phase velocity




distance = getVars('Distance')

pT00 = getVars('Peak Vertical 00 Rotation:')
pT10 = getVars('Peak Vertical 10 Rotation:')


mags = getVars('Magnitude')


p00RT = getVars('Peak Radial 00 Translation Rate:')
p00TT = getVars('Peak Transverse 00 Translation Rate:')

PGV = []
for pair in zip(p00RT, p00TT):
    PGV.append(max(pair))
PGV = np.asarray(PGV)

mos = []
for mag in mags:
    mos.append(10**((mag+6.)*1.5))


mos = np.asarray(mos)


 
cs =[]
for pair in zip(PGV, pT00):
    cs.append(pair[0]/(2*pair[1]))

for pair in zip(PGV, pT10):
    cs.append(pair[0]/(2*pair[1]))
    
c=np.mean(cs)



fig =plt.figure(2, figsize=(8,8))
ax =plt.subplot(111)
ax.scatter(np.log10(mos[distance >50.]),PGV[distance >50.]/(2*c),marker='v', label='Scaled PGV $>$ 50 km', s=35, color=(0.,0.,0.0))
ax.scatter(np.log10(mos[distance >50.]),pT00[distance >50.], marker='v', label='PG$\omega$ 00 $>$ 50 km', s=35, color=(230./255., 159./255., 0.))
ax.scatter(np.log10(mos[distance >50.]),pT10[distance >50.], marker='v', label='PG$\omega$ 10 $>$ 50 km', s=35, color=(86./255., 180./255., 33./255.))
ax.scatter(np.log10(mos[distance <=50.]),PGV[distance <=50.]/(2*c),  label='Scaled PGV $\leq$ 50 km', s=35, color=(0., 158./255., 115./255.))
ax.scatter(np.log10(mos[distance <=50.]),pT00[distance <=50.],  label='PG$\omega$ 00 $\leq$ 50 km', s=35, color=(240./255., 228./255., 66./255.))
ax.scatter(np.log10(mos[distance <=50.]),pT10[distance <=50.],  label='PG$\omega$ 10 $\leq$ 50 km', s=35, color=(0., 114./255., 178./255.))
ax.set_yscale("log")
plt.ylim((10**-9,10**-3))
#plt.semilogy(distance,np.abs(newval),'.',label='Theoretical')
plt.xlabel('Seismic Momoment (log(dyne$\cdot$cm))')
plt.ylabel('Peak Ground Rotation (Radians)')
plt.legend(loc='upper left')
plt.title('Apparent c= ' + str(int(c)) + '$\pm$' + str(int(np.std(cs))) + ' (m/s)')
plt.savefig('PGW_PGV.jpg', format='jpeg', dpi=400)
plt.savefig('PGW_PGV.pdf',format='pdf', dpi=400)
plt.clf()
plt.close()





#############################################################################################
#Apparent phase velocity as a function of distance
cFM=[]
cFS=[]
events=[]
dvec = range(10,200)
for dis in dvec:

    PGVT = PGV[distance<=dis]
    pT00T = pT00[distance<=dis]
    pT10T = pT10[distance<=dis]
    events.append(len(PGVT))

    cs =[]
    for pair in zip(PGVT, pT00T):
        cs.append(pair[0]/(2*pair[1]))

    for pair in zip(PGVT, pT10T):
        cs.append(pair[0]/(2*pair[1]))
    
    cFM.append(np.mean(cs))
    cFS.append(np.std(cs))



cFS= np.asarray(cFS)
cFM = np.asarray(cFM)
fig = plt.figure(1)
fig, ax1 = plt.subplots()
ax1.plot(np.asarray(dvec),cFM,linewidth=2,color='k',label='Estimated c')
ax1.set_xlabel('Maximum Event Distance (km)')
ax1.set_ylabel('Estimated c (m/s)')
ax1.set_xlim((10,200))
ax2= ax1.twinx()
ax2.plot(np.asarray(dvec), events, linewidth=2,color='.7', label='Events Used')
ax2.set_ylabel('Events Used',color='.7')
ax2.tick_params('y', colors='.7')
ax2.set_xlim((10,200))
plt.savefig('c_estimate.jpg',format='jpeg',dpi=400)
plt.savefig('c_estimate.pdf',format='pdf',dpi=400)
plt.clf()
plt.close()





#################################################################################################
#HV ratio versus VH ratio


VT = getVars('Peak Spectral Value 00 Translational Vertical:')
HT1 = getVars('Peak Spectral Value 00 Translational Transverse:')
HT2 = getVars('Peak Spectral Value 00 Translational Radial:')
 
 
HT = (HT1 + HT2)/2.


# HVT is our translational HV ratio

VR00 = getVars('Peak Spectral Value 00 Rotational Vertical:')
HR1 = getVars('Peak Spectral Value 00 Rotational Transverse:')
HR2 = getVars('Peak Spectral Value 00 Rotational Radial:')
HR00 = (HR1 + HR2)/2.

VR10 = getVars('Peak Spectral Value 10 Rotational Vertical:')
HR1 = getVars('Peak Spectral Value 10 Rotational Transverse:')
HR2 = getVars('Peak Spectral Value 10 Rotational Radial:')

HR10 = (HR1 + HR2)/2.


fVR00 = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
fVR10 = getVars('Peak Spectral Frequency Value 10 Rotational Vertical:')
fVT = getVars('Peak Spectral Frequency Value 00 Translational Vertical:')

fHR1 = getVars('Peak Spectral Frequency Value 00 Rotational Transverse:')
fHR2 = getVars('Peak Spectral Frequency Value 00 Rotational Radial:')
fHR00 = (fHR1 + fHR2)/2.

fHR1 = getVars('Peak Spectral Frequency Value 10 Rotational Transverse:')
fHR2 = getVars('Peak Spectral Frequency Value 10 Rotational Radial:')
fHR10 = (fHR1 + fHR2)/2.

fH1 = getVars('Frequency Peak 00 Radial:')
fH2 = getVars('Frequency Peak 00 Transverse:')

fHT = (fH1 + fH2)/2.



lf00= np.poly1d(np.polyfit(HR00,VR00,1))
lf10 = np.poly1d(np.polyfit(HR10,VR10,1))
lfT = np.poly1d(np.polyfit(VT,HT,1))


fig = plt.figure(1)
plt.plot(VR00,HR00,'.',color='g')
plt.plot(VR10,HR10,'.',color='b')
plt.plot(VT,HT,'.',color='r')
plt.plot(HR00, lf00(HR00), color='g')
plt.plot(HR10, lf10(HR10), color='b')
plt.plot(lfT(VT),VT, color='r')

#plt.show()
plt.clf()
plt.close()

#######################################################################################################
# Mean spectral ratios
VT =  getVars('Mean Spectral 00 Vertical Translational:')
HT1 =getVars('Mean Spectral 00 Radial Translational:')
HT2 =getVars('Mean Spectral 00 Transverse Translational:')
HT = (HT1 + HT2)/2.

# HVT is our translational HV ratio

VR00 = getVars('Mean Spectral 00 Vertical Rotational:')
HR1 = getVars('Mean Spectral 00 Transverse Rotational:')
HR2 = getVars('Mean Spectral 00 Radial Rotational:')
HR00 = (HR1 + HR2)/2.

VR10 = getVars('Mean Spectral 10 Vertical Rotational:')
HR1 = getVars('Mean Spectral 10 Transverse Rotational:')
HR2 = getVars('Mean Spectral 10 Radial Rotational:')
HR10 = (HR1 + HR2)/2.





lf00= np.polyfit(VR00,HR00,1)
lf10 = np.polyfit(VR10,HR10,1)
lfT = np.polyfit(VT,HT,1)


lfTp = np.poly1d(np.polyfit(VT,HT,1))
lf00p = np.poly1d(np.polyfit(VR00,HR00,1))
lf10p = np.poly1d(np.polyfit(VR10,HR10,1))


print(lfT)
print(lf00)
print(lf10)

fig = plt.figure(1, figsize=(8,8))
plt.plot(VR00,HR00,'.',color='k',label='00 Rotational')
plt.plot(VR10,HR10,'.',color=(230./255., 159./255., 0.), label='10 Rotational')
plt.plot(VT,HT,'.',color=(86./255.,180./255.,233./255.), label='Translational')
plt.plot(VR00, lf00p(VR00), color='k', label = 'Slope: ' + str(round(lf00[0],2)) + ' intercept: ' + str(round(lf00[1],2)))
plt.plot(VR10, lf10p(VR10), color=(230./255., 159./255., 0.), label = 'Slope: ' + str(round(lf10[0],2)) + ' intercept: ' + str(round(lf10[1],2)))
plt.plot(VT,lfTp(VT), color=(86./255.,180./255.,233./255.), label = 'Slope: ' + str(round(lfT[0],2)) + ' intercept: ' + str(round(lfT[1],2)))
plt.xlabel('Mean Spectral Vertical Power (dB)')
plt.ylabel('Mean Spectral Horizontal Power (dB)')
plt.legend()
plt.savefig('HVratio.jpg',format='jpeg',dpi=400)
plt.savefig('HVratio.pdf', format='pdf', dpi=400)
plt.clf()
plt.close()

########################################################################################3
# Grab some statistics  
pgwH00 = getVars('PGwRsH00:')
pgwZ00 = getVars('PGwRsZ00:')
pgwH10 = getVars('PGwRsH10:')
pgwZ10 = getVars('PGwRsZ10:')
pgVH = getVars('PGVH:')
pgVZ = getVars('PGVZ:')
mags = getVars('Magnitude')
events = getEvent()



fig = plt.figure(1, figsize=(8,8))
plt.loglog(pgwH00, pgVZ ,'.', markersize=14.,label='PGV$_Z$ vs. 00 PG$\dot{\omega}_H$', color=(0., 0., 0.))
plt.loglog(pgwH10, pgVZ ,'.', markersize=14., label='PGV$_Z$ vs. 10 PG$\dot{\omega}_H$', color=(230./255., 159./255., 0.))
plt.loglog(pgwZ00, pgVH, '+', markersize=14., label='PGV$_H$ vs. 00 PG$\dot{\omega}_Z$', color=(86./255.,180./255.,233./255.))
plt.loglog(pgwZ10, pgVH, '+', markersize=14., label='PGV$_H$ vs. 10 PG$\dot{\omega}_Z$', color=(0.,158./255.,115./255.))
plt.legend()
plt.xlabel('Peak Ground Velocity (m/s)')
plt.ylabel('Peak Ground Rotation Rate (Radians/s)')
plt.savefig('PGRVSPGV.jpg',format='jpeg')
plt.savefig('PGRVSPGV.pdf', format='pdf')
plt.clf()
plt.close()



event = events[np.argmax(pgwH00)]
print(event)
event = events[np.argmax(pgwH10)]
print(event)
event = events[np.argmax(pgwZ00)]
print(event)
event = events[np.argmax(pgwZ10)]
print(event)
event = events[np.argmax(pgVH)]
print(event)
event = events[np.argmax(pgVZ)]
print(event)
#m1 =max(pR10R)
#m2= max(pR00R)
#m3= max(pT00T)
#m4= max(pT10T)
#print(m1)
#print(m2)
#print(m3)
#print(m4)
#print(max(pV00R))
#print(max(pV10R))

#distance = getVars('Distance:')
#dm = distance[np.argmax(m3)]
#events = getEvent()
#print(dm)
#print(events[np.argmax(m3)])

##################################################################################################3
pgwH00 = getVars('PGwRsH00:')
pgwZ00 = getVars('PGwRsZ00:')
pgwH10 = getVars('PGwRsH10:')
pgwZ10 = getVars('PGwRsZ10:')
pgVH = getVars('PGVH:')
pgVZ = getVars('PGVZ:')
mags = getVars('Magnitude')
distance = getVars('Distance')






pgwH00=pgwH00[distance>10.]*10**9
pgwH10=pgwH10[distance>10.]*10**9
pgwZ10=pgwZ10[distance>10.]*10**9
pgwZ00=pgwZ00[distance>10.]*10**9
pgVH = pgVH[distance>10.]*10**9
pgVZ = pgVZ[distance>10.]*10**9
mags = mags[distance>10.]*10**9
distance = distance[distance>10.]


#pgwH00 = np.divide(pgwH00,10**(1.5*mags+6))*10**9
#pgwH10 = np.divide(pgwH10,10**(1.5*mags+6))*10**9
#pgwZ00 = np.divide(pgwZ00,10**(1.5*mags+6))*10**9
#pgwZ10 = np.divide(pgwZ10,10**(1.5*mags+6))*10**9
#pgVH = np.divide(pgVH,10**(1.5*mags+6))*10**9
#pgVZ = np.divide(pgVZ,10**(1.5*mags+6))*10**9


print(pgwH00)


def fun(x, a,c):
    return a*x+c

from scipy.optimize import curve_fit

popt0, pcov0 = curve_fit(fun, distance, pgwH00)
popt1, pcov1 = curve_fit(fun, distance, pgwH10)
popt2, pcov2 = curve_fit(fun, distance, pgwZ00)
popt3, pcov3 = curve_fit(fun, distance, pgwZ10)
popt4, pcov4 = curve_fit(fun, distance, pgVH)
popt5, pcov5 = curve_fit(fun, distance, pgVZ)


print(popt0[0])



fig = plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
plt.subplot(211)
plt.semilogy(distance, pgwH00, 's',color=(0.,0.,0.), markersize=14.,label='00 PG$\dot{\omega}_H$ Slope=' + str(round(popt0[0],1)))
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt0),'--', color=(0.,0.,0.), linewidth=2.)
plt.semilogy(distance, pgwZ00, 'v', markersize=14., label='00 PG$\dot{\omega}_Z$ Slope=' + str(round(popt1[0],1)), color=(230./255., 159./255., 0.))
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt1),'b--', linewidth=2., color=(230./255.,159./255., 0.))
plt.semilogy(distance, pgwH10,'s', markersize=14.,label='10 PG$\dot{\omega}_H$ Slope=' + str(round(popt2[0],1)), color=(86./255.,180./255., 233./255.))
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt2),'b--', linewidth=2., color=(86./255.,180./255., 233./255.))
plt.semilogy(distance, pgwZ10, 'v', markersize=14., label='10 PG$\dot{\omega}_Z$ Slope=' + str(round(popt3[0],1)), color=(0./255.,158./255., 115./255.))
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt3),'b--', linewidth=2., color=(0./255.,158./255., 115./255.))
plt.legend(loc=2)
plt.xlim([0.,200.])
plt.ylabel('Peak Ground Rotation Rate (nRadians/s)')
plt.xticks([])
plt.subplot(212)
plt.semilogy(distance, pgVH, '.', markersize=18., label='PGV$_H$ Slope=' + str(round(popt4[0],1)), color=(240./255.,240./255., 66./255.))
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt4),'--', color=(240./255.,240./255., 66./255.), linewidth=2.)
plt.semilogy(distance, pgVZ, '.', markersize=18., label='PGV$_Z$ Slope=' + str(round(popt5[0],1)), color=(0./255.,114./255., 178./255.) )
plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt5),'--', color=(0./255.,114./255., 178./255.), linewidth=2.)
plt.legend(loc=2)
plt.xlim([0.,200])
plt.xlabel('Distance from Event (km)')
plt.ylabel('Peak Ground Velocity (nm/s)')
#plt.show()
plt.savefig('DistancePGVals.jpg',format='JPEG',dpi=400)
plt.savefig('DistancePGVals.pdf', format='PDF', dpi=400)
plt.clf()
plt.close()





'''

###################################################################################################################3
##################################################################################################3
pr00Z = getVars('Peak Spectral Value 00 Rotational Vertical:')
fr00Z = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
pr00R = getVars('Peak Spectral Value 00 Rotational Radial:')
fr00R = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
pr00T = getVars('Peak Spectral Value 00 Rotational Transverse:')
fr00T = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
pr10Z = getVars('Peak Spectral Value 10 Rotational Vertical:')
fr10Z = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
pr10R = getVars('Peak Spectral Value 10 Rotational Radial:')
fr10R = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')
pr10T = getVars('Peak Spectral Value 10 Rotational Transverse:')
fr10T = getVars('Peak Spectral Frequency Value 00 Rotational Vertical:')




pt00Z = getVars('Peak Spectral Value 00 Translational Vertical:')
ft00Z = getVars('Peak Spectral Frequency Value 00 Translational Vertical:')
pt00R = getVars('Peak Spectral Value 00 Translational Radial:')
ft00R = getVars('Peak Spectral Frequency Value 00 Translational Vertical:')
pt00T = getVars('Peak Spectral Value 00 Translational Transverse:')
ft00T = getVars('Peak Spectral Frequency Value 00 Translational Vertical:')
distance= getVars('Distance')

mags = getVars('Magnitude')

#fsv=fsv[distance>20]
#pvv=pvv[distance>20]
#mags = mags[distance>20]
#distance = distance[distance>20]

#pr00Z = pr00Z[distance>20]
#pr10Z = pr10Z[distance>20]
#pr00R = pr00R[distance>20]
#pr10R = pr10R[distance>20]
#pr00T = pr00T[distance>20]
#pr10T = pr10T[distance>20]
#pt00Z = pt00Z[distance>20]
#pt00R = pt00R[distance>20]
#pt00T = pt00T[distance>20]
#fr00Z = fr00Z[distance>20]
#fr10Z = fr10Z[distance>20]
#fr00R = fr00R[distance>20]
#fr10R = fr10R[distance>20]
#fr00T = fr00T[distance>20]
#fr10T = fr10T[distance>20]

#ft00Z = ft00Z[distance>20]
#ft00R = ft00R[distance>20]
#ft00T=ft00T[distance>20]
#mags = mags[distance>20]
#distance = distance[distance>20]



from scipy.linalg import lstsq

def grabparameter(distance,powv,fs,mags):
    #powv = np.multiply(np.multiply(10**(powv/20.),np.sqrt(fs)), (2*np.pi*fs)**2)

    A = np.zeros((len(powv),3))
    B = np.zeros((len(powv),1))

    for idx, val in enumerate(zip(distance,powv,fs,mags)):
        dis = val[0]
        p = val[1]
        f=val[2]
        mag=val[3]
        A[idx,0]=-np.pi*f*dis/(3.5) 
        A[idx,1]=10**(1.5*mag+6.)
        A[idx,2]=-np.log(dis)
        B[idx,0]=p/np.log10(np.e)

    c,_,_,_ = lstsq(A,B)
    Q = c[0,0]
    b = c[1,0]
    g = c[2,0] 
    return Q,b,g
    
    
def newfun(Q,b,g,distance):
    #vals = np.log(distance)**(-g) 
    vals = np.log(distance)**(-g) + b*10**(1.5*2+6) -Q*np.pi*distance*5/(3.5)   
    #vals = 20*np.log10(vals
    return vals
    
Q00Z,b00Z,g00Z = grabparameter(distance, pr00Z, fr00Z,mags)
print(str(Q00Z) + ' ' + str(b00Z) + ' ' + str(g00Z))
Q,b,g = grabparameter(distance, pr10Z, fr10Z,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pr00R, fr00R,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pr10R, fr10R,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pr00T, fr00T,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pr10T, fr10T,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pt00Z, ft00Z,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pt00R, ft00R,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
Q,b,g = grabparameter(distance, pt00T, ft00T,mags)
print(str(Q) + ' ' + str(b) + ' ' + str(g))
fig = plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
plt.subplot(211)

blah = 10**(np.sqrt(fr00Z)*pr00Z/20.)/distance
blah= blah*(10**(1.5*mags+6))
disv = np.linspace(10,200,200)
#pr00Z = np.multiply(np.multiply(10**(pr00Z/20.),np.sqrt(fr00Z)), (2*np.pi*fr00Z)**2)
plt.plot(distance, pr00Z -pt00R, '.', markersize=14., color='C1', label='P')
#plt.plot(distance, pr00R -pt00R, '.', markersize=14., color='C2', label='P')
#plt.plot(distance, pr00T -pt00T, '.', markersize=14., color='C3', label='P')

plt.show()
import sys
sys.exit()
plt.plot(distance, pr10Z, '.', markersize=14., color='C2', label='P')
plt.plot(distance, pr00R, 'v', markersize=14., color='C3', label='P')
plt.plot(distance, pr10R, 'v', markersize=14., color='C1', label='P')
plt.plot(distance, pr00T, 's', markersize=14., color='C2',label='P')
plt.plot(distance, pr10T, 's', markersize=14., color='C3', label='P')
plt.subplot(212)
plt.plot(distance, pt00Z, '.', markersize=14., color='C1', label='P')
plt.plot(distance, pt00R, 'v', markersize=14., color='C2', label='P')
plt.plot(distance, pt00T, 's', markersize=14., color='C3', label='P')
plt.show()


#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt0),'--', color='C2', linewidth=2.)
#plt.semilogy(distance, pgwZ00, 'v', markersize=14., label='00 PG$\dot{\omega}_Z$ Slope=' + str(round(popt1[0]*1000,1)), color='C3')
#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt1),'b--', linewidth=2., color='C3')
#plt.semilogy(distance, pgwH10,'v', markersize=14.,label='10 PG$\dot{\omega}_H$ Slope=' + str(round(popt2[0]*1000,1)), color='C4')
#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt2),'b--', linewidth=2., color='C4')
#plt.semilogy(distance, pgwZ10, 'v', markersize=14., label='10 PG$\dot{\omega}_Z$ Slope=' + str(round(popt3[0]*1000,1)), color='C5')
#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt3),'b--', linewidth=2., color='C5')
#plt.legend(loc=2)
#plt.xlim([10.,200.])
#plt.ylabel('Corrected Peak Ground Rotation Rate (nRadians/s)')
#plt.xticks([])
#plt.subplot(212)
#plt.semilogy(distance, pgVH, '.', markersize=18., label='PGV$_H$ Slope=' + str(round(popt4[0]*1000,1)))
#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt4),'--', color='C0', linewidth=2.)
#plt.semilogy(distance, pgVZ, '.', markersize=18., label='PGV$_Z$ Slope=' + str(round(popt5[0]*1000,1)) )
#plt.semilogy(np.linspace(0,200,200), fun(np.linspace(0,200,200),*popt5),'--', color='C1', linewidth=2.)
#plt.legend(loc=2)
#plt.xlim([10.,200])
#plt.xlabel('Distance from Event (km)')

#pgwH00=pgwH00[distance>10.]
#pgwH10=pgwH10[distance>10.]
#pgwZ10=pgwZ10[distance>10.]
#pgwZ00=pgwZ00[distance>10.]
#pgVH = pgVH[distance>10.]
#pgVZ = pgVZ[distance>10.]
#mags = mags[distance>10.]
#distance = distance[distance>10.]
'''

##########################################################################################
#Plot of moment and phase velocity




distance = getVars('Distance')

pT00 = getVars('Peak Transverse 00 Rotation:')
pT10 = getVars('Peak Transverse 10 Rotation:')


mags = getVars('Magnitude')


PGV = getVars('Peak Vertical 00 Translation Rate:')




mos = []
for mag in mags:
    mos.append(10**((mag+6.)*1.5))


mos = np.asarray(mos)


 
cs =[]
for pair in zip(PGV, pT00):
    cs.append(pair[0]/(2*pair[1]))

for pair in zip(PGV, pT10):
    cs.append(pair[0]/(2*pair[1]))
    
c=np.mean(cs)



fig =plt.figure(2, figsize=(8,8))
ax =plt.subplot(111)
ax.scatter(np.log10(mos[distance >50.]),PGV[distance >50.]/(2*c),marker='v', label='Scaled PGV $>$ 50 km', s=35, color=(0.,0.,0.0))
ax.scatter(np.log10(mos[distance >50.]),pT00[distance >50.], marker='v', label='PG$\omega$ 00 $>$ 50 km', s=35, color=(230./255., 159./255., 0.))
ax.scatter(np.log10(mos[distance >50.]),pT10[distance >50.], marker='v', label='PG$\omega$ 10 $>$ 50 km', s=35, color=(86./255., 180./255., 33./255.))
ax.scatter(np.log10(mos[distance <=50.]),PGV[distance <=50.]/(2*c),  label='Scaled PGV $\leq$ 50 km', s=35, color=(0., 158./255., 115./255.))
ax.scatter(np.log10(mos[distance <=50.]),pT00[distance <=50.],  label='PG$\omega$ 00 $\leq$ 50 km', s=35, color=(240./255., 228./255., 66./255.))
ax.scatter(np.log10(mos[distance <=50.]),pT10[distance <=50.],  label='PG$\omega$ 10 $\leq$ 50 km', s=35, color=(0., 114./255., 178./255.))
ax.set_yscale("log")
plt.ylim((10**-9,10**-3))
#plt.semilogy(distance,np.abs(newval),'.',label='Theoretical')
plt.xlabel('Seismic Momoment (log(dyne$\cdot$cm))')
plt.ylabel('Peak Ground Rotation (Radians)')
plt.legend(loc='upper left')
plt.title('Apparent c= ' + str(int(c)) + '$\pm$' + str(int(np.std(cs))) + ' (m/s)')
plt.savefig('Rayleigh_PGW_PGV.jpg', format='jpeg', dpi=400)
plt.savefig('Rayleigh_PGW_PGV.pdf',format='pdf', dpi=400)
plt.clf()
plt.close()





