#!/usr/bin/env python
from obspy.core import UTCDateTime, read, Stream
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal import PPSD
from matplotlib.mlab import csd

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib as mpl
import math
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)
from scipy.optimize import fmin

def computeresp(resp,delta,lenfft):
    respval, freq = paz_to_freq_resp(np.asarray(resp['poles']),np.asarray(resp['zeros']),1.,t_samp = delta, 
        nfft=lenfft,freq = True)

    print(resp['sensitivity'])
    respval *= resp['sensitivity']
    respval = np.absolute(respval*np.conjugate(respval))

    respval = respval[1:]
    freq = freq[1:]
    return respval, freq



# Here are the poles/zeros for each of the ATA sensors
paz1202X = {'gain': 1., 'zeros': [0.,0., -9.07085], 
            'poles': [-4.31551, -4.3363, -4.2547, -2.05994*10**4, 
            -2.24048*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.50309*10**18}
paz1202Y = {'gain': 1., 'zeros': [0., 0., -0.635013], 
            'poles': [-2.1234, -2.131, -2.08896, -3.01964*10**4, 
            -3.29513*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*3.124641*10**18}
paz1202Z = {'gain': 1., 'zeros': [0., 0., -8.51744], 
            'poles': [-4.08456, -4.07133, -3.989, -9.47394*10**3, 
            -1.16777*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*3.63243*10**17}
            
paz1274X = {'gain': 1., 'zeros': [0., 0., -7.71072], 
            'poles': [-4.16958, -4.15483, -4.07121, -6.40679*10**4, 
            -6.46124*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.37317*10**19}
paz1274Y = {'gain': 1., 'zeros': [0., 0., -8.52596], 
            'poles': [-4.42831, -4.41308, -4.3247, -2.24554*10**4, 
            -2.26434*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.58848*10**18}
paz1274Z = {'gain': 1., 'zeros': [0., 0., -9.79158], 
            'poles': [-4.60581, -4.63541, -4.54258, -2.70303*10**4, 
            -2.47819*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.197*10**18}

paz1275X = {'gain': 1., 'zeros': [0., 0., -12.1646], 
            'poles': [-6.72298, -6.76686, -6.5603, -2.22233*10**4, 
            -2.21882*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.66805*10**18}
paz1275Y = {'gain': 1., 'zeros': [0., 0., -12.1404], 
            'poles': [-5.38773, -5.3692, -5.26061, -2.75699*10**4, 
            -2.50238*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.29388*10**18}
paz1275Z = {'gain': 1., 'zeros': [0., 0., -98.6279], 
            'poles': [-7.30373, -4.68949, -1.01088*10**2, -1.52476*10**4, 
            -2.06738*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.06106*10**18}




fig = plt.figure(1,figsize=(8,10))
plt.subplots_adjust(hspace=0.001)
# Do PDF calculation
if True:
    #sens =[]
    #sens.append({'sensor': 'GS_OKR1', 'X': paz1275X, 'Y': paz1275Y, 'Z': paz1275Z})
    ##sens.append({'sensor': 'XX_TST5', 'X': paz1202X, 'Y': paz1202Y, 'Z': paz1202Z})
    #sens.append({'sensor': 'GS_OKR1', 'X': paz1274X, 'Y': paz1274Y, 'Z': paz1274Z})
    sens =[]
    sens.append({'sensor': 'XX_TST4', 'X': paz1275X, 'Y': paz1275Y, 'Z': paz1275Z})
    #sens.append({'sensor': 'XX_TST5', 'X': paz1202X, 'Y': paz1202Y, 'Z': paz1202Z})
    sens.append({'sensor': 'XX_TST6', 'X': paz1274X, 'Y': paz1274Y, 'Z': paz1274Z})
    #stime = UTCDateTime('2017-075T00:00:00.0')
    #etime = UTCDateTime('2017-075T13:30:00.0')
    #etime = stime + 12.*60.*60.
    stime = UTCDateTime('2017-075T00:00:00.0')
    #etime = UTCDateTime('2017-075T13:30:00.0')
    etime = stime + 12.*60.*60.


    for idx, sen in enumerate(sens):
        #print(sen)
        #datastr='/msd/' + sen['sensor'] +'/2017/12*/' + str(int(idx)) +'0_EJ*.seed'
        st = read('/msd/' + sen['sensor'] +'/2017/075/00_EH*.seed')

        
        #st = read(datastr)
        print(st)
        st.trim(stime,etime)
        newmean =[]
        for tr in st:
            if tr.stats.channel == 'EH1':
                comp = 'X'
            elif tr.stats.channel == 'EH2':
                comp ='Y'
            elif tr.stats.channel == 'EH0':
                comp = 'Z'
            #print(tr)
            
            #inst1resp, freq = computeresp(sen[comp],tr.stats.delta,2**12)
            #fig = plt.figure(1)
            #plt.semilogx(freq, 10.*np.log10(inst1resp))
            #plt.show()
            
            h,f = paz_to_freq_resp(sen[comp]['poles'], sen[comp]['zeros'], sen[comp]['sensitivity'], 1./200., 180000, freq=True)
            
            
            sen[comp]['sensitivity'] *= 1./np.abs(h[np.abs(f-10.).argmin()])
            sen[comp]['sensitivity'] *= 2000.

            h,f = paz_to_freq_resp(sen[comp]['poles'], sen[comp]['zeros'], sen[comp]['sensitivity'], 1./200., 180000, freq=True)
            
            nlab = tr.id
            nlab = nlab.replace('H','J')
            nlab = nlab.replace('XX','GS')
            nlab = nlab.replace('TST4.00','OKR1.00')
            nlab = nlab.replace('TST6.00','OKR1.10')
            nlab = nlab.replace('EJ0','EJZ')
            nlab = nlab.replace('EJ1','EJN')
            nlab = nlab.replace('EJ2','EJE')
            
            plt.subplot(211)
            plt.title('Frequency Response of ATA Proto-SMHD Sensors')
            plt.semilogx(f,20.*np.log10(np.abs(h)), label=nlab)
            plt.xlim((.01, 100))
            plt.xticks([])
            plt.ylabel('Amplitude Response (V/radians/s)')
            plt.subplot(212)
            plt.semilogx(f,np.unwrap(np.angle(h))*180./math.pi, label=nlab)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Phase Response (degrees)')
            plt.xlim((.01, 100))
            plt.ylim((0.,180.))
            #plt.legend()
    #stimestr = str(stime.year) + '-' + str(stime.month).zfill(2) + '-' + str(stime.day).zfill(2)
    
    plt.xlabel('Frequency (Hz)')
    #plt.ylabel('ORD (dB rel. 1 $(radian/s)^2$)')
    plt.xlim((.01, 100))
#plt.tight_layout()
plt.subplot(211)
plt.legend(loc=9, ncol=2, bbox_to_anchor=(0.5,0.3), fontsize=14)
plt.savefig('ATARESP.png',format='png', dpi=400)
plt.savefig('ATARESP.pdf',format='pdf',dpi=400)

#plt.show()


    #fig = plt.figure(2,figsize=(12,12))
    #newstd = np.std(newmean)
    #plt.semilogx(1./per,np.mean(newmean))
    #plt.show()

