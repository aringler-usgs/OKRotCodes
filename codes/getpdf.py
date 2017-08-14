#!/usr/bin/env python
from obspy.core import UTCDateTime, read, Stream
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal import PPSD
from matplotlib.mlab import csd

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib as mpl
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







# Do PDF calculation
if True:
    sens =[]
    sens.append({'sensor': 'GS_OKR1', 'X': paz1275X, 'Y': paz1275Y, 'Z': paz1275Z})
    #sens.append({'sensor': 'XX_TST5', 'X': paz1202X, 'Y': paz1202Y, 'Z': paz1202Z})
    sens.append({'sensor': 'GS_OKR1', 'X': paz1274X, 'Y': paz1274Y, 'Z': paz1274Z})

    #stime = UTCDateTime('2017-075T00:00:00.0')
    #etime = UTCDateTime('2017-075T13:30:00.0')
    #etime = stime + 12.*60.*60.

    for idx, sen in enumerate(sens):
        print(sen)
        datastr='/tr1/telemetry_days/' + sen['sensor'] +'/2017/2017_12*/' + str(int(idx)) +'0_EJ*.seed'
        print(datastr)
        
        st = read(datastr)
        print(st)
        #st.trim(stime,etime)
        newmean =[]
        for tr in st:
            if tr.stats.channel == 'EJ1':
                comp = 'X'
            elif tr.stats.channel == 'EJ2':
                comp ='Y'
            elif tr.stats.channel == 'EJZ':
                comp = 'Z'
            #print(tr)
            
            #inst1resp, freq = computeresp(sen[comp],tr.stats.delta,2**12)
            #fig = plt.figure(1)
            #plt.semilogx(freq, 10.*np.log10(inst1resp))
            #plt.show()
            
            h,f = paz_to_freq_resp(sen[comp]['poles'], sen[comp]['zeros'], sen[comp]['sensitivity'], 1./200., 1800, freq=True)
            
            
            sen[comp]['sensitivity'] *= 1./np.abs(h[np.abs(f-10.).argmin()])
            sen[comp]['sensitivity'] *= 2000.*(2**26/40.)

            
            
            
            ppsd = PPSD(tr.stats,sen[comp], db_bins=(-180.,180,1), ppsd_length=1800.0)
            ppsd.add(tr)
            #ppsd.plot(show=True, show_noise_models=False, show_mean=True)
            per, mean = ppsd.get_mean()
            
            mean = 10.**(mean/20.)
            mean *= (2*np.pi*1./per)**-1
            mean = 20.*np.log10(mean)
            newmean.append(mean)
            fig = plt.figure(1,figsize=(12,12))
            plt.semilogx(1./per, mean, label=tr.id)
        #print(st)
        #st.decimate(2)
        #st.detrend()
        #for idx, tr in enumerate(st):
            #t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)/(60.*60.)
            #fig =plt.figure(2,figsize=(12,12))
            #plt.subplot(3,1, idx+1)
            #plt.plot(t, tr.data/(2000.*(2**26/40)), label=tr.id, color='k')
            #plt.xlim((min(t),max(t)))
            #if idx ==1:
                #plt.ylabel('Vel. (Rad/s)')
            #if idx == 2:
                #plt.xlabel('Time (Hours)')
            #plt.legend()
        #plt.tight_layout()
        #plt.show()
        #plt.savefig('Tseries' + st[0].stats.station + '.png',format='png')
        #plt.clf()
    stimestr = str(stime.year) + '-' + str(stime.month).zfill(2) + '-' + str(stime.day).zfill(3)
    plt.title('Mean ATA Noise Estimates ' + stimestr + ' duration: ' + str(int(tr.stats.endtime - stime)) + ' s')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD (dB rel. 1 (rad/s)^2/Hz)')
    plt.xlim((.01, 100))
    plt.legend()
    plt.savefig('ATANoise.png',format='png')

    fig = plt.figure(2,figsize=(12,12))
    newstd = np.std(newmean)
    plt.semilogx(1./per,np.mean(newmean))
    plt.show()

