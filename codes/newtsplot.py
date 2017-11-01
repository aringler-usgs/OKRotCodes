#!/usr/bin/env python
from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
from obspy.signal.invsim import paz_to_freq_resp
from matplotlib.mlab import csd
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from obspy.imaging.maps import plot_basemap
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
##########################################
# To do




# Pick out peak frequency in spectra
# Difference from horizontal to vertical at peak frequency
# Pick frequency of biggest difference between Horizontal and Vertical ratios

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=12)


phs1=[]
phs2=[]

debug = True
length=4096

def cp(tr1,tr2,lenfft,lenol,delta):
    sr = 1/delta
    cpval,fre = csd(tr1.data,tr2.data,NFFT=lenfft,Fs=sr,noverlap=lenol,scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:].real
    return cpval, fre

#OKR1 10
paz10 = [{'gain': 1., 'zeros': [0., 0., -7.71072], 
            'poles': [-4.16958, -4.15483, -4.07121, -6.40679*10**4, 
            -6.46124*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.37317*10**19},
        {'gain': 1., 'zeros': [0., 0., -8.52596], 
            'poles': [-4.42831, -4.41308, -4.3247, -2.24554*10**4, 
            -2.26434*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.58848*10**18},
        {'gain': 1., 'zeros': [0., 0., -9.79158], 
            'poles': [-4.60581, -4.63541, -4.54258, -2.70303*10**4, 
            -2.47819*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.197*10**18}]

#OKR1 00
paz00 = [{'gain': 1., 'zeros': [0., 0., -12.1646], 
            'poles': [-6.72298, -6.76686, -6.5603, -2.22233*10**4, 
            -2.21882*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.66805*10**18},
    {'gain': 1., 'zeros': [0., 0., -12.1404], 
            'poles': [-5.38773, -5.3692, -5.26061, -2.75699*10**4, 
            -2.50238*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.29388*10**18},
    {'gain': 1., 'zeros': [0., 0., -98.6279], 
            'poles': [-7.30373, -4.68949, -1.01088*10**2, -1.52476*10**4, 
            -2.06738*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.06106*10**18}]

#OKR1 10
paz10A = [{'gain': 1., 'zeros': [0., -7.71072], 
            'poles': [-4.16958, -4.15483, -4.07121, -6.40679*10**4, 
            -6.46124*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.37317*10**19},
        {'gain': 1., 'zeros': [0., -8.52596], 
            'poles': [-4.42831, -4.41308, -4.3247, -2.24554*10**4, 
            -2.26434*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.58848*10**18},
        {'gain': 1., 'zeros': [0., -9.79158], 
            'poles': [-4.60581, -4.63541, -4.54258, -2.70303*10**4, 
            -2.47819*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.197*10**18}]

#OKR1 00
paz00A = [{'gain': 1., 'zeros': [0.,  -12.1646], 
            'poles': [-6.72298, -6.76686, -6.5603, -2.22233*10**4, 
            -2.21882*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.66805*10**18},
    {'gain': 1., 'zeros': [0.,  -12.1404], 
            'poles': [-5.38773, -5.3692, -5.26061, -2.75699*10**4, 
            -2.50238*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.29388*10**18},
    {'gain': 1., 'zeros': [0., -98.6279], 
            'poles': [-7.30373, -4.68949, -1.01088*10**2, -1.52476*10**4, 
            -2.06738*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.06106*10**18}]
#OKR1 10
paz10D = [{'gain': 1., 'zeros': [0., 0.,0., -7.71072], 
            'poles': [-4.16958, -4.15483, -4.07121, -6.40679*10**4, 
            -6.46124*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.37317*10**19},
        {'gain': 1., 'zeros': [0., 0.,0., -8.52596], 
            'poles': [-4.42831, -4.41308, -4.3247, -2.24554*10**4, 
            -2.26434*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.58848*10**18},
        {'gain': 1., 'zeros': [0., 0.,0., -9.79158], 
            'poles': [-4.60581, -4.63541, -4.54258, -2.70303*10**4, 
            -2.47819*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.197*10**18}]

#OKR1 00
paz00D = [{'gain': 1., 'zeros': [0., 0., 0., -12.1646], 
            'poles': [-6.72298, -6.76686, -6.5603, -2.22233*10**4, 
            -2.21882*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.66805*10**18},
    {'gain': 1., 'zeros': [0., 0., 0., -12.1404], 
            'poles': [-5.38773, -5.3692, -5.26061, -2.75699*10**4, 
            -2.50238*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*2.29388*10**18},
    {'gain': 1., 'zeros': [0., 0., 0., -98.6279], 
            'poles': [-7.30373, -4.68949, -1.01088*10**2, -1.52476*10**4, 
            -2.06738*10**4, -1.59155*10**6], 
            'sensitivity': (2**24/40)*1.06106*10**18}]
pazTC = {'zeros': [0.j, 0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],
    'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 754.3*6.29327*10**5}

pazTCD = {'zeros': [0.j, 0.j, 0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],
    'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 754.3*6.29327*10**5}
        
pazTCA = {'zeros': [0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],
    'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 754.3*6.29327*10**5}

stime = UTCDateTime('2017-117T00:00:00.0')
etime = UTCDateTime('2017-158T00:00:00.0')
client = Client("IRIS")

stalat = 36.47819
stalon = -98.742240


def writeEventInfo(event,fevent, debug=False):
    if debug:
        print(event)
    eventtime = event.origins[0].time
    eventtimeStr = str(eventtime.year) + ' ' + str(eventtime.julday).zfill(3) + ' ' + \
                    str(eventtime.hour).zfill(2) + ':' + str(eventtime.minute).zfill(2) + ':' + \
                    str(eventtime.second).zfill(2)
    mag = event.magnitudes[0].mag
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude,event.origins[0].longitude)
    dis *= 1./1000.
    fevent.write('Event time: ' + eventtimeStr + '\n')
    fevent.write('Magnitude: ' + str(mag) + ' ' + event.magnitudes[0].magnitude_type + '\n')
    fevent.write('Depth: ' + str(event.origins[0].depth) + '\n')
    fevent.write('Latitude:' + str(event.origins[0].latitude)+ '\n')
    fevent.write('Longitude:' + str(event.origins[0].longitude) + '\n')
    if debug:
        print('Azimuth1: ' + str(azi))
        print('Back azimuth: ' + str(bazi))
        print('distance: ' + str(dis))
    fevent.write('Azimuth: ' + str(azi) + '\n')
    fevent.write('Back azimuth: ' + str(bazi) + '\n')
    fevent.write('Distance: ' + str(dis) + '\n')
    
    return bazi

def grabdata(eventtime, loc, sta, bazi, fevent):
    if sta == 'OKR1':
        chan = 'EJ'
    else:
        chan = 'HH'
    path = '/msd/GS_' + sta + '/' + str(eventtime.year) + '/' + \
        str(eventtime.julday).zfill(3) + '/' + \
        loc + '_' + chan + '*.512.seed'
        
    st = read(path, starttime=eventtime-140., endtime=eventtime+120.)
    if chan == 'EJ':    
        st.decimate(2)
        st[0].stats.channel = 'EJN'
        st[0].data *= -1
        st[1].stats.channel = 'EJE'
        st[1].data *= -1
    else:
        st[0].stats.channel = 'HHN'
        st[1].stats.channel = 'HHE'
        # Time correction
        if eventtime.julday <= 143.:
            for tr in st:
                tr.stats.starttime += 1. 
    
    stD = st.copy()
    stA = st.copy()
    if (loc == '00') and (sta == 'OKR1'):
        pazD = paz00D
        paz = paz00
        pazA = paz00A
        word='Rotation'
    elif (loc =='10') and (sta == 'OKR1'):
        pazD = paz10D
        paz = paz10
        pazA = paz10A
        word ='Rotation'
    else:
        pazD = pazTCD
        paz = pazTC
        pazA = pazTCA
        word = 'Translation'
    if sta == 'OKR1':
        for idx, tr in enumerate(stD):
            tr.simulate(paz_remove=pazD[idx])
            
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=paz[idx])
            
        for idx, tr in enumerate(stA):
            tr.simulate(paz_remove=pazA[idx])
    else:
        for idx, tr in enumerate(stD):
            tr.simulate(paz_remove=pazTCD)
            
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=pazTC)
        for idx, tr in enumerate(stA):
            tr.simulate(paz_remove=pazTCA)
    st.rotate('NE->RT', bazi)
    stD.rotate('NE->RT', bazi)
    stA.rotate('NE->RT', bazi)
    signal = st.copy()
    noise = st.copy()
    stT= st.copy()
    
    stT.taper(.5)
    stD.taper(.5)
    stA.taper(.5)
    stD.filter('bandpass',freqmin = 1., freqmax=20.)
    stT.filter('bandpass',freqmin=1., freqmax=20.)
    stA.filter('bandpass', freqmin=1., freqmax=20.)
    stD.trim(starttime=eventtime-20.)
    stT.trim(starttime=eventtime-20.)
    stA.trim(starttime=eventtime-20.)
    fevent.write('Peak Radial ' + loc + ' ' + word + ': ' + str(abs(stD[0].max())) + '\n')
    fevent.write('Peak Transverse ' + loc + ' ' + word + ': ' + str(abs(stD[1].max())) + '\n')
    fevent.write('Peak Vertical ' + loc + ' ' + word + ': ' + str(abs(stD[2].max())) + '\n')
    fevent.write('Peak Radial ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[0].max())) + '\n')
    fevent.write('Peak Transverse ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[1].max())) + '\n')
    fevent.write('Peak Vertical ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[2].max())) + '\n')
    fevent.write('Peak Radial ' + loc + ' ' + word + ' RateRate: ' + str(abs(stA[0].max())) + '\n')
    fevent.write('Peak Transverse ' + loc + ' ' + word + ' RateRate: ' + str(abs(stA[1].max())) + '\n')
    fevent.write('Peak Vertical ' + loc + ' ' + word + ' RateRate: ' + str(abs(stA[2].max())) + '\n')
    stD.filter('bandpass',freqmin=1.,freqmax=20.)
    if sta == 'OKR1':
        sttemp = signal.copy()
        sttemp.filter('bandpass',freqmin=1.,freqmax=20.)
        PGV = np.max(np.sqrt((sttemp[0].data**2 +  sttemp[1].data**2)/2.))
        fevent.write('Peak ' + loc + ' PGWZ: ' + str(np.abs(np.max(stD[2].data))) + '\n')
        fevent.write('Peak ' + loc + ' PGWT: ' + str(np.abs(np.max(stD[1].data))) + '\n')
        fevent.write('PGwRsH' + loc + ': ' + str(PGV) + '\n')
        fevent.write('PGwRsZ' + loc + ': ' + str(np.abs(np.max(sttemp[2].data))) + '\n')
    else:
        sttemp = signal.copy()
        sttemp.filter('bandpass',freqmin=1.,freqmax=20.)
        PGV = np.max(np.sqrt((sttemp[0].data**2 +  sttemp[1].data**2)/2.))
        PGA = np.max(np.sqrt((stA[0].data**2 +  stA[1].data**2)/2.))
        fevent.write('PGVH: ' + str(PGV) + '\n')
        fevent.write('PGVZ: ' + str(np.abs(np.max(sttemp[2].data))) + '\n')
        fevent.write('PGAH: ' + str(PGA) + '\n')
        fevent.write('PGAZ: ' + str(np.abs(np.max(stA[2].data))) + '\n')

    noise.trim(endtime = eventtime-20)
    noiseRaw = noise.copy()
    signal.trim(starttime=eventtime)
    noise.filter('bandpass', freqmin=1., freqmax=5.)
    noise.taper(.5)
    signal.filter('bandpass', freqmin=1., freqmax=5.)
    signal.taper(.5)
    for pair in zip(noise, signal):
        SNR = pair[1].std()/pair[0].std()
        fevent.write('SNR ' + loc + ' ' + pair[0].stats.channel + ': ' + str(SNR) + '\n')
        st.trim(starttime=eventtime)
    # st is the signal without filtering
    # noiseRaw is the noise without filtering
    return st, noiseRaw




def plotTS(event, eveidx, st00R, st10R, st00):
    eventtime = event.origins[0].time
    eventtimeStr = str(eventtime.year) + ' ' + str(eventtime.julday).zfill(3) + ' ' + \
        str(eventtime.hour).zfill(2) + ':' + str(eventtime.minute).zfill(2) + ':' + \
        str(eventtime.second).zfill(2)
    mag = event.magnitudes[0].mag
    magstr = event.magnitudes[0].magnitude_type
    if 'Lg' in magstr:
        magstr = 'mb_{Lg}'
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude,event.origins[0].longitude)
    dis *= 1./1000.
    for comp in ['Vertical', 'Radial', 'Transverse']:
        fig = plt.figure(1)
        ax = plt.subplot(311)
        
        plt.title('Time Series ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' + magstr + '$=' + str(mag))
        
        if comp =='Vertical':
            idx = 2
        elif comp =='Radial':
            idx = 0
        elif comp == 'Transverse':
            idx = 1
        lim = 2.*np.abs((st00R[idx].max()))*10**6
        print(lim)
        t = np.arange(st00R[idx].stats.npts)/st00R[idx].stats.sampling_rate
        plt.plot(t,st00R[idx].data*10**6,label='00 Rotational ' + comp ,color='k')
        plt.legend(loc='upper right')
        plt.ylabel('$\mu$Radians/s')
        plt.xlim((0,120))
        plt.ylim((-lim,lim))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width, box.height * 0.8])
        ax= plt.subplot(312)
        t = np.arange(st10R[idx].stats.npts)/st10R[idx].stats.sampling_rate
        plt.plot(t,st10R[idx].data*10**6, label='10 Rotational '+ comp,color='k')
        plt.xlim((0,120))
        
        plt.ylim((-lim,lim))
        plt.ylabel('$\mu$Radians/s')
        plt.legend(loc='upper right')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width, box.height * 0.8])
        ax= plt.subplot(313)
        t = np.arange(st00[idx].stats.npts)/st00[idx].stats.sampling_rate
        plt.plot(t,st00[idx].data*10**6, label='00 Translational ' + comp ,color='k')
        plt.legend(loc='upper right')
        plt.xlabel('Time (s)')
        plt.xlim((0,120))
        plt.ylabel('$\mu$m/s')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3,
            box.width, box.height * 0.8])
        plt.savefig(comp + 'TimeSeries_Event' + str(eveidx+1) + 'CHECK.jpg', format='jpeg', dpi=400)
        plt.clf()
        plt.close()
    return



def plotNEW(event, eveidx, st00R, st10R, st00):
    eventtime = event.origins[0].time
    eventtimeStr = str(eventtime.year) + ' ' + str(eventtime.julday).zfill(3) + ' ' + \
        str(eventtime.hour).zfill(2) + ':' + str(eventtime.minute).zfill(2) + ':' + \
        str(eventtime.second).zfill(2)
    mag = event.magnitudes[0].mag
    magstr = event.magnitudes[0].magnitude_type
    if 'Lg' in magstr:
        magstr = 'mb_{Lg}'
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude,event.origins[0].longitude)
    dis *= 1./1000.
    
    fig = plt.figure(1,figsize=(16,6))
    plt.suptitle('Time Series ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' + magstr + '$' + str(mag))
    plt.subplots_adjust(hspace=0.001)
    #plt.subplots_adjust(wspace=0.002)
    ax = plt.subplot(321)
    plt.suptitle('Time Series ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' + magstr + '$' + str(mag))

    lim = 2.*np.abs((st00R[0].max()))*10**6
    print(lim)
    t = np.arange(st00R[0].stats.npts)/st00R[0].stats.sampling_rate
    plt.plot(t,st00R[0].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[0].data*10**6,label='10 Rotational ', linewidth=1.)
    plt.text(2.5, .3, 'Radial')
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.xticks([])
    plt.xlim((0,120))
    ax = plt.subplot(323)
    plt.text(2.5, .3, 'Transverse')
    t = np.arange(st00R[1].stats.npts)/st00R[1].stats.sampling_rate
    plt.plot(t,st00R[1].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[1].data*10**6,label='10 Rotational ', linewidth=1.)
    plt.ylabel('$\mu$Radians/s')
    plt.xlim((0,120))
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.xticks([])
    ax = plt.subplot(325)
    t = np.arange(st00R[2].stats.npts)/st00R[1].stats.sampling_rate
    plt.plot(t,st00R[2].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[2].data*10**6,label='10 Rotational ', linewidth=1.)
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.text(2.5, .5, 'Vertical')
    plt.xlim((0,120))
    plt.xlabel('Time (s)')
    #plt.xticks([])
    ax = plt.subplot(322)
    #plt.title('Time Series ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' + magstr + '$' + str(mag))
    plt.xlabel('Time (s)')
    lim = 2.*np.abs((st00R[0].max()))*10**6
    print(lim)
    t = np.arange(st00R[0].stats.npts)/st00R[0].stats.sampling_rate
    plt.plot(t,st00R[0].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[0].data*10**6,label='10 Rotational ', linewidth=1.)
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.xticks([])
    plt.yticks([])
    plt.xlim((35.,45.))
    ax = plt.subplot(324)
    t = np.arange(st00R[1].stats.npts)/st00R[1].stats.sampling_rate
    plt.plot(t,st00R[1].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[1].data*10**6,label='10 Rotational ', linewidth=1.)
    #plt.ylabel('$\mu$Radians/s')
    plt.xlim((35.,45.))
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.xticks([])
    plt.yticks([])
    ax = plt.subplot(326)
    t = np.arange(st00R[2].stats.npts)/st00R[1].stats.sampling_rate
    plt.plot(t,st00R[2].data*10**6,label='00 Rotational ', linewidth=1.)
    plt.plot(t,st10R[2].data*10**6,label='10 Rotational ', linewidth=1.)
    ax.axvspan(35., 45., alpha=0.5, color='.5')
    plt.xlim((35.,45.))
    plt.xlabel('Time (s)')
    plt.yticks([])
    #plt.tight_layout()
    #fig.tight_layout(rect=[0,0,.8,1]) 
    plt.savefig('NEWTimeSeries_Event' + str(eveidx+1) + 'CHECK.jpg', format='jpeg', dpi=400)
    plt.clf()
    plt.close()
    
    
    #plt.show()
    #plt.legend(loc='upper right')
    #
    #plt.xlim((0,120))
    #plt.ylim((-lim,lim))
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.3,
             #box.width, box.height * 0.8])
    #ax= plt.subplot(312)
    #t = np.arange(st10R[idx].stats.npts)/st10R[idx].stats.sampling_rate
    #plt.plot(t,st10R[idx].data*10**6, label='10 Rotational '+ comp,color='k')
    #plt.xlim((0,120))
    
    #plt.ylim((-lim,lim))
    #plt.ylabel('$\mu$Radians/s')
    #plt.legend(loc='upper right')
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.3,
             #box.width, box.height * 0.8])
    #ax= plt.subplot(313)
    #t = np.arange(st00[idx].stats.npts)/st00[idx].stats.sampling_rate
    #plt.plot(t,st00[idx].data*10**6, label='00 Translational ' + comp ,color='k')
    #plt.legend(loc='upper right')
    #plt.xlabel('Time (s)')
    #plt.xlim((0,120))
    #plt.ylabel('$\mu$m/s')
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.3,
        #box.width, box.height * 0.8])
    #plt.savefig(comp + 'TimeSeries_Event' + str(eveidx+1) + 'CHECK.jpg', format='jpeg', dpi=400)
    #plt.clf()
    #plt.close()
    return




# Map of events
cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=2.,endtime=etime)

inventory = client.get_stations(network="GS", station="OKR1")

figmap = cat.plot(projection="local", resolution="i", show=False, title='Color by Depth (km)')
figmap = inventory.plot(projection="local", show=False, fig=figmap, size=350, color="r", legend=None)

figmap.savefig('Eventsplot.jpg',format='jpeg')
figmap.clf()

print(len(cat))


# Now go through each event


for eveidx, event in enumerate(cat):
    #if eveidx != 114:
    #    continue
    print('On event: ' + str(eveidx+1))
    fevent = open('OTHERResults_Event_' + str(eveidx+1),'w')
    bazi = writeEventInfo(event,fevent)
    
    try:
    #if True:
        st00R, n00R = grabdata(event.origins[0].time, '00', 'OKR1', bazi, fevent)
    except:
        print('Unable to get 00 R data')
    try:
    
        st10R, n10R = grabdata(event.origins[0].time, '10', 'OKR1', bazi, fevent)
    except:
        print('Unable to get 10 R data')
    try:
    #if True:
        st00, n00 = grabdata(event.origins[0].time, '00', 'OK038', bazi, fevent)
    except:
        print('Unable to get 00 translational data')
    
    # Now we have all the data for the events
    



    try:
    #if True:
        st00R.filter('bandpass', freqmin=1., freqmax=5.)
        st00R.taper(.5)
        st00.filter('bandpass', freqmin=1., freqmax=5.)
        st00.taper(.5)
        st10R.filter('bandpass',freqmin=1., freqmax=5.)
        st10R.taper(.5)
        #try:
        if True:
            plotNEW(event, eveidx, st00R, st10R, st00)
    except:
        print('Unable to plot T series')
    fevent.close()
    
    
    
    
    
    
    
    
    
    
    
    

