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
    path = '/tr1/telemetry_days/GS_' + sta + '/' + str(eventtime.year) + '/' + \
        str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '/' + \
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
    if (loc == '00') and (sta == 'OKR1'):
        pazD = paz00D
        paz = paz00
        word='Rotation'
    elif (loc =='10') and (sta == 'OKR1'):
        pazD = paz10D
        paz = paz10
        word ='Rotation'
    else:
        pazD = pazTCD
        paz = pazTC
        word = 'Translation'
    if sta == 'OKR1':
        for idx, tr in enumerate(stD):
            tr.simulate(paz_remove=pazD[idx])
            
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=paz[idx])
    else:
        for idx, tr in enumerate(stD):
            tr.simulate(paz_remove=pazTCD)
            
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=pazTC)
    st.rotate('NE->RT', bazi)
    stD.rotate('NE->RT', bazi)
    signal = st.copy()
    noise = st.copy()
    stT= st.copy()
    
    stT.taper(.5)
    stD.taper(.5)
    stD.filter('bandpass',freqmin = 1., freqmax=20.)
    stT.filter('bandpass',freqmin=1., freqmax=20.)
    stD.trim(starttime=eventtime-20.)
    stT.trim(starttime=eventtime-20.)
    fevent.write('Peak Radial ' + loc + ' ' + word + ': ' + str(abs(stD[0].max())) + '\n')
    fevent.write('Peak Transverse ' + loc + ' ' + word + ': ' + str(abs(stD[1].max())) + '\n')
    fevent.write('Peak Vertical ' + loc + ' ' + word + ': ' + str(abs(stD[2].max())) + '\n')
    fevent.write('Peak Radial ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[0].max())) + '\n')
    fevent.write('Peak Transverse ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[1].max())) + '\n')
    fevent.write('Peak Vertical ' + loc + ' ' + word + ' Rate: ' + str(abs(stT[2].max())) + '\n')
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
        fevent.write('PGVH: ' + str(PGV) + '\n')
        fevent.write('PGVZ: ' + str(np.abs(np.max(sttemp[2].data))) + '\n')
    

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

def makePSD(st,noise, event, eveidx, debug=True):
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
    ps=[]
    ns=[]
    fig =plt.figure(1)
    for tr in st:
        (p,f)=cp(tr,tr,length,0,tr.stats.delta)
        p = konno_ohmachi_smoothing(p,f, bandwidth=15.)
        plt.semilogx(f,10.*np.log10(p),label='Signal: ' + tr.id)
        ps.append(10*np.log10(p))
    for tr in noise:
        (n,f)=cp(tr,tr,length,0,tr.stats.delta)
        n = konno_ohmachi_smoothing(n,f, bandwidth=15.)
        plt.semilogx(f,10.*np.log10(n),label='Noise: ' + tr.id)
        ns.append(10.*np.log10(n))
    plt.xlim((.5,50))
    plt.xlabel('Frequency (Hz)')
    if st[0].stats.station == 'OKR1':
        plt.ylabel('dB (rel. 1 $(radian/s^2)^2/Hz$)')
        word='R'
    else:
        plt.ylabel('dB (rel. 1 $(m/s^2)^2/Hz$)')
        word ='T'
    plt.legend()
    plt.title(eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' +  magstr + '$=' + str(mag))
    plt.savefig('SmoothPSD' + word + st[0].stats.location +'_Event' + str(eveidx + 1) + '.jpg', format='jpeg', dpi=400)
    plt.clf()
    plt.close()
    return ps, ns, f


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
        plt.savefig(comp + 'TimeSeries_Event' + str(eveidx+1) + '.jpg', format='jpeg', dpi=400)
        plt.clf()
        plt.close()
    return



def MakePeakPlot(event,eveidx):
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
    path1 = '/tr1/telemetry_days/GS_OKR1/' + str(eventtime.year) + '/' + \
        str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '/' + \
        '*_EJZ.512.seed'
        
    st1 = read(path1, starttime=eventtime-140., endtime=eventtime+120.)
    
    path2='/tr1/telemetry_days/GS_OK038/' + str(eventtime.year) + '/' + \
        str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '/' + \
        '00_HH*.512.seed'
    st2 = read(path2, starttime=eventtime-140., endtime=eventtime+120.)
    st2[0].stats.channel = 'HHN'
    st2[1].stats.channel = 'HHE'
    # Time correction
    if eventtime.julday <= 143.:
        for tr in st2:
            tr.stats.starttime += 1. 
    
    for tr in st2:
        tr.simulate(paz_remove=pazTC)
    st1[0].simulate(paz00D[2])
    st1[1].simulate(paz10D[2])

    
    st2.rotate('NE->RT', bazi)

    
    st2.taper(.5)
    st1.taper(.5)
    st2.filter('bandpass',freqmin = 1., freqmax=20.)
    st1.filter('bandpass',freqmin=1., freqmax=20.)
    st2.trim(starttime=eventtime-20.)
    st1.trim(starttime=eventtime-20.)
    lim = 3.*np.abs(st1[0].max())*10**6
    fig = plt.figure(1)
    ax = plt.subplot(311)
    plt.title('Time Series ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' + magstr + '$=' + str(mag))
    t = np.arange(st1[0].stats.npts)/st1[0].stats.sampling_rate
    t+= -20.
    plt.plot(t,st1[0].data*10**6,label='00 Rotational Vertical',color='k')
    plt.ylim((-lim,lim))
    plt.legend(loc='upper left', fontsize=8)
    plt.ylabel('$\mu$Radians')
    plt.xlim((0,120))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
        box.width, box.height * 0.8])
    ax= plt.subplot(312)
    t = np.arange(st1[1].stats.npts)/st1[1].stats.sampling_rate
    t+= -20.
    plt.plot(t,st1[1].data*10**6, label='10 Rotational Vertical',color='k')
    plt.xlim((0,120))
    plt.ylim((-lim,lim))
    plt.ylabel('$\mu$Radians')
    plt.legend(loc='upper left', fontsize=8)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
        box.width, box.height * 0.8])
    ax= plt.subplot(313)
    t = np.arange(st2[0].stats.npts)/st2[0].stats.sampling_rate
    t+= -20.
    plt.plot(t,(np.sqrt((st2[0].data**2 +  st2[1].data**2))/2.)*10**6, label='00 Horizontal' ,color='k')
    plt.legend(loc='upper left')
    plt.xlabel('Time (s)')
    plt.xlim((0,120))
    plt.ylabel('$\mu$m/s')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
        box.width, box.height * 0.8])
    plt.savefig('PGVROTTimeSeries_Event' + str(eveidx+1) + '.jpg', format='jpeg', dpi=400)
    plt.clf()
    plt.close()
    

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
    if eveidx != 114:
        continue
    print('On event: ' + str(eveidx+1))
    fevent = open('Results_Event_' + str(eveidx+1),'w')
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
        st00, n00 = grabdata(event.origins[0].time, '00', 'OK038', bazi, fevent)
    except:
        print('Unable to get 00 translational data')
    
    # Now we have all the data for the events
    
    try:

        ps00R,ns00R, f = makePSD(st00R, n00R, event, eveidx)
    except:
        print('Unable to make 00 R psd')

    try:
        ps10R,ns10R, f = makePSD(st10R, n10R, event, eveidx)
    except:
        print('Unable to make 10 R psd')
    
    try:
        ps00T, ns10T, f = makePSD(st00, n00, event, eveidx)
    except:
        print('Unable to make 00 T psd')


    try:

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
        # Now take spectral ratios and plot
        fig = plt.figure(1)
        for idx, comp in enumerate(['Radial', 'Transverse', 'Vertical']):
            plt.semilogx(f, ps00T[idx] - ps00R[idx], label='00 ' + comp)
            plt.semilogx(f, ps00T[idx] - ps10R[idx], label='10 ' + comp)
        plt.xlim((.5,50))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Ratio (dB)')
        plt.legend()
        plt.title('Spectral Ratio ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' +  magstr + '$=' + str(mag))
        plt.savefig('PSDRatio_Event' + str(eveidx + 1) + '.jpg', format='jpeg', dpi=400)
        plt.clf()
        plt.close()
    except:
        plt.clf()
        plt.close()
        print('Unable to compute splectral ratios')

    try:
        minfre =1.
        maxfre = 10.
        for idx, comp in enumerate(['Radial','Transverse', 'Vertical']):
            p00DIFF = ps00T[idx] -ps00R[idx]
            p10DIFF = ps00T[idx] -ps10R[idx]
            p00TTemp = ps00T[idx][(minfre <= f) & (maxfre >= f)]
            p00RTemp = ps00R[idx][(minfre <= f) & (maxfre >= f)]
            p10RTemp = ps10R[idx][(minfre <= f) & (maxfre >= f)]
            p00DIFF = p00DIFF[(minfre <= f) & (maxfre >= f)]
            p10DIFF = p10DIFF[(minfre <= f) & (maxfre >= f)]
            fT = f[(minfre <= f) & (maxfre >= f)]
            fevent.write('Mean Spectral Ratio 00 ' + comp + ': ' + str(np.mean(p00DIFF)) + '\n')
            fevent.write('Peak Frequency Ratio 00 ' + comp + ': ' + str(fT[np.argmax(p00DIFF)]) + '\n')
            fevent.write('Mean Spectral Ratio 00 : '+ str(np.mean(p00TTemp)) + '\n')
            fevent.write('Mean Spectral Ratio 10 ' + comp + ': ' + str(np.mean(p10DIFF)) + '\n')
            fevent.write('Peak Frequency Ratio 10 ' + comp + ': ' + str(fT[np.argmax(p10DIFF)]) + '\n')
            fevent.write('Mean Spectral 00 ' + comp + ': '+ str(np.mean(p00TTemp)) + '\n')
            fevent.write('Mean Spectral 00 ' + comp +  ' Translational: ' + str(np.mean(p00TTemp)) + '\n')
            fevent.write('Mean Spectral 00 ' + comp +  ' Rotational: ' + str(np.mean(p00RTemp)) + '\n')
            fevent.write('Mean Spectral 10 ' + comp +  ' Rotational: ' + str(np.mean(p10RTemp)) + '\n')
            fevent.write('Frequency Peak 00 ' + comp + ': '+ str(fT[np.argmax(p00TTemp)]) + '\n')
            fevent.write('Frequency Peak 00 ' + comp + ' Rotational: '  + str(fT[np.argmax(p00RTemp)]) + '\n')
            fevent.write('Frequency Peak 10 ' + comp + ' Rotational: '  + str(fT[np.argmax(p10RTemp)]) + '\n')
            fevent.write('Peak Spectral Value 00 Rotational ' + comp + ': ' + str(np.max(p00RTemp)) + '\n')
            fevent.write('Peak Spectral Value 10 Rotational ' + comp + ': ' + str(np.max(p10RTemp)) + '\n')
            fevent.write('Peak Spectral Value 00 Translational ' + comp + ': ' + str(np.max(p00TTemp)) + '\n')
            fevent.write('Peak Spectral Frequency Value 00 Rotational ' + comp + ': ' + str(fT[np.argmax(p00RTemp)]) + '\n')
            fevent.write('Peak Spectral Frequency Value 10 Rotational ' + comp + ': ' + str(fT[np.argmax(p10RTemp)]) + '\n')
            fevent.write('Peak Spectral Frequency Value 00 Translational ' + comp + ': ' + str(fT[np.argmax(p00TTemp)]) + '\n')
    except:
        print('problem with spectral calc')
    
    
    
    
    try:
        MakePeakPlot(event,eveidx)
    except:
        print('Problem with PG plot')
    
    try:
        ###################################### Make a spectral ratio plot for each event
        fig = plt.figure(1, figsize=(8,8))
        plt.semilogx(f, ps00T[2], label='Vertical Velocity')
        plt.semilogx(f, (ps00T[1] + ps00T[0])/2., label='Mean Horizontal Velocity')
        plt.semilogx(f, ps00R[2], label='00 Vertical Rotation Rate')
        plt.semilogx(f, (ps00R[1] + ps00R[0])/2., label='00 Mean Horizontal Rotation Rate')
        plt.semilogx(f, ps10R[2], label='10 Vertical Rotation Rate')
        plt.semilogx(f, (ps10R[1] + ps10R[0])/2., label='10 Mean Horizontal Rotation Rate')
        plt.axvspan(1.,10., alpha=.5, color='0.75')
        plt.title('PSD ' + eventtimeStr + ' Distance:' + str(int(dis)) + ' km $' +  magstr + '$=' + str(mag))
        plt.xlim((0.5,50))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power (dB)')
        plt.legend(loc='upper right')

        plt.savefig('PSDMean_Event' + str(eveidx + 1) + '.jpg', format='jpeg', dpi=400)
        plt.clf()
        plt.close()
    except:
        print('Problem with spectral plot')
        
    
    
    
    # Last thing to do is filter and plot

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
            plotTS(event, eveidx, st00R, st10R, st00)
    except:
        print('Unable to plot T series')
    fevent.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ######################################################################
    ## compute phase velocity
    
    #try:
    ##if True:
        #minfre=.1
        #maxfre=20.
        #ph = .5*(10**(pTC[0]/20.)*(2*np.pi*f))/(10**(p00[2]/20))
        #ph2 =.5*(10**(pTC[0]/20.)*(2*np.pi*f))/(10**(p10[2]/20))
        #f2 = f[(minfre <= f) & (maxfre >= f)]
        #ph = ph[(minfre <= f) & (maxfre >= f)]
        #ph2 = ph2[(minfre <= f) & (maxfre >= f)]
        #fig = plt.figure(9)
        #plt.semilogx(f2,ph/1000., label='00 Phase Velocity')
        
        #plt.semilogx(f2,ph2/1000., label='10 Phase Velocity')
        #plt.xlim((.5,20))
        #plt.xlabel('Frequency (Hz)')
        #plt.ylabel('Phase Velocity (km/s)')
        #plt.legend()
        #plt.title(eventtimeStr + ' Distance:' + str(int(dis)) + ' km ' +  event.magnitudes[0].magnitude_type + '=' + str(mag))
        #plt.savefig('PhaseVel_' + str(eveidx+1) + '.jpg', format='jpeg')
        #plt.clf()
        
        
        
        
        #phs1.append(ph)
        #phs2.append(ph2)
        #fgood = f2
        
        
        
        
        
        #nminfre = 1.
        #nmaxfre = 3.
        #ph = ph[(nminfre <=f2) &(nmaxfre >= f2)]
        #ph2 = ph2[(nminfre <=f2) &(nmaxfre >= f2)]
        #fevent.write('Mean Phase Velocity 00: ' + str(np.mean(ph)) + '\n')
        #fevent.write('Mean Phase Velocity 10: ' + str(np.mean(ph2)) + '\n')
        #fevent.write('Max Phase Velocity 00: ' + str(np.max(ph)) + '\n')
        #fevent.write('Max Phase Velocity 10: ' + str(np.max(ph2)) + '\n')
        #fevent.write('Peak Phase Frequency 00: ' + str(f2[np.argmax(ph)]) + '\n')
        #fevent.write('Peak Phase Frequency 10: ' + str(f2[np.argmax(ph2)]) + '\n')
        
        
        
        
        
        
    #except:
        #print('Problem with phase velocity')
    
    
    
    

   
    
#fig = plt.figure(10)
#for idx, ele in enumerate(zip(phs1,phs2)):
    #if idx==0:
        #m1 = ele[0]
        #m2 = ele[1]
    #else:
        #m1 += ele[0]
        #m2 += ele[1]
    #plt.semilogx(fgood,ele[0]/1000., alpha=.3, color='.5')
    #plt.semilogx(fgood,ele[1]/1000., alpha=.3,color='.5')

#plt.semilogx(fgood, (m1/float(len(phs1)))/1000.,color='k')
#plt.semilogx(fgood, (m2/float(len(phs2)))/1000.,color='k')
#plt.xlim((.5,20))
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Phase Velocity (km/s)')
#plt.savefig('PHASEALL.jpg',format='jpeg')
