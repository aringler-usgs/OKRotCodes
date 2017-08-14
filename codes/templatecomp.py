#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
from obspy.signal.invsim import paz_to_freq_resp
from matplotlib.mlab import csd
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from obspy.imaging.maps import plot_basemap
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
from obspy.signal.cross_correlation import xcorr_3c, xcorr
from itertools import combinations, permutations
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



debug = True
length=4096

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


def grabdata(eventtime, loc, sta, bazi):
    if sta == 'OKR1':
        chan = 'EJ'
    else:
        chan = 'HH'
    path = '/tr1/telemetry_days/GS_' + sta + '/' + str(eventtime.year) + '/' + \
        str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '/' + \
        loc + '_' + chan + '*.512.seed'
        
    st = read(path, starttime=eventtime, endtime=eventtime+120.)
    st.merge()
    st.sort(['channel'])
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

    if (loc == '00') and (sta == 'OKR1'):
        paz = paz00
        word='Rotation'
    elif (loc =='10') and (sta == 'OKR1'):
        paz = paz10
        word ='Rotation'
    else:

        paz = pazTC
        word = 'Translation'
    if sta == 'OKR1':
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=paz[idx])
    else:
        for idx, tr in enumerate(st):
            tr.simulate(paz_remove=pazTC)

    st.rotate('NE->RT', bazi)

    #st.taper(.5)
    st.filter('bandpass',freqmin = 2., freqmax=8.)

    return st


def crossCorr(event1, event2, fcomp, debug =False):
    eventtime1 = event1.origins[0].time
    mag = event1.magnitudes[0].mag
    (dis1,azi, bazi1) = gps2dist_azimuth(stalat, stalon, event1.origins[0].latitude,event1.origins[0].longitude)
    dis1 *= 1./1000.
    
    eventtime2 = event2.origins[0].time
    mag = event2.magnitudes[0].mag
    (dis2,azi, bazi2) = gps2dist_azimuth(stalat, stalon, event2.origins[0].latitude,event2.origins[0].longitude)
    dis2 *= 1./1000.
    
    (evedis,azi, bazi) = gps2dist_azimuth(event1.origins[0].latitude, event1.origins[0].longitude, event2.origins[0].latitude,event2.origins[0].longitude)
    evedis *= 1./1000.
    
    if debug:
        print('Distance 1: ' + str(dis1) + ' Distance 2: ' + str(dis2))
    
    stas=['00,OK038','00,OKR1','10,OKR1']
    Rot1 = 0.
    Rot2 = 0.
    for idx, sta in enumerate(stas):
        loc = sta.split(',')[0]
        sta = sta.split(',')[1]
        # Now cross-correlate the pairs
        try:
        #if True:
            st1 = grabdata(eventtime1, loc, sta, bazi1)
            st2 = grabdata(eventtime2, loc, sta, bazi2)
            lag, val = xcorr_3c(st1, st2, 1000, components=['Z', 'R', 'T'])
            string = sta + ', ' + loc + ', ' + str(dis1) + ', ' + str(dis2) + ', ' + str(evedis) + ', ' + str(val) + '\n'
            fcomp.write(string)
            if idx == 0:
                Rot1 = val
                Rot2 = val
            if idx == 1:
                Rot1 += val
                string = '6' +  sta + ', ' + loc + ', ' + str(dis1) + ', ' + str(dis2) + ', ' + str(evedis) + ', ' + str(Rot1/2.) + '\n'
                fcomp.write(string)
            if idx == 2:
                Rot2 += val
                string = '6' +  sta + ', ' + loc + ', ' + str(dis1) + ', ' + str(dis2) + ', ' + str(evedis) + ', ' + str(Rot2/2.) + '\n'
                fcomp.write(string)
        except:
            print('Unable to correlate ' + sta + ' ' + loc)
    
    return

def crossCorrTR(event):
    eventtime = event.origins[0].time
    mag = event.magnitudes[0].mag
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude,event.origins[0].longitude)
    dis *= 1./1000.
    st1 = grabdata(eventtime, '00', 'OK038', bazi)
    st2 = grabdata(eventtime, '00', 'OKR1', bazi)
    for ele in permutations([0,1,2],2):
        print(ele)
        lag, val = xcorr(st1[ele[0]], st2[ele[1]], 1000)
        string =st1[ele[0]].stats.channel + ', ' + st2[ele[1]].stats.channel + ', ' + str(dis) + ', ' + str(val) + '\n'
        print(string)
        
        
        


# Map of events
cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=2., endtime=etime)







fcomp = open('corrs','w')
idx = 1
for ele in combinations(cat,2):
    print('On events pair :'  + str(idx))
    idx += 1 
    event1 = ele[0]
    event2 = ele[1]
    (evedis,azi, bazi) = gps2dist_azimuth(event1.origins[0].latitude, event1.origins[0].longitude, event2.origins[0].latitude,event2.origins[0].longitude)
    evedis *= 1./1000.
    
    if evedis < 20.:
        crossCorr(event1, event2, fcomp)
    
    #except:
    #    print('Problem with events')
    #sys.exit()
    
fcomp.close()






