#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client
import sys
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.signal.cross_correlation import xcorr
from itertools import combinations, permutations
import multiprocessing
from functools import partial
##########################################


def grabdata(starttime, endtime, loc, sta):
    if sta == 'OKR1':
        chan = 'EJ'
    else:
        chan = 'HH'
    path = '/tr1/telemetry_days/GS_' + sta + '/' + str(starttime.year) + '/' + \
        str(starttime.year) + '_' + str(starttime.julday).zfill(3) + '/' + \
        loc + '_' + chan + '*Z.512.seed'
    st = read(path, starttime=starttime, endtime=endtime)
    #st.merge()
    st.decimate(5)
    if chan == 'EJ':    
        st.decimate(2)
    else:
        # Time correction
        if starttime.julday <= 143.:
            for tr in st:
                tr.stats.starttime += 1. 
    st.filter('bandpass',freqmin = 2., freqmax=8.)
    return st

def docorr(st, triple):
    template=ele[0]
    event=ele[1]
    idx= ele[2]
    vals = []
    for stC in st.slide(120., 0.5):
        try:
            lag, val = xcorr(stC[0],template[0],100)
            vals.append(val)
        except:
            vals.append(0.)
    fstring = 'TemplateResults_' + stC[0].stats.location + '_' 
    fstring += stC[0].stats.station + '_' + str(idx) + '_' 
    fstring += str(st[0].stats.starttime.julday)
    f=open(fstring,'w')
    f.write('lat:' + str(event.origins[0].latitude) + '\n')
    f.write('lon:' + str(event.origins[0].longitude) + '\n')
    f.write('time:' + str(event.origins[0].time) + '\n')
    f.write('Magnitude:' + str(event.magnitudes[0].mag) + '\n')
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude, event.origins[0].longitude)
    f.write('distance:' + str( dis/1000.) + '\n')
    f.write('bazi:' + str(bazi) + '\n')
    f.write('data\n')
    for val in vals:
        f.write(str(val) + '\n')
    f.close()
    return 
        
 
debug = True
       
stime = UTCDateTime('2017-118T00:00:00.0')
etime = UTCDateTime('2017-128T00:00:00.0')
client = Client("IRIS")

stalat = 36.47819
stalon = -98.742240



stas = ['OK038','OKR1']
locs = ['00','10']



# Grab all the events
deg20km = kilometer2degrees(20.)
cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=deg20km)


for ele in cat:
    print(ele)
sys.exit()



if debug:
    print('We have ' + str(len(cat)) + ' initial events')
# Now populate the templates as 120 second events

for sta, loc in zip(stas, locs):
    print('One station: ' + loc + ' ' + sta)
    templates = []
    goodevents =[]
    for idx, event in enumerate(cat):
        # We have our template event
        time = event.origins[0].time
        try:
            st=grabdata(time, time + 120., loc, sta)
            templates.append(st)
            goodevents.append(event)
        except:
            print('ditching event number: ' + str(idx))

    if debug:
        print('Number of templates: ' + str(len(templates)))
        print('Number of events: ' + str(len(goodevents))) 


    ctime=stime
    sts = []
    while ctime <= etime:
        if debug:
            print('On current time of: ' + str(ctime))
        sts.append(grabdata(ctime,ctime+24.*60.*60.,loc, sta))
        ctime += 24.*60.*60.
    print('Here are the number of streams: ' + str(len(sts)))
    eles = zip(templates, goodevents, range(len(templates)))
    for ele in eles:
        print('On ele of: ' + str(ele[2]))
        pool = multiprocessing.Pool(10)
        pool.map(partial(docorr, triple=ele), sts)
        pool.close()


