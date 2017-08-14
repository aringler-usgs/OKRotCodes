#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
from obspy.signal.invsim import paz_to_freq_resp
from matplotlib.mlab import csd
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.imaging.maps import plot_basemap
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
from obspy.signal.cross_correlation import xcorr_3c, xcorr
from itertools import combinations, permutations, product
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




def grabdata(eventtime, etime, loc, sta, bazi):
    if sta == 'OKR1':
        chan = 'EJ'
    else:
        chan = 'HH'
    path = '/tr1/telemetry_days/GS_' + sta + '/' + str(eventtime.year) + '/' + \
        str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '/' + \
        loc + '_' + chan + '*.512.seed'
        
    st = read(path, starttime=eventtime, endtime=etime)
    st.merge()
    st.sort(['channel'])
    if chan == 'EJ':    

        st.decimate(2)
        st.decimate(5)
        st[0].stats.channel = 'EJN'
        st[0].data *= -1
        st[1].stats.channel = 'EJE'
        st[1].data *= -1
    else:

        st.decimate(5)
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
    stimeC = max([tr.stats.starttime for tr in st])
    etimeC = min([tr.stats.endtime for tr in st])
    st.trim(starttime=stimeC, endtime=etimeC)
    for tr in st:
        if tr.stats.npts > 2400:
            tr.data = tr.data[:2400]
        if tr.stats.npts < 2400:
            tr.data = np.append(tr.data, np.zeros(2400 - tr.stats.npts))
    st.rotate('NE->RT', bazi)

    st.filter('bandpass',freqmin = 2., freqmax=8.)

    return st

def addCorrs(corrs, idx, event):
    corrs[idx].update({'event':idx})
    corrs[idx].update({'evelat': event.origins[0].latitude})
    corrs[idx].update({'evelon': event.origins[0].longitude})
    corrs[idx].update({'time': event.origins[0].time})
    corrs[idx].update({'mag': event.magnitudes[0].mag})
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude, event.origins[0].longitude)
    corrs[idx].update({'distance': dis/1000.})
    corrs[idx].update({'bazi': bazi})
    corrs[idx].update({'data':[]})
    return corrs
    
def doTemplate(templates,loc,sta, corrs, stime, etime):
    for idx, template in enumerate(templates):
        currenttime = stime
        try:
        #if True:    
            st = grabdata(currenttime, etime, loc, sta, corrs[idx]['bazi'])
        except:
            print('Didnt get the data')
            continue
            
        # We have the data now lets go through it
        while currenttime <= etime:
            try:
            #if True:
                stC = st.copy()
                stC.trim(currenttime, currenttime + 120.)

                lag, val = xcorr_3c(template, stC, 100, components=['Z', 'R', 'T'])
                corrs[idx]['data'].append(val)

            except:
                print('Have a problem')

       
            currenttime += 0.5
            #print('Here is our current time:' + str(currenttime))
    return corrs
    
def writeresults(corrs, n, loc, sta, stime2):
    for idx in range(n):
        f=open('TemplateResults_' + loc + '_' + sta + '_' + str(idx) + '_' + str(stime2.julday) ,'w')
        f.write('lat:' + str(corrs[idx]['evelat']) + '\n')
        f.write('lon:' + str(corrs[idx]['evelon']) + '\n')
        f.write('time:' + str(corrs[idx]['time']) + '\n')
        f.write('Magnitude:' + str(corrs[idx]['mag']) + '\n')
        f.write('distance:' + str(corrs[idx]['distance']) + '\n')
        f.write('bazi:' + str(corrs[idx]['bazi']) + '\n')
        f.write('data\n')
        for ele in corrs[idx]['data']:
            f.write(str(ele) + '\n')
        f.close()
    return

def doALL(time, debug= True):
    etime = time + 2.*60.*60
    doTemplate(templates, loc, sta, corrs, time, etime)
        
    return


stime = UTCDateTime('2017-118T00:00:00.0')
etime = UTCDateTime('2017-158T00:00:00.0')
client = Client("IRIS")

stalat = 36.47819
stalon = -98.742240

loc = '00'
sta ='OKR1'

     

# Grab all the events
deg20km = kilometer2degrees(20.)
cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=deg20km)


# Now populate the templates as 120 second events
templates00T = []
templates00R = []
templates10R =[]
corrs = [dict() for x in range(len(cat))]

for idx, event in enumerate(cat):
    # We have our template event
    corrs = addCorrs(corrs,idx,event)
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, event.origins[0].latitude, event.origins[0].longitude)
    time = event.origins[0].time
    print(time)
    time.precision=1.
    print(time)
    try:
    #if True:
        st=grabdata(time, time + 120., '00', 'OK038', bazi)
        
        templates00R.append(st)
        st=grabdata(time, time + 120., '00', 'OKR1', bazi)
        templates00T.append(st)
        st=grabdata(time, time + 120., '10', 'OKR1', bazi)
        templates10R.append(st)
    except:
        print('ditching an event')


# We have a list of streams
# We have a list of corrs



# We now have a list of streams and a list of dictionaries



from matplotlib.mlab import griddata
xl=[]
yl=[]
vals00T=[]
vals00R=[]
vals10R=[]
good=0
for ele in product(range(len(templates00T)),repeat=2):
    if True:
        print(templates00T[ele[0]])
        print(templates00T[ele[1]])
        lag, val = xcorr_3c(templates10R[ele[0]], templates10R[ele[1]], 100, components=['Z', 'R', 'T'])
        
        vals10R.append(val)
        xl.append(ele[0])
        yl.append(ele[1])

        lag, val = xcorr_3c(templates00T[ele[0]], templates00T[ele[1]], 100, components=['Z', 'R', 'T'])
        vals00T.append(val)
        lag, val = xcorr_3c(templates00R[ele[0]], templates00R[ele[1]], 100, components=['Z', 'R', 'T'])
        vals00R.append(val)
        good +=1
    #except:
    #    pass
from scipy.cluster import hierarchy
from scipy.spatial import distance
   
xl=np.asarray(xl)
yl=np.asarray(yl)
vals00T=np.asarray(vals00T)
vals00R=np.asarray(vals00R)
vals10R=np.asarray(vals10R)
#check = vals00T-vals00R  
X = xl.reshape((len(templates00T),len(templates00T)))
Y = yl.reshape((len(templates00T), len(templates00T)))
Z = vals00T.reshape((len(templates00T),len(templates00T)))
#Z = check.reshape((len(templates00T),len(templates00T)))
fig = plt.figure(1,figsize=(12,12))


diss = distance.squareform(1-Z)
t= 0.3
linkage=hierarchy.linkage(diss)
clusters = hierarchy.fcluster(linkage, t, criterion="distance")

fig = plt.figure(1)
hierarchy.dendrogram(linkage, color_threshold=0.3)

#plt.subplot(311)
#plt.pcolormesh(X,Y,Z)
#plt.subplot(312)
#Z = vals00R.reshape((len(templates00T),len(templates00T)))
#plt.pcolormesh(X,Y,Z)
#plt.subplot(313)
#Z = vals10R.reshape((len(templates00T),len(templates00T)))
#plt.pcolormesh(X,Y,Z)
plt.show()
