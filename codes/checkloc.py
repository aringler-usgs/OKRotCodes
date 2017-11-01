#!/usr/bin/env python
from obspy.core import UTCDateTime, read, Stream
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

import numpy as np

##########################################


import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


stime = UTCDateTime('2017-117T00:00:00.0')
etime = UTCDateTime('2017-158T00:00:00.0')
client = Client("IRIS")


stalat = 36.47819
stalon = -98.742240

# Map of events
cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=2.)
                        
for eve in cat:
    print(eve.origins[0])
