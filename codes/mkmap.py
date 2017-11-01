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



fig = plt.figure(1, figsize=(8, 10))
m = Basemap(projection='merc', lon_0=stalon, lat_0=stalat, llcrnrlon=stalon-5, 
            llcrnrlat=stalat-2, urcrnrlon=stalon+3, urcrnrlat=stalat+10, epsg=4269, resolution="i")


#m.drawcoastlines()
m.fillcontinents(color='.9')
#m.drawcounties()
m.drawparallels([36.478], labels=[True, False, True, False])
m.drawmeridians([-98.74],labels=[True, False, True, False])
m.drawmapboundary()
m.drawstates(linewidth=1., zorder=3.)
ax=plt.gca()


lats=[]
lons=[]
mags=[]
for event in cat:
    lats.append(event.origins[0].latitude)
    lons.append(event.origins[0].longitude)
    mags.append(event.magnitudes[0].mag)

min_size = 2
max_size = 30
min_size_ = min(mags) - 1
max_size_ = max(mags) + 1
frac = [(0.2 + (_i - min_size_)) / (max_size_ - min_size_) for _i in mags]
size_plot = [(_i * (max_size - min_size)) ** 2 for _i in frac]

mags2 = [2., 2.5, 3., 3.5, 4., 4.5]
frac2 = [(0.2 + (_i - min_size_)) / (max_size_ - min_size_) for _i in mags2]
size_plot2 = [(_i * (max_size - min_size)) ** 2 for _i in frac2]

mymap=plt.get_cmap('viridis')
mycolor= mymap(frac2)
x,y = m(stalon-1.,stalat-1.5)
x2,y2 = m(stalon+2.5,stalat+1.5)
from matplotlib.patches import Polygon
p = Polygon([(x,y),(x,y2),(x2,y2),(x2,y)], facecolor='w',edgecolor='k') 
plt.gca().add_patch(p) 
x,y =m(lons,lats)
sc= m.scatter(x, y, s=size_plot, c=mags, zorder =3, cmap='viridis')
print(mags)
#cb=plt.colorbar(sc, orientation="horizontal",fraction=0.046, pad=0.04)
#cb.set_label('Event Magnitude ($M_L$)')
for mag, siz, f2,cc in zip(mags2,size_plot2, frac2, mycolor):
    sc2 = m.scatter([0.] ,[0.],s=siz ,c=cc,  cmap='viridis', zorder=3, label='$M_L$' + str(mag))
plt.legend(loc=9, ncol=3, bbox_to_anchor=(0.5,-0.01))

#for mag in mags2:
#    sc2 = m.scatter(0.,0.,mag,label='$M_L$' + str(int(mag)))
#plt.legend()

#cb.set_ticks([2.,3.,4.,5.])
#cb.set_ticklabels([2.,3.,4.,5.])
x,y =m(stalon, stalat)

m.scatter(x,y, 200, color="r", edgecolor="r", marker="v", zorder=3)
x,y =m(stalon+.15, stalat-.05)
#plt.text(x,y,'OKR1', va="top", family="monospace", weight="bold", color="r")


from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
axins = zoomed_inset_axes(ax, 2, loc=2) # zoom = 4
x,y =m(lons,lats)
sc= m.scatter(x, y, s=size_plot, c=mags, zorder =3)

x,y = m(stalon-1.,stalat-1.5)
x2,y2 = m(stalon+2.5,stalat+1.5)
axins.set_xlim(x,x2)
axins.set_ylim(y,y2)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax, axins, loc1=1, loc2=3, fc='r', ec="k")
x,y =m(stalon, stalat)
m.scatter(x,y, 200, color="r", edgecolor="r", marker="v", zorder=3)
x,y =m(stalon+.1, stalat-.2)
plt.text(x,y,'OKR1', va="top", family="monospace", weight="bold", color="r")
m.drawcounties()


m.drawmapboundary()
m.drawstates(linewidth=1., zorder=3.)














plt.savefig('Map.jpg',format='jpeg',dpi=400)
plt.savefig('Map.pdf',format='pdf',dpi=400)

plt.show()











#stime = UTCDateTime('2017-117T00:00:00.0')
#client = Client("IRIS")



## Map of events
#cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat, 
                        #longitude=stalon, maxradius=2.)

#inventory = client.get_stations(network="GS", station="OKR1")

#figmap = cat.plot(projection="local", resolution="i", show=False, title='Color by Depth (km)', label=None)
#figmap = inventory.plot(projection="local", show=False, fig=figmap, size=350, color="r", legend=None)
#figmap.show()
#figmap.savefig('Eventsplot.jpg',format='jpeg')
#figmap.clf()
