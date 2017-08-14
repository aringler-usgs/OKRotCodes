#!/usr/bin/env python
import io, urllib
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial import distance
url = "https://examples.obspy.org/dissimilarities.npz"
with io.BytesIO(urllib.urlopen(url).read()) as fh, np.load(fh) as data:
    dissimilarity = data['dissimilarity']
print(dissimilarity.size)
