#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

f = open('TemplateResults_00_OK038_0_118','r')
dataTime=False
vals =[]
for line in f:
    if dataTime:
        val = float(line)

        vals.append(val)
 
    if 'data' in line:
        dataTime= True
        
f.close()

f = open('TemplateResults_00_OKR1_0_118','r')
dataTime=False
vals2 =[]
for line in f:
    if dataTime:
        val = float(line)
        
        vals2.append(val)
      
    if 'data' in line:
        dataTime= True
        
f.close()

f = open('TemplateResults_10_OKR1_0_118','r')
dataTime=False
vals3 =[]
for line in f:
    if dataTime:
        val = float(line)
        
        vals3.append(val)
      
    if 'data' in line:
        dataTime= True
f.close()

fig = plt.figure(1)
plt.plot(vals)
plt.plot(vals2)
plt.plot(vals3)
plt.ylim((0.,1.))

plt.show()
