"""
Extract max By amplitude and position of maximum
from all frames of a given project
Creates optional plots of both as a function of time
Console output may be further processed in gnuplot.

Overplotted curves are linear and exponential, respectively,
based on first and last data point only, i.e. disregarding all
intermediate values!

Usage:  python wavefit.py <pname> [-showplot]'
"""

import numpy as np
import os
import re
import sys
import h5py
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def string_has_patt (string, patt):    
    regexp = re.compile (patt)
    number = len (regexp.findall(string))
    return (number == 1)

datadir = 'data'       # <-- adapt path to data if necessary

flag_pl = False
if len(sys.argv) == 3:
    if sys.argv[2] == '-showplot':
        flag_pl = True
    else:
        print ('Warning: argument '+ sys.argv[2] +' not recognized.')
else:
    if len(sys.argv) != 2:
        print ('Usage: python wavefit.py <pname> [-showplot]')
        exit()
pname = sys.argv[1]

floatfiles = os.listdir (datadir +'/'+ pname +'_float')

for filename in floatfiles:
    if (string_has_patt (filename, pname +'_flt_step.*\.h5\Z')) == False:
        floatfiles.remove (filename)

n_floats = len (floatfiles)
if n_floats == 0:
    print ('ERROR: No float files found.')
    exit(1)
else:
    pass
    #print (str (n_floats) +' h5 files found.')

floatfiles.sort()
dir = 0  # {0,1,2} = {x,y,z}

quan = 'B_z'
results = []
profiles = []
for frame in floatfiles:
    #print ('opening '+ frame +'...')
    h5file = h5py.File (datadir +'/'+ pname +'_float/'+ frame, 'r')
    h5group = h5file ['fluid0']
    dset = h5group [quan]
    xb = dset.attrs ["origin"][dir]
    dx = dset.attrs ["delta"][dir]
    n_cells = dset.shape[2-dir]
    
    dataset = h5group [quan]
    xVal = dataset [0,0,:]
    xPos = np.arange (0.5*dx, n_cells*dx, dx, dtype='f')

    print (len(xVal), len(xPos))
    h5group = h5file ['Data']
    time =  float (h5group.attrs['time'])
    mpos = xPos[xVal.argmax()]
    ampl = abs (xVal.max())
    results.append([time, mpos, ampl])
    profiles.append(xVal)

results.sort()
Ltime = []
Lmpos = []
Lampl = []

k_fr = 2.*math.pi/(dx*n_cells)     # estimated k value
k_rd = int (100*k_fr +0.5) / 100.  # round to within 0.01
lamb = 2.*math.pi / k_rd           # wavelength

for i in xrange(len(results)-1):
    while (results[i+1])[1] < (results[i])[1]:  # drop in t range?
        (results[i+1])[1] += lamb               # lift

for re in results:
    Ltime.append(re[0])
    Lmpos.append(re[1])
    Lampl.append(re[2])

vph = (Lmpos[0]-Lmpos[-1]) / (Ltime[0]-Ltime[-1])
xoo =  Lmpos[0]

damp = -math.log(Lampl[-1]/Lampl[0]) / (Ltime[-1]-Ltime[0])

print ('#  '+ pname)
print ('##  k='+ str(k_rd) +'  v_phase=' + str(vph) +'  damping='+ str(damp))

tran = np.arange(Ltime[0], Ltime[-1], (Ltime[-1]-Ltime[0])/200.)

print '# time     position  '+ quan +' amplitude'
print '# ----------------------------'
for data in results:
    print ' {:f}  {:f}  {:f}'.format(data[0], data[1], data[2])

if (flag_pl):  # -> show plot
    fig = plt.figure (tight_layout=True)
    fig.suptitle ("Debug plots for test run '" + pname + "'")
    gs = gridspec.GridSpec (2, 2)
    
    ax_prof = fig.add_subplot (gs[0, :])
    ax_prof.set_ylabel (quan)
    ax_prof.set_xlabel ('x')
    ax_prof.set_xlim ([xPos[0], xPos[-1]])
    for prof in profiles:
        ax_prof.plot (xPos, prof)

    ax_xpos = fig.add_subplot (gs[1, 0])
    ax_damp = fig.add_subplot (gs[1, 1])
    ax_xpos.set (xlabel='time', ylabel='pos. of max '+quan)
    ax_damp.set (xlabel='time', ylabel=quan+' amplitude')
    
    ax_xpos.plot (Ltime, Lmpos, 'r.')
    ax_damp.plot (Ltime, Lampl, 'r.')

    ax_xpos.plot ([Ltime[0], Ltime[-1]], [Lmpos[0], Lmpos[-1]], 'b')
    ax_damp.plot (tran, Lampl[0]*np.exp(-damp*tran), 'b')

    ax_xpos.set_xlim ([Ltime[0], Ltime[-1]])
    ax_damp.set_xlim ([Ltime[0], Ltime[-1]])
 
    plt.show()
