import h5py
import math
import matplotlib.pyplot as plt

def get_derive_time(var, times, i_tstep) :
    deriv = (var[i_tstep+1] - var[i_tstep-1])/(times[i_tstep+1] - times[i_tstep-1])
    return deriv

h5file = h5py.File("Kepler.h5","r")
h5group = h5file['Orbit']

times_num = h5group["time"][:]
xA_num = h5group["xA"]
yA_num = h5group["yA"]
xB_num = h5group["xB"]
yB_num = h5group["yB"]

vAx_num = h5group["vAx"]
vAy_num = h5group["vAy"]
vBx_num = h5group["vBx"]
vBy_num = h5group["vBy"]

# verify velocities by time-derivative of positions (use centred
# second-order estimate)
n_tstep = len(vAx_num)
print("number of time steps:",n_tstep)

d_xA_num_dt = []
d_yA_num_dt = []
d_xB_num_dt = []
d_yB_num_dt = []
d_xA_num_dt.append(0)
d_yA_num_dt.append(0)
d_xB_num_dt.append(0)
d_yB_num_dt.append(0)

vA_ave = 0.
vB_ave = 0.
steps = 0.

for i_tstep in range(1,n_tstep-1,1) :
    steps += 1.
    d_xA_num_dt.append( get_derive_time(xA_num, times_num, i_tstep) )
    d_yA_num_dt.append( get_derive_time(yA_num, times_num, i_tstep) )
    vA_ave += math.sqrt(d_xA_num_dt[i_tstep]**2 + d_yA_num_dt[i_tstep]**2)

    d_xB_num_dt.append( get_derive_time(xB_num, times_num, i_tstep) )
    d_yB_num_dt.append( get_derive_time(yB_num, times_num, i_tstep) )
    vB_ave += math.sqrt(d_xB_num_dt[i_tstep]**2 + d_yB_num_dt[i_tstep]**2)

    #    print( i_tstep,d_xA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep] )

vA_ave /= steps    
vB_ave /= steps

eps = 1.e-2
for i_tstep in range(1,n_tstep-1,1) :
    if( abs(d_xA_num_dt[i_tstep] - vAx_num[i_tstep]) > abs(eps*vA_ave) ) :
        print(" Largest deviations for vAx at",i_tstep)
        print(abs(d_xA_num_dt[i_tstep] - vAx_num[i_tstep]),d_xA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_yA_num_dt[i_tstep] - vAy_num[i_tstep]) > abs(eps*vA_ave) ) :
        print(" Largest deviations for vAy at",i_tstep)
        print(abs(d_yA_num_dt[i_tstep] - vAy_num[i_tstep]),d_yA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_xB_num_dt[i_tstep] - vBx_num[i_tstep]) > abs(eps*vB_ave) ) :
        print(" Largest deviations for vBx at",i_tstep)
        print(abs(d_xB_num_dt[i_tstep] - vBx_num[i_tstep]),d_xB_num_dt[i_tstep],vBx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_yB_num_dt[i_tstep] - vBy_num[i_tstep]) > abs(eps*vB_ave) ) :
        print(" Largest deviations for vBy at",i_tstep)
        print(abs(d_yB_num_dt[i_tstep] - vBy_num[i_tstep]),d_yB_num_dt[i_tstep],vBx_num[i_tstep],vBy_num[i_tstep])  


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(xA_num, yA_num, color='r',linestyle='', marker='o', markersize=1)
ax.plot(xB_num, yB_num, color='g',linestyle='', marker='o', markersize=1)
#ax.plot(times_num, times_num**2, color='r',linestyle='', marker='o', markersize=1)
ax.set_aspect('equal')
figName = "Orbits.pdf"
fig.savefig(figName, transparent=True)
