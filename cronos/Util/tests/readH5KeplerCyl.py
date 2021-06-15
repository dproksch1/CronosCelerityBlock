import h5py
import math
import matplotlib.pyplot as plt

def get_derive_time(var, times, i_tstep) :
    deriv = (var[i_tstep+1] - var[i_tstep-1])/(times[i_tstep+1] - times[i_tstep-1])
    return deriv

h5file = h5py.File("KeplerCyl.h5","r")
h5group = h5file['Orbit']

times_num = h5group["time"][:]
rhoA_num = h5group["rhoA"]
phiA_num = h5group["phiA"]
rhoB_num = h5group["rhoB"]
phiB_num = h5group["phiB"]

vAx_num = h5group["vAx"]
vAy_num = h5group["vAy"]
vBx_num = h5group["vBx"]
vBy_num = h5group["vBy"]

# verify velocities by time-derivative of positions (use centred
# second-order estimate)
n_tstep = len(vAx_num)
print("number of time steps:",n_tstep)

d_rhoA_num_dt = []
d_phiA_num_dt = []
d_rhoB_num_dt = []
d_phiB_num_dt = []
d_rhoA_num_dt.append(0)
d_phiA_num_dt.append(0)
d_rhoB_num_dt.append(0)
d_phiB_num_dt.append(0)

vA_ave = 0.
vB_ave = 0.
steps = 0.

for i_tstep in range(1,n_tstep-1,1) :
    steps += 1.
    d_rhoA_num_dt.append( get_derive_time(rhoA_num, times_num, i_tstep) )
    d_phiA_num_dt.append( get_derive_time(phiA_num, times_num, i_tstep) )
    vA_ave += math.sqrt(d_rhoA_num_dt[i_tstep]**2 + d_phiA_num_dt[i_tstep]**2)

    d_rhoB_num_dt.append( get_derive_time(rhoB_num, times_num, i_tstep) )
    d_phiB_num_dt.append( get_derive_time(phiB_num, times_num, i_tstep) )
    vB_ave += math.sqrt(d_rhoB_num_dt[i_tstep]**2 + d_phiB_num_dt[i_tstep]**2)

    #    print( i_tstep,d_rhoA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep] )

vA_ave /= steps    
vB_ave /= steps

eps = 1.e-2
for i_tstep in range(1,n_tstep-1,1) :
    if( abs(d_rhoA_num_dt[i_tstep] - vAx_num[i_tstep]) > abs(eps*vA_ave) ) :
        print(" Largest deviations for vAx at",i_tstep)
        print(abs(d_rhoA_num_dt[i_tstep] - vAx_num[i_tstep]),d_rhoA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_phiA_num_dt[i_tstep] - vAy_num[i_tstep]) > abs(eps*vA_ave) ) :
        print(" Largest deviations for vAy at",i_tstep)
        print(abs(d_phiA_num_dt[i_tstep] - vAy_num[i_tstep]),d_phiA_num_dt[i_tstep],vAx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_rhoB_num_dt[i_tstep] - vBx_num[i_tstep]) > abs(eps*vB_ave) ) :
        print(" Largest deviations for vBx at",i_tstep)
        print(abs(d_rhoB_num_dt[i_tstep] - vBx_num[i_tstep]),d_rhoB_num_dt[i_tstep],vBx_num[i_tstep],vAy_num[i_tstep])
    if( abs(d_phiB_num_dt[i_tstep] - vBy_num[i_tstep]) > abs(eps*vB_ave) ) :
        print(" Largest deviations for vBy at",i_tstep)
        print(abs(d_phiB_num_dt[i_tstep] - vBy_num[i_tstep]),d_phiB_num_dt[i_tstep],vBx_num[i_tstep],vBy_num[i_tstep])  


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.polar(phiA_num, rhoA_num, color='r',linestyle='', marker='o', markersize=1)
plt.polar(phiB_num, rhoB_num, color='g',linestyle='', marker='o', markersize=1)
#ax.plot(times_num, times_num**2, color='r',linestyle='', marker='o', markersize=1)
ax.set_aspect('equal')
figName = "OrbitsCyl.pdf"
fig.savefig(figName, transparent=True)
