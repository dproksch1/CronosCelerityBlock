import math
import matplotlib.pyplot as plt

vals = []
res = []

vals.append(0.00360391)
vals.append(0.000893484)
vals.append(0.000221779)
vals.append(5.53181e-05)
vals.append(1.38143e-05)
vals.append(3.45807e-06)

res.append(16)
res.append(32)
res.append(64)
res.append(128)
res.append(256)
res.append(512)

order = []
number = len(res)

for ires in range(number-1) :
    order_loc = math.log(vals[ires]/vals[ires+1])/math.log(2)
    order.append(order_loc)


num_sec = 100
res_min = res[0]
res_max = res[len(res)-1]

ratio = res_max/res_min

err_res = []
err_val = []

for inum in range(num_sec) :
    res_loc = res_min*ratio**(inum/(num_sec-1.))
    err_loc = 0.5/res_loc**2.
    print str(res_loc)
    err_res.append(res_loc)
    err_val.append(err_loc)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(res,vals, color='g',marker="o",linestyle="");
ax.plot(err_res,err_val)
ax.set_yscale('log')
ax.set_xscale('log')
fig.savefig("order.pdf")



