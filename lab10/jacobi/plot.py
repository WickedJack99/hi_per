import pylab as plt


ser_init = plt.loadtxt("serial_init.asc")
par_init = plt.loadtxt("parallel_init.asc")
ser_data = plt.loadtxt("result_serial.asc")
par_data = plt.loadtxt("result_parallel.asc")
# diff = par_data - ser_data

plt.figure()
plt.subplot(221)
plt.title("Serial initial condition")
plt.contourf(ser_init, cmap="jet", levels=100)
plt.subplot(222)
plt.title("Parallel initial condition")
plt.contourf(par_init, cmap="jet", levels=100)
plt.subplot(223)
plt.title("Serial result")
plt.contourf(ser_data, cmap="jet", levels=100)
plt.subplot(224)
plt.title("Parallel result")
plt.contourf(par_data, cmap="jet", levels=100)

plt.savefig("result.png")
