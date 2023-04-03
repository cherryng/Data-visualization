import matplotlib as mpl
import pylab as plt
import numpy as np


dms = np.arange(0,100)

mpl.rcParams.update({'font.size': 16,'font.family': 'serif'})

fig = plt.figure(figsize=[6,6])
ax = fig.add_subplot(1, 1, 1)


f1 = 0.045
f2 = 0.8
delay = 4.15*dms*(f1**-2-f2**-2)/1000.
ax.plot(dms,delay/60.,label='NenuFAR (45-80MHz)')

f1 = 0.11
f2 = 0.188
delay = 4.15*dms*(f1**-2-f2**-2)/1000.
ax.plot(dms,delay/60.,label='LOFAR')
print('nenufar max delay', np.max(delay))
f1 = 0.4
f2 = 0.8
delay = 4.15*dms*(f1**-2-f2**-2)/1000.
ax.plot(dms,delay/60.,label='CHIME')


f1 = 0.01
f2 = 0.8
delay = 4.15*dms*(f1**-2-f2**-2)/1000.
ax.plot(dms,delay/60.,label='NenuFAR (10-80MHz)')


plt.legend(loc=0)
plt.grid()
plt.xlabel('Dispersion Measure (cm-3pc)')
plt.ylabel('Time delay b/w top and bottom of band (min)')
#plt.show()
plt.savefig('Dispersive-delay.png',bbox_inches='tight')
