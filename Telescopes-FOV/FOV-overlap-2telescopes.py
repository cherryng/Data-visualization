import numpy as np
import math
import time
import matplotlib as mpl
import pylab as plt
import datetime
from datetime import date
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

mpl.rcParams.update({'font.size': 20})  #,'font.family': 'serif'})

D2R = math.pi/180.
R2D = 180./math.pi
TAU = 2*math.pi

time_bin = 0.05 #3min
times = np.arange(0,24.1,time_bin)

inst_lat_PKS = -(32+59/60.+52/3600.)
inst_long_PKS = -(148+15/60.+47/3600.)

inst_lat_MK = -(30+43/60.+16/3600.)
inst_long_MK = -(21+24/60.+40/3600.)

def find_track(inst_lat, inst_long, ra, dec):
    alts=[]
    azs=[]
    for time in times:
        LST = time + inst_long / 15.
        while (LST < 0):
            LST = LST + 24
        LST = np.mod(LST, 24)
        hour_angle = LST * 15. - ra
        alt = np.sin(dec*D2R)*np.sin(inst_lat*D2R)+np.cos(dec*D2R)*np.cos(inst_lat*D2R)*np.cos(hour_angle*D2R)
        alt = math.asin(alt)
        alts.append(alt*R2D)
        az = (np.sin(dec*D2R) - np.sin(alt)*np.sin(inst_lat*D2R))/(np.cos(alt)*np.cos(inst_lat*D2R))
        az = math.acos(az)
        if(np.sin(hour_angle*D2R) >= 0):
            az = TAU - az
        azs.append(az*R2D)
    return(alts, azs)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1)

ra = 0
decs = np.arange(0,-91,-1)

vis_hours = np.zeros(len(decs))

for j in range(len(decs)):
    dec = decs[j]
    alts_PKS, azs_PKS = find_track(inst_lat_PKS, inst_long_PKS, ra, dec)
    alts_MK, azs_MK = find_track(inst_lat_MK, inst_long_MK, ra, dec)
    for i in range(len(times)):
        if alts_PKS[i] > 30 and alts_MK[i] > 20:
            vis_hours[j] = vis_hours[j] + 1
    vis_hours[j] = vis_hours[j]*time_bin
            
ax.plot(decs, vis_hours, c='k')
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlabel('Declination (deg)')
ax.set_ylabel('Overlap hours (hr)')

DEC_plot =  -62-40/60.-46/3600.
ax.axvline(x=DEC_plot, c='b', label="Prox Cen")
DEC_plot = -29-28.1/3600.
ax.axvline(x=DEC_plot, c='r', label="Galactic Centre")

plt.legend(loc=0)
plt.grid()
plt.grid(b=True, which='major', color='gray', linestyle='-')
plt.grid(b=True, which='minor', color='gainsboro', linestyle='-',alpha=0.5)


plt.title("Overlap hours b/w MeerKAT and Parkes")
plt.show()
#plt.savefig("Overlap-sky-MK-PK.png")
