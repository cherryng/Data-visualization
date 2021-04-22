import numpy as np
import math
import time
import matplotlib as mpl
import pylab as plt
import datetime
from datetime import date
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

mpl.rcParams.update({'font.size': 18}) 

D2R = math.pi/180.
R2D = 180./math.pi
TAU = 2*math.pi

time_bin = 0.05 #Step every 3min
times = np.arange(0,24.1,time_bin)

inst_lat_PKS = -(32+59/60.+52/3600.)
inst_long_PKS = -(148+15/60.+47/3600.)

inst_lat_MK = -(30+43/60.+16/3600.)
inst_long_MK = -(21+24/60.+40/3600.)

inst_lat_VLA =34+4/60.+43/3600.
inst_long_VLA = 107+37/60.+4/3600.

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
decs = np.arange(90,-91,-1)

vis_hours = np.zeros(len(decs))
for j in range(len(decs)):
    dec = decs[j]
    alts_PKS, azs_PKS = find_track(inst_lat_PKS, inst_long_PKS, ra, dec)
    alts_MK, azs_MK = find_track(inst_lat_MK, inst_long_MK, ra, dec)
    alts_VLA, azs_VLA = find_track(inst_lat_VLA, inst_long_VLA, ra, dec)
    for i in range(len(times)):
        if alts_PKS[i] > 30 and alts_MK[i] > 20:
        #if alts_PKS[i] > 30 and alts_VLA[i] > 8:
            vis_hours[j] = vis_hours[j] + 1
    vis_hours[j] = vis_hours[j]*time_bin
            
ax.plot(decs, vis_hours, c='k')
plt.fill_between(decs,vis_hours,color='gainsboro')

ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_xlabel('Declination (deg)')
ax.set_ylabel('Overlap hours (hr)')

DEC_plot =  -62-40/60.-46/3600.
ax.axvline(x=DEC_plot, c='r', label="Prox Cen",zorder=10)
DEC_plot = -29-28.1/3600.
ax.axvline(x=DEC_plot, c='b', label="Galactic Centre",zorder=10)
DEC_plot = 38+47/60.+1.2802/3600.
#ax.axvline(x=DEC_plot, c='g', label="Vega")

legend = plt.legend(loc=0,fontsize=15)
frame = legend.get_frame()
frame.set_facecolor('white')
plt.grid()
plt.grid(b=True, which='major', color='gray', linestyle='-')
plt.grid(b=True, which='minor', color='gainsboro', linestyle='-',alpha=0.5)

#ax.set_xlim(-90,90)
ax.set_xlim(-90,0)

plt.title("Overlap hours b/w MeerKAT and Parkes")
#plt.title("Overlap hours b/w VLA and Parkes")
#plt.show()
plt.savefig("Overlap-sky-MK-PK.png",bbox_inches='tight')
#plt.savefig("Overlap-sky-VLA-PKS.png",bbox_inches='tight')
