import numpy as np
import sys, os
import math
#import time
#from operator import itemgetter
#import operator
import matplotlib as mpl
import pylab as plt
from matplotlib.pyplot import cm

#Fixed parameters-----------------------------------------
colors = cm.get_cmap('magma', 110)
d2r = math.pi/180.
r2d = 180./math.pi
#Sidereal rate is the rotation of the Earth, in radian per hr.
Rsid = (2*(math.pi) ) / (23.93447)

#Telescope specific, here set for CHIME/CHORD-------------
LAT = float(49.4911)
Long = 119.5886
Diff_Long = Long / 15. #In hour
Direction = "W"
Declimit = 20-90+LAT 
Bin = float(sys.argv[1])
FOV = Bin /  180. * (math.pi) #In radian

#Input file containing RA/DEC of sources-------------------
INF = "psrcat_radec.dat"
RA = np.genfromtxt(INF).transpose()[1]
DEC = np.genfromtxt(INF).transpose()[2]

LST_rise = []
LST_set = []
ALT = []
DEC2 = []

for j in range(len(RA)): 
    if ( DEC[j] > Declimit) : #Only consider those above horizon of 20 deg

        #Drift time is always the same for a particular pulsar
        try:
            DriftAngle = math.acos((np.cos(FOV) - (np.sin(DEC[j]*d2r))**2 ) / (np.cos(DEC[j]*d2r))**2 )
            TIME = DriftAngle / Rsid
        except ValueError:
            TIME = 24

        #Find the LST rannge (based on RA only)
        LST = RA[j]/15.
        r = LST - TIME/2.
        s = LST + TIME/2.

        #[TO-DO] Take care of some pulsars rising just before 24:00 and setting after        
        ALT.append(90 - LAT + DEC[j])
        DEC2.append(DEC[j])
        LST_rise.append(r)
        LST_set.append(s)

        #If circumpol, list it again after 12 hrs
        if ( DEC[j] >= (90-LAT)) :
            if (RA[j] <= 180):
                RAapp = RA[j] + 180
            else:
                RAapp = RA[j] - 180
            LST = RAapp/15.
            r = LST - TIME/2.
            s = LST + TIME/2.
            ALT.append(270 - LAT - DEC[j])
            DEC2.append(DEC[j])
            LST_rise.append(r)
            LST_set.append(s)

#Make the plot----------------------------------------------
mpl.rcParams.update({'font.size': 13,'font.family': 'serif'})
fig = plt.figure(figsize=[12,6])
ax = plt.subplot2grid((5,1),(2,0), rowspan=3,colspan=1)
plt.subplots_adjust(wspace = 0.9,left = 0.06, right=0.94)

for j in range(len(LST_rise)):
    ax.plot( (LST_rise[j],LST_set[j]),(ALT[j],ALT[j]),lw=2,c=colors(int(DEC2[j]+20)),zorder=6)

plt.axhline(y=180-LAT, c='k') #NCP
plt.axhline(y=180-LAT*2, c='k',ls=':') #lower limit of circumpolar region

ax.set_xlim(0,24)
ax.set_ylim([20,180])
ax.set_xlabel("LST (hr)")
ax.set_ylabel("Altitude (deg)")

axY2 = ax.twinx()
axY2.set_ylim([20-90+LAT,180-90+LAT])
axY2.set_yticks(np.arange(-20, 140, 20))
axY2.set_yticklabels(['-20','0','20','40','60','80','80','60'])

ax.set_xticks(np.arange(0, 24, 1))
minor_ticks = np.arange(0, 24, 1)
ax.set_xticks(minor_ticks, minor=True)

ax2 = plt.subplot2grid((5,1),(0,0), rowspan=2,colspan=1)
TimeOne = np.sort(LST_rise)

PerAlt = np.arange(20,180-LAT,Bin)
for d in range(len(PerAlt)):
    CountAll = []
    for i in range(len(TimeOne)):
        c = 0
        for j, (r, s, a) in enumerate(zip(LST_rise,LST_set,ALT)):
            if ((r <= TimeOne[i]) and (s >= TimeOne[i])) and (a<=PerAlt[d]+Bin/2. and a>PerAlt[d]-Bin/2.) :
                c = c+1
        CountAll.append(c)
    ax2.plot(TimeOne,CountAll,c=colors(int(PerAlt[d]-90+LAT+20)))

ax2.set_xlim(0,24)
ax2.set_xticks(np.arange(0, 24.1, 1))
ax2.set_xticks(minor_ticks, minor=True)


from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
ax2.yaxis.set_minor_locator(MultipleLocator(5))

ax2.set_xticklabels([])
ax2.set_ylabel("#psr")
ax2.grid()
ax2.set_title("CHORD only known pulsars: "+str(int(Bin))+" deg beam width")
ax2.grid(b=True, which='minor', color='gainsboro', linestyle='-',alpha=0.5)

plt.figtext(0.2,0.43,"Circumpolar",fontdict={'fontsize':13},horizontalalignment='center')
plt.figtext(0.2,0.4,"Region",fontdict={'fontsize':13},horizontalalignment='center')

plt.figtext(0.98,0.32,"Declination ($^{\circ}$)",rotation=-90,horizontalalignment='center',verticalalignment='center')

#plt.show()
plt.savefig("CHORD-known-beamwidth_"+str(int(Bin))+"deg.png",transparent=False, bbox_inches='tight')

