import matplotlib as mpl
import numpy as np
import pylab as plt
import ephem
import time
import sys

today=time.strftime("%x")

plt.rcParams.update({'font.size':14})
plt.figure(figsize=[10,5])
ax = plt.axes(projection="hammer")

#setup axes and tick labels
plt.text(0.5,-0.05,r"$l$",
         fontsize=16,
         verticalalignment="top",
         horizontalalignment="center",
         transform=ax.transAxes,
         alpha=1,weight='bold' )
plt.ylabel(r"$b$",fontsize=16,rotation=0,
           alpha=1)
locs,labs = plt.xticks()
labs = [r"$%s^{\circ}$"%(-1*ii) for ii in range(-150,210,30)]
plt.xticks(locs,labs)
locs,labs = plt.yticks()
labs = [r"$%s^{\circ}$"%(ii) for ii in range(-75,90,15)]
plt.yticks(locs,labs)
plt.grid()

#transform 0-360 gal to -pi to pi
def trans(val):
    if val<180:
        return -1*np.pi*val/180.
    else:
        return (360-val)*np.pi/180.

def transrad(val):
    print val
    if val > np.pi:
        return -1*(2*np.pi-val)
    else:
        return val

x = np.genfromtxt("known-pulsars.dat").transpose()
x[0] = np.array([trans(i) for i in x[0]])
x[1]*=(np.pi/180.)
plt.scatter(x[0],x[1],s=10,marker="x",c="Orange",zorder=1000,edgecolor="OrangeRed",alpha=0.7, label="Previously known pulsars")

htruN = np.genfromtxt("new-HTRU-North.dat", dtype={'names': ('ra', 'dec'),'formats': ('S12', 'S12')}, unpack=True).transpose()
pos =np.array([(j.lon,j.lat) for j in [ephem.Galactic(ephem.Equatorial(i[0],i[1])) for i in htruN]]).transpose()
nx = np.array([transrad(i) for i in pos[0]])
plt.scatter(-nx,pos[1],marker="s",s=45,color="y",edgecolor="k",zorder=2000,label="HTRU-North")

htruNm = np.genfromtxt("new-HTRU-North-MSP.dat", dtype={'names': ('ra', 'dec'),'formats': ('S12', 'S12')}, unpack=True).transpose()
pos =np.array([(j.lon,j.lat) for j in [ephem.Galactic(ephem.Equatorial(i[0],i[1])) for i in htruNm]]).transpose()
nx = np.array([transrad(i) for i in pos[0]])
plt.scatter(-nx,pos[1],marker="s",s=45,color="black",edgecolor="k",zorder=2100,label="HTRU-North MSPs")

htruS = np.genfromtxt("new-HTRU-South.dat", dtype={'names': ('ra', 'dec'),'formats': ('S12', 'S12')}, unpack=True).transpose()
pos =np.array([(j.lon,j.lat) for j in [ephem.Galactic(ephem.Equatorial(i[0],i[1])) for i in htruS]]).transpose()
nx = np.array([transrad(i) for i in pos[0]])
plt.scatter(-nx,pos[1],marker=".",s=150,color="red",edgecolor="k",zorder=2000, label="HTRU-South")

htruSm = np.genfromtxt("new-HTRU-South-MSP.dat", dtype={'names': ('ra', 'dec'),'formats': ('S12', 'S12')}, unpack=True).transpose()
pos =np.array([(j.lon,j.lat) for j in [ephem.Galactic(ephem.Equatorial(i[0],i[1])) for i in htruSm]]).transpose()
nx = np.array([transrad(i) for i in pos[0]])
plt.scatter(-nx,pos[1],marker=".",s=150,color="black",edgecolor="k",zorder=2100, label="HTRU-South MSPs")

#Plot FRBs (optional)
if len(sys.argv) > 1:
    frbFile = sys.argv[1]
    frb = np.genfromtxt(frbFile, dtype={'names': ('ra', 'dec'),'formats': ('S12', 'S12')}, unpack=True).transpose()
    pos =np.array([(j.lon,j.lat) for j in [ephem.Galactic(ephem.Equatorial(i[0],i[1])) for i in frb]]).transpose()
    nx = np.array([transrad(i) for i in pos[0]])
    plt.scatter(-nx,pos[1],marker="*",s=150,color="red",edgecolor="k",zorder=2000, label="FRBs")


plt.legend(bbox_to_anchor=(1.105,1.08),  scatterpoints = 1, prop={'size':12})
plt.figtext(0.38, 0.95, "HTRU discoveries (as of "+str(today)+")")
plt.show()

