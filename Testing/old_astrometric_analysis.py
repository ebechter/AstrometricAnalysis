'''
Calculate astrometric motion of star over given time frame using JPL Horizons ephermerides 
and SIMBAD proper motion and parallax 
'''
from __future__ import print_function
import matplotlib
matplotlib.use('agg')
import csv
import jplephem
import de421
from jplephem.spk import SPK
kernel = SPK.open('de430.bsp')
import numpy as np
import matplotlib.pyplot as pl
import astropy.units as u
from astropy.constants import G
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime
import matplotlib.ticker as ticker
import warnings
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
warnings.filterwarnings('ignore', category=UserWarning)
startTime = datetime.now()

def convert_coords(array):
	RADEC = SkyCoord(array[0],array[1],frame='icrs')
	RA = RADEC.ra.rad
	DEC = RADEC.dec.rad
	return RA, DEC

def proj_RA_DEC(RA, DEC, pm_ra,pm_dec,prlx,epoch,JD):

    """
    Takes the star's ra, dec, parallax, proper motion and a given epoch.

    Returns the parallax and proper motion shift in star's position, in milliarcseconds (sexagesimal seconds for raoff; decimal seconds for decoff)

    INPUTS:
    RA          -- RA position of star, in radians
    DEC         -- DEC position of star, in radians
    parallax    -- parallax of star, in MILLIarcseconds (1000.0*1/distance_pc)
    epoch       -- epoch (decimal years) to compute parallax offset (scalar or monotonically increasing vector)
    
    OUTPUTS:

    prop_DEC    -- Declination of each point along a time baseline     
    prop_RA     -- Right Ascension of each point along a time baseline
    prop_par_DEC -- parallax shift in star's position, in milliarcseconds
    prop_par_RA      -- parallax shift in star's position, in milliarcseconds

    """


    x,y,z = earthPos(JD) # get earth geocenter ephemerides using jplephem package


    prop_DEC = DEC_J2000 + pm_dec * epoch      # Calculate DEC first here
    if np.size(prop_DEC)>1:
        prop_RA  = RA  + pm_ra * epoch/np.cos(prop_DEC[0]) # rads
    else: 
        prop_RA  = RA  + pm_ra * epoch/np.cos(prop_DEC) # rads
    
    prop_par_RA  =  prop_RA +  prlx/np.cos(prop_DEC)*(x*np.sin(prop_RA) - y*np.cos(prop_RA))
    prop_par_DEC = prop_DEC +  prlx*(x*np.cos(prop_RA)*np.sin(prop_DEC) + y*np.sin(prop_RA)*np.sin(prop_DEC) - z* np.cos(prop_RA))
    
    return prop_DEC,prop_RA,prop_par_RA,prop_par_DEC

def earthPos(JD):

    eph = jplephem.Ephemeris(de421)
    barycenter = eph.position('earthmoon', JD)/  1.49597870700e+8
    moonvector = eph.position('moon', JD) /  1.49597870700e+8
    earthPos = (barycenter - moonvector * eph.earth_share) 
    x = earthPos[0,:]
    y = earthPos[1,:]
    z = earthPos[2,:]
    return x,y,z

def mas2rad(mas):
    return mas*2*np.pi/1000/3600/360  # rads

def rad2mas(rad):
    return rad*360/2/np.pi*3600*1000

def Sec(theta):
    return 1/np.cos(theta)


# Set these values for each system
target = 'HAT-P-7'
params = 'TRENDS_ast_params.txt'
data_name = target+ '_ast.txt'
npoints =25            # number of points that makes up the "tornado path"


## Read in Data
labels = np.loadtxt(params, delimiter=',', dtype=np.str, usecols=[0])
ICRS   = np.loadtxt(params, delimiter=',', dtype=np.str,usecols=[1,2])
values = np.loadtxt(params, delimiter=',', dtype=np.float,usecols=[3,4,5,6,7,8])
data   = np.loadtxt(data_name, delimiter=',', dtype=np.float)
a = np.where(np.char.find(labels, target) > -1) # string comparison
ind = a[0][0]                  # save index and strip extra array things. not the cleanest code


# Convert input params into radians 
RA_J2000, DEC_J2000 = convert_coords(ICRS[ind,:]) # convert coordinates to radians
values = mas2rad(values)   # convert everything else to radians now

 
# Assign SIMBAD parameters individual names	
prlx    = values[0][0]   # mas
dprlx   = values[0][1]   # mas
pm_ra   = values[0][2]   # mas/yr
dpm_ra  = values[0][3]   # mas/yr
pm_dec  = values[0][4]   # mas/yr
dpm_dec = values[0][5]   # mas/yr

# Assign measured astrometry parameters individual names 
D0_NS   = data[0][3]  # NS initial positions  - mas
D0_EW   = data[0][1]  # EW initial positions  - mas
JD      = data[:,0]   # full JD
NS      = data[:,3]   # NS data
EW      = data[:,1]   # EW data
dEW     = data[:,2]   # EW error
dNS     = data[:,4]   # NS error


# Convert JD into decimal year
t = Time(JD,format = 'jd').decimalyear

# create time baseline for tornado path
# dt = np.max(t)-np.min(t)
# tl = np.linspace(np.min(t),np.max(t),npoints)
tl = np.linspace(np.floor(np.min(t)), np.ceil(np.max(t)), npoints)
JDL = Time(tl,format = 'decimalyear').jd


# Print out system parameters... dont really need this anymore
print("""Astrometric parameters in radians:

    Name : {0}
    RA_J2000  = {1}
    DEC_J2000 = {2}
    Prlx   = {3} +/- {4}
    pm_ra  = {5} +/- {6} 
    pm_dec = {7} +/- {8} 
    D_NS   = {9} 
    D_EW  = {10} 

""".format(labels[ind],RA_J2000,DEC_J2000,prlx,dprlx,pm_ra,dpm_ra,pm_dec,dpm_dec,D0_NS,D0_EW))


# position = kernel[0,3].compute(JD)/ 1.49597870700e+8
# print(position)
# Use jplephem with de421 to determine position of Earth geocenter wrt SS Barycenter
# earthPos   = np.loadtxt('HAT-P-7.txt', delimiter='\t', dtype=np.float)
# print(earthPos)
# print(earthPos)



# Combine proper motion and parallax with RA and DEC to project star path in time

# Calculate initial positions
x0,y0,z0=earthPos(JD[0])
_,_,RA_start,DEC_start = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra,pm_dec,prlx,t[0],JD[0])#,x0,y0,z0)



# Calculate full curves
x,y,z = earthPos(JDL)
prop_DEC,prop_RA,prop_par_RA,prop_par_DEC = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra,pm_dec,prlx,tl,JDL)#,x,y,z)

# pl.figure()
# pl.plot(prop_DEC,tl)
# pl.show()




# Error Propagation

dDEC = mas2rad(dNS[0])
dRA = mas2rad(dEW[0])

# RA Error
sigma_ra_2 = Sec(DEC_J2000 + pm_dec*tl)**2 *(dprlx**2 * (y*np.cos(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)) - x*np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)))**2 +\
dpm_ra**2*prlx**2*tl**2*Sec(DEC_J2000 + pm_dec*tl)**2 * (x*np.cos(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)) + y*np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)))**2 +\
dpm_dec**2*prlx**2*tl**2*(np.cos(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl))*(-y + pm_ra*tl*x*Sec(DEC_J2000 + pm_dec*tl))   +  (x + pm_ra*tl*y*Sec(DEC_J2000 + pm_dec*tl))*\
np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)))**2*np.tan(DEC_J2000 + pm_dec*tl)**2)+\
prlx**2*Sec(DEC_J2000 + pm_dec*tl)**2   *( dRA**2*(x*np.cos(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)) + y*np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)))**2 + dDEC**2*\
(np.cos(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl))*(-y + pm_ra*tl*x*Sec(DEC_J2000 + pm_dec*tl)) + (x + pm_ra*tl*y*Sec(DEC_J2000 + pm_dec*tl))*\
np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)))**2*np.tan(DEC_J2000 + pm_dec*tl)**2)

sigma_ra = np.sqrt(rad2mas(np.sqrt(sigma_ra_2))**2+np.min(dEW)**2)


# DEC Error

sigma_dec_2=dprlx**2 *(np.cos(RA_J2000 + pm_ra* tl* Sec(DEC_J2000 + pm_dec* tl))* (-z + x* np.sin(DEC_J2000 + pm_dec* tl)) + y *np.sin(DEC_J2000 + pm_dec* tl)\
*np.sin(RA_J2000 + pm_ra *tl *Sec(DEC_J2000 + pm_dec* tl)))**2 + dRA**2 *prlx**2 *(y *np.cos(RA_J2000 + pm_ra* tl *Sec(DEC_J2000 + pm_dec* tl))\
*np.sin(DEC_J2000 + pm_dec* tl) + (z* - x *np.sin(DEC_J2000 + pm_dec* tl))*np.sin(RA_J2000 + pm_ra *tl *Sec(DEC_J2000 + pm_dec* tl)))**2\
+dpm_ra**2 *prlx**2* tl**2* (z *Sec(DEC_J2000 + pm_dec* tl) *np.sin(RA_J2000 +pm_ra*tl *Sec(DEC_J2000 + pm_dec*tl)) + \
(y* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)))*np.tan(DEC_J2000 + pm_dec*tl))**2\
+dDEC**2*prlx**2*(np.cos(DEC_J2000 +pm_dec*tl)* (x *np.cos(RA_J2000 + pm_ra*tl *Sec(DEC_J2000 + pm_dec*tl)) +y* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)))\
+pm_ra*tl* np.tan(DEC_J2000 +pm_dec*tl)* (z* Sec(DEC_J2000 + pm_dec*tl)*np.sin(RA_J2000 + pm_ra*tl*Sec(DEC_J2000 + pm_dec*tl)) + \
(y* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl))) *np.tan(DEC_J2000 + pm_dec*tl)))**2 +\
dpm_dec**2* prlx**2* tl**2* (np.cos(DEC_J2000 +pm_dec*tl)* (x* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) +y *np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)))\
+pm_ra*tl *np.tan(DEC_J2000 + pm_dec*tl)* (z* Sec(DEC_J2000 + pm_dec*tl)* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) + (y* np.cos(RA_J2000 + pm_ra*tl\
*Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl))) *np.tan(DEC_J2000 + pm_dec*tl)))**2

sigma_dec = np.sqrt(rad2mas(np.sqrt(sigma_dec_2))**2+np.min(dNS)**2)

# Total change in RA and DEC in mas

dprop_par_RA_mas =  rad2mas(prop_par_RA - RA_start)
dprop_par_DEC_mas = rad2mas(prop_par_DEC - DEC_start)
# dprop_par_RA_mas =  ((prop_par_RA-prop_par_RA[0])*u.rad).to(u.arcsec).value *1000  
# dprop_par_DEC_mas = ((prop_par_DEC-prop_par_DEC[0])*u.rad).to(u.arcsec).value *1000 


# Convert RA and DEC to NS,EW
prop_par_EW_mas=-dprop_par_RA_mas*np.cos(prop_DEC)  # EW motion of a background object relative to star (mas)
prop_par_NS_mas=-dprop_par_DEC_mas                   # NS motion of a background object relative to star (mas)

# Start from initial NS,EW offset (start from the first data point)
NS_vector = D0_NS + prop_par_NS_mas
EW_vector = D0_EW + prop_par_EW_mas

hfont = {'fontname':'Helvetica'}

# plt.title('title',**csfont)
# plt.xlabel('xlabel', **hfont)
majorLocator   = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(0.2)

'''
trying to fix the fill between plots. 
'''

fig, axarr=pl.subplots(2,sharex=True)
axarr[1].xaxis.set_major_locator(majorLocator)
axarr[1].xaxis.set_major_formatter(majorFormatter)
axarr[1].xaxis.set_minor_locator(minorLocator)

axarr[0].plot(tl, NS_vector, color="black", lw=2, alpha=1)
# Plot North Offset 1 and 2 sigma errors
axarr[0].fill_between(tl, NS_vector, NS_vector+sigma_dec,alpha=0.3,color='black',linewidth=0)
axarr[0].fill_between(tl, NS_vector, NS_vector-sigma_dec,alpha=0.3,color='black',linewidth=0)
axarr[0].fill_between(tl, NS_vector, NS_vector+2*sigma_dec,alpha=0.2,color='black',linewidth=0)
axarr[0].fill_between(tl, NS_vector, NS_vector-2*sigma_dec,alpha=0.2,color='black',linewidth=0)
axarr[0].errorbar(t,NS, yerr=dNS, fmt=".k") 

axarr[0].text(0.8, 0.8,target , fontweight='bold', fontsize = 25, horizontalalignment='center', verticalalignment='center', transform=axarr[0].transAxes)

# pl.errorbar(t, data, yerr=rverr, fmt=".k")
axarr[0].set_ylabel("North offset (mas)",fontweight='bold',fontsize=16)
axarr[0].yaxis.set_major_locator(MaxNLocator(prune='both'))
axarr[0].locator_params(axis = 'y', nbins = 6)  
axarr[1].plot(tl, EW_vector, color="black", lw=2, alpha=1)

# Plot East Offset 1 and 2 sigma errors
axarr[1].fill_between(tl, EW_vector, EW_vector+sigma_ra,alpha=0.3,color='black',linewidth=0)
axarr[1].fill_between(tl, EW_vector, EW_vector-sigma_ra,alpha=0.3,color='black',linewidth=0)
axarr[1].fill_between(tl, EW_vector, EW_vector+2*sigma_ra,alpha=0.2,color='black',linewidth=0)
axarr[1].fill_between(tl, EW_vector, EW_vector-2*sigma_ra,alpha=0.2,color='black',linewidth=0)
axarr[1].errorbar(t,EW, yerr=dEW, fmt=".k") 
# axarr[1].plot(tl, EW_vector-sigma_ra, color="black", lw=1, alpha=0.9)
axarr[1].set_ylabel("East offset (mas)",fontweight='bold',fontsize=16)
axarr[1].set_xlabel("Epoch (yrs)",fontweight='bold',fontsize=16)
axarr[1].yaxis.set_major_locator(MaxNLocator(prune='both'))
# axarr[1].get_xaxis().get_major_formatter().set_useOffset(False)
axarr[1].locator_params(axis = 'y', nbins = 6) #(or axis = 'y') 
pl.xlim([np.min(tl),np.max(tl)])
fig.tight_layout()
fig.subplots_adjust(hspace=0.001) # no horizontal space between figures

pl.savefig(target + '_plot_NSEW.png')



# print(NS_vector,EW_vector)
# Plot PA and SEP now
# print(NS_vector**2,EW_vector**2)
PA =np.abs((180./np.pi)*np.arctan2(NS_vector,EW_vector))+90

SEP = np.sqrt(NS_vector**2+EW_vector**2)
# print(PA,SEP)


fig, axarr=pl.subplots(2,sharex=True)
axarr[1].xaxis.set_major_locator(majorLocator)
axarr[1].xaxis.set_major_formatter(majorFormatter)
axarr[1].xaxis.set_minor_locator(minorLocator)

axarr[0].plot(tl, PA, color="black", lw=2, alpha=1)
# Plot North Offset 1 and 2 sigma errors
# axarr[0].fill_between(tl, NS_vector, NS_vector+sigma_dec,alpha=0.3,color='black',linewidth=0)
# axarr[0].fill_between(tl, NS_vector, NS_vector-sigma_dec,alpha=0.3,color='black',linewidth=0)
# axarr[0].fill_between(tl, NS_vector, NS_vector+2*sigma_dec,alpha=0.2,color='black',linewidth=0)
# axarr[0].fill_between(tl, NS_vector, NS_vector-2*sigma_dec,alpha=0.2,color='black',linewidth=0)
# axarr[0].errorbar(t,NS, yerr=dNS, fmt=".k") 

axarr[0].text(0.8, 0.8,target , fontweight='bold', fontsize = 25, horizontalalignment='center', verticalalignment='center', transform=axarr[0].transAxes)

# pl.errorbar(t, data, yerr=rverr, fmt=".k")
axarr[0].set_ylabel("PA (rad)",fontweight='bold',fontsize=16)
axarr[0].yaxis.set_major_locator(MaxNLocator(prune='both'))
axarr[0].locator_params(axis = 'y', nbins = 6)  
axarr[1].plot(tl, SEP, color="black", lw=2, alpha=1)

# Plot East Offset 1 and 2 sigma errors
# axarr[1].fill_between(tl, EW_vector, EW_vector+sigma_ra,alpha=0.3,color='black',linewidth=0)
# axarr[1].fill_between(tl, EW_vector, EW_vector-sigma_ra,alpha=0.3,color='black',linewidth=0)
# axarr[1].fill_between(tl, EW_vector, EW_vector+2*sigma_ra,alpha=0.2,color='black',linewidth=0)
# axarr[1].fill_between(tl, EW_vector, EW_vector-2*sigma_ra,alpha=0.2,color='black',linewidth=0)
# axarr[1].errorbar(t,EW, yerr=dEW, fmt=".k") 
# axarr[1].plot(tl, EW_vector-sigma_ra, color="black", lw=1, alpha=0.9)
axarr[1].set_ylabel("Sep (mas)",fontweight='bold',fontsize=16)
axarr[1].set_xlabel("Epoch (yrs)",fontweight='bold',fontsize=16)
axarr[1].yaxis.set_major_locator(MaxNLocator(prune='both'))
# axarr[1].get_xaxis().get_major_formatter().set_useOffset(False)
axarr[1].locator_params(axis = 'y', nbins = 6) #(or axis = 'y') 
pl.xlim([np.min(tl),np.max(tl)])
fig.tight_layout()
fig.subplots_adjust(hspace=0.001) # no horizontal space between figures

pl.savefig(target + '_plot_PASEP.png')


pl.figure()
pl.plot(EW_vector,NS_vector,color='black')
pl.errorbar(EW[0],NS[0],yerr=dEW[0],xerr=dNS[0],fmt='ok',markersize=12)
pl.errorbar(EW[1],NS[1],yerr=dEW[1],xerr=dNS[1],fmt='dk',markersize=12)
# pl.errorbar(EW[2],NS[2],yerr=dEW[2],xerr=dNS[2],fmt='sk',markersize=12)
[]
pl.xlabel('East Offset (mas)')
pl.ylabel('North Offset (mas)')
pl.tight_layout()
pl.savefig('test.png')

