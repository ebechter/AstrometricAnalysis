'''
Sample proper motion, parallax, and initial offset from a gaussian 
feed into main code to generate tracks
use percentile to take 68 and 95 positions on each timestep. 

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


def proj_RA_DEC(RA, DEC, pm_ra, pm_dec, prlx, epoch, JD, D0_NS, D0_EW,date1):
	'''
    RA_start    -- RA position of star, in radians
    DEC_start   -- DEC position of star, in radians 
    date1 		-- JD of first data point
	'''
	# first calculate projections due to proper motion only
	prop_RA,prop_DEC = pm_offset(RA,DEC,pm_ra,pm_dec,epoch)
	# then calculate projections due to parallax only 
	par_RA,par_DEC = parallax_offset(prlx, prop_RA,prop_DEC,JD)	

	# Add together
	delta_RA = prop_RA+par_RA
	delta_DEC = prop_DEC + par_DEC
	# get starting positions in RA and DEC
	date1_jd = Time(date1,format = 'decimalyear').jd
	RA_start, DEC_start = get_initial_RA_DEC(RA, DEC, pm_ra,pm_dec,prlx,date1,date1_jd)

	# Subtract starting position
	delta_RA =  rad2mas(delta_RA  - RA_start)
	delta_DEC = rad2mas(delta_DEC - DEC_start)

	# Convert RA and DEC to NS,EW
	delta_EW= -delta_RA*np.cos(prop_DEC)  # EW motion of a background object relative to star (mas)
	delta_NS= -delta_DEC                   # NS motion of a background object relative to star (mas)

	# Start from initial NS,EW offset (start from the first data point)
	NS_vector = D0_NS + delta_NS
	EW_vector = D0_EW + delta_EW
	return NS_vector,EW_vector


def pm_offset(RA,DEC,pm_ra,pm_dec,epoch):
	'''
	Calculates the simple RA and DEC offset for each point in epoch vector

	INPUT:
    RA          -- RA position of star, in radians
    DEC         -- DEC position of star, in radians 
    pm_ra 		-- proper motion in RA direction rad/yr
    pm_dec 		-- proper motion in DEC direction rad/yr
    epoch 	    -- decimalyear vector or single value
	
	OUTPUT: 
	prop_RA 	-- resulting offset due to RA proper motion (X(t)= X0 + V*t)
	prop_DEC 	-- resulting offset due to DEC proper motion 
	'''
	prop_DEC = DEC + pm_dec * epoch      # Calculate DEC first here

	if np.size(prop_DEC) > 1:
		prop_RA  = RA  + pm_ra * epoch/np.cos(prop_DEC[0])  # rads
	else: 

		prop_RA  = RA  + pm_ra * epoch/np.cos(prop_DEC) 	# rads
	return prop_RA, prop_DEC


def parallax_offset(prlx, prop_RA,prop_DEC,JD):
	
	'''
	Calculates the parallactic offset in RA and DEC directions

	INPUT:
    prlx    	-- parallax of star, in MILLIarcseconds (1000.0*1/distance_pc)

    JD 	    	-- full Julian date 
	
	OUTPUT: 
	par_RA
	par_DEC

	'''
	x,y,z = earthPos(JD) # get earth geocenter ephemerides using jplephem package

	par_RA  = prlx/np.cos(prop_DEC)*(x*np.sin(prop_RA) - y*np.cos(prop_RA))
	par_DEC = prlx*(x*np.cos(prop_RA)*np.sin(prop_DEC) + y*np.sin(prop_RA)*np.sin(prop_DEC) - z* np.cos(prop_RA))
    	
	return par_RA,par_DEC

def get_initial_RA_DEC(RA, DEC, pm_ra,pm_dec,prlx,epoch,JD):
	'''
	This is a bit redundant as I need to repeat proj_RA_DEC function to a certain point to get the correct
	zero (starting) point. This is the same calculation as the whole track but for just our first data point
	'''

	# first calculate projections due to proper motion only
	prop_RA,prop_DEC = pm_offset(RA,DEC,pm_ra,pm_dec,epoch)

	# then calculate projections due to parallax only 
	par_RA,par_DEC = parallax_offset(prlx, prop_RA,prop_DEC,JD)	

	# Add together
	RA_start = prop_RA + par_RA
	DEC_start = prop_DEC + par_DEC
	return RA_start, DEC_start


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

target = 'HAT-P-7'
params = 'TRENDS_ast_params.txt'
data_name = target+ '_ast.txt'
npoints =100            # number of points that makes up the "tornado path"
ntracks = 2000

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
tl = np.linspace(np.floor(np.min(t)), np.ceil(np.max(t)), npoints)

JDL = Time(tl,format = 'decimalyear').jd

NS_vector,EW_vector = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra, pm_dec, prlx, tl, JDL, D0_NS, D0_EW,t[0])




# Computing the errors - North first
Nstart_track = np.random.normal(D0_NS , scale = np.min(dNS), size = ntracks)
Estart_track = np.random.normal(D0_EW , scale = np.min(dEW), size = ntracks)
prlx_track   = np.random.normal(prlx  , scale = dprlx      , size = ntracks)
pm_dec_track = np.random.normal(pm_dec, scale = dpm_dec    , size = ntracks)
pm_ra_track  = np.random.normal(pm_ra , scale = dpm_ra     , size = ntracks)



EWtrackarray = np.zeros([ntracks, npoints])
NStrackarray = np.zeros([ntracks, npoints])
# print(np.shape(EWtrackarray))
for i in range(0,ntracks):
	# print(i,pm_ra_track[i])
	NS_track,EW_track = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra_track[i], pm_dec_track[i], prlx_track[i], tl, JDL,Nstart_track[i], Estart_track[i],t[0])
	EWtrackarray[i,:]=EW_track
	NStrackarray[i,:]=NS_track

### Now do percentile calculation 


sigma_ew = np.asarray(map(lambda v: (v[1]-v[0]),zip(*np.percentile(EWtrackarray, [50, 68], axis=0))))
sigma_ns = np.asarray(map(lambda v: (v[1]-v[0]),zip(*np.percentile(NStrackarray, [50, 68], axis=0))))

# sigma_2ew = 2.*sigma_ew
# sigma_2ew = map(lambda v: ((v[1]-v[0])/2.),zip(*np.percentile(EWtrackarray, [50, 95], axis=0)))

# EW_1sigmap = np.percentile(EWtrackarray,84, axis=0)
# EW_1sigman = EW_vector-(np.abs(EW_1sigmap - EW_vector))


# EW_2sigmap = EW_vector+(2*np.abs(EW_1sigmap-EW_vector))
# # EW_2sigmap = np.percentile(EWtrackarray,95, axis=0)

# EW_2sigman = EW_vector-(np.abs(EW_2sigmap - EW_vector))


# print(D0_NS)
# pl.figure()
# pl.plot(Nstart,'.k')
# pl.show()




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
axarr[0].fill_between(tl, NS_vector-sigma_ns,NS_vector+sigma_ns,alpha=0.25,color='black',linewidth=0)
axarr[0].fill_between(tl, NS_vector-2*sigma_ns,NS_vector+2*sigma_ns,alpha=0.2,color='black',linewidth=0)
axarr[0].errorbar(t,NS, yerr=dNS, fmt=".k") 
axarr[0].text(0.8, 0.8,target , fontweight='bold', fontsize = 25, horizontalalignment='center', verticalalignment='center', transform=axarr[0].transAxes)
axarr[0].set_ylabel("North offset (mas)",fontweight='bold',fontsize=16)
axarr[0].yaxis.set_major_locator(MaxNLocator(prune='both'))
axarr[0].locator_params(axis = 'y', nbins = 6)  


axarr[1].plot(tl, EW_vector, color="black", lw=1.5, alpha=1)
axarr[1].fill_between(tl, EW_vector-sigma_ew,EW_vector+sigma_ew,alpha=0.25,color='black',linewidth=0)
axarr[1].fill_between(tl, EW_vector-2*sigma_ew,EW_vector+2*sigma_ew,alpha=0.2,color='black',linewidth=0)
axarr[1].errorbar(t,EW, yerr=dEW, fmt=".k") 
axarr[1].set_ylabel("East offset (mas)",fontweight='bold',fontsize=16)
axarr[1].set_xlabel("Epoch (yrs)",fontweight='bold',fontsize=16)
axarr[1].yaxis.set_major_locator(MaxNLocator(prune='both'))
# axarr[1].get_xaxis().get_major_formatter().set_useOffset(False)
axarr[1].locator_params(axis = 'y', nbins = 6) #(or axis = 'y') 
pl.xlim([np.min(tl),np.max(tl)])
fig.tight_layout()
fig.subplots_adjust(hspace=0.001) # no horizontal space between figures

pl.savefig(target + '_plot_func.png')



print('Runtime: ',datetime.now() - startTime)
