'''
Calculate astrometric motion of star over given time frame using JPL Horizons ephermerides 
and SIMBAD proper motion and parallax. Plot the astrometric star path with 1 and 2 sigma uncertanties
with measured data overplotted. 

Monte Carlo error propagation:
Sample proper motion, parallax, and initial position (in North and East) from a gaussian 
to generate the simulated paths. Use percentile to take 68 and 95 positions on each timestep. 

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
from astropy.constants import G
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime
import matplotlib.ticker as ticker
import warnings
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# import seaborn as sns
warnings.filterwarnings('ignore', category=FutureWarning) # ignore warnings
# warnings.filterwarnings('ignore', category=FutureWarning) # ignore warnings

startTime = datetime.now() # time the code

def convert_coords(array):
	'''
	convert ra from hours,min,sec to radians
	and dec from deg min sec to radians 
	'''

	RADEC = SkyCoord(array[0],array[1],frame='icrs')
	RA = RADEC.ra.rad
	DEC = RADEC.dec.rad
	return RA, DEC


def proj_RA_DEC(RA, DEC, pm_ra, pm_dec, prlx, epoch, JD, D0_NS, D0_EW,date1):
	'''
    Combine proper motion and parallax to produce a stellar "track" across the Sky
    Most input variables are defined in their subfunctions
    INPUT:		
    date1 		-- the date of the first piece of data in decimalyear

    OUTPUT:
    NS_vector 	-- final offsets in delta North direction along the epoch
    EW_vector 	-- final offsets in delta East direction along the epoch

	'''

	# first calculate projections due to proper motion only
	prop_RA,prop_DEC = pm_offset(RA,DEC,pm_ra,pm_dec,epoch)
	# then calculate projections due to parallax only 
	par_RA,par_DEC = parallax_offset(prlx, prop_RA,prop_DEC,JD)	

	# Add together
	delta_RA = prop_RA + par_RA
	delta_DEC = prop_DEC + par_DEC

	# get starting positions in RA and DEC
	date1_jd = Time(date1,format = 'decimalyear').jd
	RA_start, DEC_start = get_initial_RA_DEC(RA, DEC, pm_ra,pm_dec,prlx,date1,date1_jd)

	# Subtract starting position

	delta_RA =  rad2mas(delta_RA- RA_start)
	delta_DEC = rad2mas(delta_DEC- DEC_start)


	# Convert RA and DEC to NS,EW
	delta_EW= -delta_RA*np.cos(prop_DEC)  # EW motion of a background object relative to star (mas)
	delta_NS= -delta_DEC                  # NS motion of a background object relative to star (mas)

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

	par_RA  = prlx * 1./np.cos(prop_DEC) *(x*np.sin(prop_RA) - y*np.cos(prop_RA))
	par_DEC = prlx*(x*np.cos(prop_RA)*np.sin(prop_DEC) + y*np.sin(prop_RA)*np.sin(prop_DEC) - z* np.cos(prop_RA))
	par_RA = mas2rad(par_RA)
	par_DEC = mas2rad(par_DEC)
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
	'''
	Calculate x,y,z position of geocenter relative to the barycenter at a Julian Date 
	in celestial sphere coordinates:
	x - points towards the mean equinox of epoch (E.g. J2000)
	x,y plane - along the celestial equator, projected from Earth's equator
	z - points towards the North Celestial Pole

	Note: JPL Horizons web interface returns ecliptic coordinates by default. 
		  This is wrong! 
	'''
	eph = jplephem.Ephemeris(de421)
	barycenter = eph.position('earthmoon', JD)/  1.49597870700e+8
	moonvector = eph.position('moon', JD) /  1.49597870700e+8
	earthPos = (barycenter - moonvector * eph.earth_share) 
	x = earthPos[0,:]
	y = earthPos[1,:]
	z = earthPos[2,:]
	return x,y,z

def mas2rad(mas):
	'''
	convert milliarcseconds to radians
	'''
	return mas*2*np.pi/1000/3600/360  # rads

def rad2mas(rad):
	'''
	convert radians to milliarcseconds
	'''
	return rad*360/2/np.pi*3600*1000

def Sec(theta):
	'''
	short hand secant operation on theta
	'''
	return 1/np.cos(theta)


## INPUTS: 

target    = 'HD224983AB'
params    = 'CompendiumStars.txt'
# params    = 'CompendiumStars.txt'

data_name = target + '/' + target+'.txt'
npoints   =1000            # number of points that makes up the "tornado path"
ntracks   =1000			# How many tracks used to calculate sigma from monte carlo
# triple = True 

## Read in Data
labels = np.loadtxt(params, delimiter=',', dtype=np.str, usecols=[0])
ICRS   = np.loadtxt(params, delimiter=',', dtype=np.str,usecols=[1,2,3,4])
values = np.loadtxt(params, delimiter=',', dtype=np.float,usecols=[5,6,7,8,9,10])
data   = np.loadtxt(data_name, delimiter=',', dtype=np.float)
# Find target row in params text file
# if triple == True:
a = np.where(np.char.find(labels, target[:-2]) > -1) # string comparison
# else:
# a = np.where(np.char.find(labels, target) > -1) # string comparison

ind = a[0][0]                  # save index and strip extra array things. not the cleanest..
# double check targets match
assert target[:-2] == labels[ind] ,'target star isnt listed in compendium'


# Convert input params into radians 
RA_J2000, DEC_J2000 = convert_coords(ICRS[ind,0:2]) # convert coordinates to radians
dRA, dDEC = convert_coords(ICRS[ind,2:4]) # convert coordinates to radians

values = mas2rad(values)   # convert everything else to radians now
 
platescale = 9.942

# Assign SIMBAD parameters individual names	
prlx    = rad2mas(values[ind][0])   # mas
dprlx   = rad2mas(values[ind][1])   # mas
pm_ra   = values[ind][2]   # rad/yr
dpm_ra  = values[ind][3]   # rad/yr
pm_dec  = values[ind][4]   # rad/yr
dpm_dec = values[ind][5]   # rad/yr

# Assign measured astrometry parameters individual names 
D0_NS   = data[0][3]*platescale  # NS initial positions  - mas
D0_EW   = data[0][1]*platescale  # EW initial positions  - mas
JD      = data[:,0]   # full JD
NS      = data[:,3]*platescale   # NS data
EW      = data[:,1]*platescale   # EW data
dEW     = data[:,2]*platescale   # EW error
dNS     = data[:,4]*platescale   # NS error


#### TESTING AREA %%%%%
# Assign SIMBAD parameters individual names	
# prlx    = values[ind][0]   # mas
# prlx = 2.99
# dprlx   = 0.42   # mas
# pm_ra   = mas2rad(-18.041 )   # mas/yr
# dpm_ra  = mas2rad(2.381)   # mas/yr
# pm_dec  = mas2rad(8.074)   # mas/yr
# dpm_dec = mas2rad(0.944)   # mas/yr
# # print(prlx)
# # 1/0
# # Assign measured astrometry parameters individual names 
# D0_NS   = 2.8 # NS initial positions  - mas
# D0_EW   = 3858  # EW initial positions  - mas
# JD      = [2456135.50,2456443.50]   # full JD
# NS      = [2.8,-1.3] # NS data
# EW      = [3858,3859.3]   # EW data
# dEW     = [1.8, 1.8]   # EW error
# dNS     = [1.7,1.7]   # NS error


# Convert JD into decimal year
t = Time(JD,format = 'jd').decimalyear
tl = np.linspace(np.floor(np.min(t)), np.ceil(np.max(t)), npoints)
# tl = np.linspace(2012, 2015, npoints)

# tl = np.linspace(np.min(t)-0.15, np.max(t)+0.15, npoints)
# tl = np.linspace(np.min(t), np.max(t), npoints)

# print()
JDL = Time(tl,format = 'decimalyear').jd
# JDL2 = Time(tl2,format = 'decimalyear').jd

NS_vector,EW_vector = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra, pm_dec, prlx, tl, JDL, D0_NS, D0_EW,t[0])

# SEP = np.sqrt(EW_vector**2+NS_vector**2)

# SEPA=np.sqrt(EWtrackarray**2 + NStrackarray**2) 




# # pl.plot(SEPA.T,color="black", lw=1, alpha=0.1)
# pl.plot(tl,SEP,color="blue", lw=1, alpha=1)
# pl.savefig(target +'/'+ target + 'test_plot_NSEW.pdf')


# NS2_vector,EW2_vector = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra, pm_dec, prlx, tl2, JDL2, D0_NS, D0_EW,t[0])


# fig, axarr=pl.subplots(1)
# axarr.annotate("", xy=(EW2_vector[-1],NS2_vector[-1]), xycoords='data',
#             	xytext=(EW2_vector[-10],NS2_vector[-10]), textcoords='data',
#             	arrowprops=dict(arrowstyle="-|>",
#                             connectionstyle="arc3"),)
# # axarr.annotate("", xy=(EW_vector[npoints/2],NS_vector[npoints/2]), xycoords='data',
# #             	xytext=(EW_vector[npoints/2-10],NS_vector[npoints/2-10]), textcoords='data',
# #             	arrowprops=dict(arrowstyle="-|>",
# #                             connectionstyle="arc3"),)            
# # axarr.arrow(EW_vector[-10],NS_vector[-10], EW_vector[-1],NS_vector[-1])
# axarr.plot(EW2_vector, NS2_vector, color="black", lw=1.5, alpha=1)
# axarr.errorbar(EW,NS, yerr=dEW,xerr=dNS, fmt=".k")
# # pl.xlim([-600,100])
# # pl.ylim([-900,-200])
# # pl.axis('equal')
# pl.gca().invert_xaxis()
# pl.xlabel('East Offset (mas)')
# pl.ylabel('North Offset (mas)')
# fig.savefig(target +'/'+target+ '_plot_all.pdf')


# Computing the errors - draw samples from a normal distribution for each uncertain value
RA_track = np.random.normal(RA_J2000 , scale = dRA, size = ntracks)
DEC_track = np.random.normal(DEC_J2000 , scale = dDEC, size = ntracks)
Nstart_track = np.random.normal(D0_NS , scale = dNS[0], size = ntracks)
Estart_track = np.random.normal(D0_EW , scale = dEW[0], size = ntracks)
prlx_track   = np.random.normal(prlx  , scale = dprlx      , size = ntracks)
pm_dec_track = np.random.normal(pm_dec, scale = dpm_dec    , size = ntracks)
pm_ra_track  = np.random.normal(pm_ra , scale = dpm_ra     , size = ntracks)

# create empty array of for tracks
EWtrackarray = np.zeros([ntracks, npoints])
NStrackarray = np.zeros([ntracks, npoints])

# fill array
for i in range(0,ntracks):
	# NS_track,EW_track = proj_RA_DEC(RA_J2000, DEC_J2000, pm_ra_track[i], pm_dec_track[i], prlx_track[i],\
	# tl, JDL,Nstart_track[i], Estart_track[i],t[0])
	NS_track,EW_track = proj_RA_DEC(RA_track[i], DEC_track[i], pm_ra_track[i], pm_dec_track[i], prlx_track[i],\
	tl, JDL,Nstart_track[i], Estart_track[i],t[0])


	EWtrackarray[i,:]=EW_track
	NStrackarray[i,:]=NS_track

# Now the percentile calculation 
sigma_ew = 2*np.asarray(map(lambda v: (v[1]-v[0]),zip(*np.percentile(EWtrackarray, [50, 68], axis=0))))
sigma_ns = 2*np.asarray(map(lambda v: (v[1]-v[0]),zip(*np.percentile(NStrackarray, [50, 68], axis=0))))


# SEP = np.sqrt(EW_vector**2+NS_vector**2)
# sig = np.sqrt(sigma_ew**2+sigma_ns**2)
# SEParray = np.sqrt(EWtrackarray**2+NStrackarray**2).T

# pl.figure
# pl.plot(tl,SEP,'-k')
# pl.plot(tl,SEP+sig)
# pl.plot(tl,SEP-sig)

# pl.show()
# 1/0

'''
Plotting
'''
csfont = {'fontname':'Gill Sans MT'}

# set 
majorLocator   = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocatorx   = MultipleLocator(0.1)

# minorLocatory  =  MultipleLocator(1)
minorLocator1   = AutoMinorLocator(n=2)

fig, axarr=pl.subplots(2,sharex=True)
# fig.suptitle(target, fontweight='normal', fontsize = 16,**csfont )
axarr[1].xaxis.set_major_locator(majorLocator)
axarr[1].xaxis.set_major_formatter(majorFormatter)
axarr[1].xaxis.set_minor_locator(minorLocatorx)
# axarr[1].yaxis.set_minor_locator(minorLocatory)

# Plot North offset with 1 and 2 sigma errors
axarr[0].plot(tl, NS_vector, color= 'teal', lw=1, alpha=1)
axarr[0].fill_between(tl, NS_vector-sigma_ns,NS_vector+sigma_ns,alpha=0.3,color='teal',linewidth=0)
axarr[0].fill_between(tl, NS_vector-2*sigma_ns,NS_vector+2*sigma_ns,alpha=0.2,color='teal',linewidth=0)
axarr[0].errorbar(t,NS, yerr=dNS, fmt="ok",markersize = 2.5,markerfacecolor='black',markeredgewidth=1,elinewidth=1,capsize=2) 
axarr[0].set_ylabel("$\Delta$ North (mas)",fontweight='normal',fontsize=18,labelpad = 2,**csfont)
axarr[0].yaxis.set_major_locator(MaxNLocator(prune='both'))
axarr[0].locator_params(axis = 'y', nbins = 8)  
axarr[0].yaxis.set_minor_locator(minorLocator1)
axarr[0].yaxis.set_tick_params(which='minor',width=1, length=3, color='k')
axarr[0].yaxis.set_tick_params(which='major',width=1, length=6, color='k')
axarr[0].xaxis.set_tick_params(which='minor',width=1, length=3, color='k')
axarr[0].xaxis.set_tick_params(which='major',width=1, length=6, color='k')
axarr[0].set_title(target, fontweight='normal', fontsize = 18,**csfont )

# Plot East offset with 1 and 2 sigma errors
axarr[1].plot(tl, EW_vector, color="teal", lw=1, alpha=1)
axarr[1].fill_between(tl, EW_vector-sigma_ew,EW_vector+sigma_ew,alpha=0.3,color='teal',linewidth=0)
axarr[1].fill_between(tl, EW_vector-2*sigma_ew,EW_vector+2*sigma_ew,alpha=0.2,color='teal',linewidth=0)
axarr[1].errorbar(t,EW, yerr=dEW, fmt="ok",markersize = 2.5,markerfacecolor='black',markeredgewidth=1,elinewidth=1,capsize=2)
axarr[1].set_ylabel("$\Delta$ East (mas)",fontweight='normal',fontsize=18,labelpad = 4,**csfont)
axarr[1].set_xlabel("Date (years)",fontweight='normal',fontsize=18,labelpad = -2,**csfont)
axarr[1].yaxis.set_major_locator(MaxNLocator(prune='both'))
# axarr[1].get_xaxis().get_major_formatter().set_useOffset(False)
axarr[1].locator_params(axis = 'y', nbins = 8) #(or axis = 'y') 
axarr[1].yaxis.set_minor_locator(AutoMinorLocator(n=2))
# axarr[1].xaxis.set_tick_params(width=1,length=5)
# axarr[0].xaxis.set_tick_params(width=1,length=5)
# axarr[1].yaxis.set_tick_params(width=1,length=5)
# axarr[0].yaxis.set_tick_params(width=1,length=5)

# Add label for the title
# axarr[1].text(0.85, 0.85, target , fontweight='normal', fontsize = 22, horizontalalignment='center',\
# verticalalignment='center', transform=axarr[0].transAxes,**csfont)

axarr[1].yaxis.set_tick_params(which='minor',width=1, length=3, color='k')
axarr[1].yaxis.set_tick_params(which='major',width=1, length=6, color='k')
axarr[1].xaxis.set_tick_params(which='minor',width=1, length=3, color='k')
axarr[1].xaxis.set_tick_params(which='major',width=1, length=6, color='k')
for axis in ['top','bottom','left','right']:
	axarr[1].spines[axis].set_linewidth(1.5)
	axarr[0].spines[axis].set_linewidth(1.5)

# Xlimits and  formatting 
pl.xlim([np.min(tl),np.max(tl)])
fig.tight_layout()
fig.subplots_adjust(hspace=0.001) # no horizontal space between figures
pl.savefig(target +'/'+ target + '_dNdE2.pdf')
# pl.savefig('Figures' +'/'+ target + '_dNdE.pdf')



# pl.figure()
# sns.plot(tl,np.sqrt(EW_vector**2+NS_vector**2))
# # pl.savefig('separation.png')
# pl.show()



print('Runtime: ',datetime.now() - startTime)
