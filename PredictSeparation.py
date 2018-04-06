from __future__ import print_function
import os
import matplotlib
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 
matplotlib.use('agg')
import emcee
import csv
import corner
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import astropy.units as u
from astropy.constants import G
from astropy.time import Time
from scipy import stats
from datetime import datetime
import matplotlib.ticker as ticker
# import seaborn
import warnings
from matplotlib.ticker import LinearLocator, MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
warnings.filterwarnings('ignore', category=UserWarning)
startTime = datetime.now()

# Reproducible results!
# np.random.seed(1244)

# Definitions
def radvel(t, k, w, e, p, tp, off):
    f = calc_true_anom(t, e, p, tp)
    return k*(e*np.cos(w)+np.cos(w+f))+off

def radvel_doff(t, k, w, e, p, tp, off, doff):
    f = calc_true_anom(t, e, p, tp)
    return k*(e*np.cos(w)+np.cos(w+f))+off+doff

def calc_true_anom(t, e, p, tp):
    M = (2*np.pi*(t-tp)/p)
    print(M)
    print(e)
    E = calc_e_fast(e,M)
    return np.mod(2*np.pi+2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2)),2*np.pi)

def calc_e_fast(e, M):
    Tol = 1e-6;
    E0 = M;    
    E1 = E0 + (M + e*np.sin(E0)-E0)/(1-e*np.cos(E0));
    print(np.absolute(E1-E0))
    while np.absolute(E1-E0) > Tol:
        E0 = E1;
        E1 = E0 + (M + e*np.sin(E0)-E0)/(1-e*np.cos(E0));
    return E1

def ast(t,a,p,e,inc,w,Om,tp,d):
    f = calc_true_anom(t,e,p,tp)
    # a = ((((p*secpyr)**2)* G.value * (mt*msoltokg) * 1/(4*np.pi**2) )**(1./3))*mtoAU # in meters
    dDec,dRa = calc_sep(a,e,inc,w,Om,f,d) 
    return dDec,dRa 

def calc_sep(a,e,inc,w,Om,f,d):
    r = a*(1-e**2)/(1+e*np.cos(f)) # needs to be in AU
    x = r*np.cos(f)
    y = r*np.sin(f)
    z = np.zeros((1,np.size(r)))
    Om = Om-np.pi/2
    """P1,P2,P3 are 3 orbital rotation matrices: see Solar System Dynamics Chapter 2 """

    P1 = np.array( [[np.cos(w),-np.sin(w), 0.],
                    [np.sin(w), np.cos(w), 0.],
                    [0.       , 0.       , 1.]])

    P2 = np.array( [[1., 0.         , 0.         ],
                    [0., np.cos(inc),-np.sin(inc)],
                    [0., np.sin(inc), np.cos(inc)]])

    P3 = np.array( [[np.cos(Om),-np.sin(Om), 0.],
                    [np.sin(Om), np.cos(Om), 0.],
                    [0.        , 0.        , 1.]])

    p_rot = np.dot(np.dot(P3,P2),P1)    # matrix multiplication
    coord = np.zeros((3,np.size(y)))        
    coord[0][:]=x
    coord[1][:]=y
    XYZ = np.dot(p_rot,coord)           # apply rotation to coordinates
    X = XYZ[0][:]/d                     # isolate x
    Y = XYZ[1][:]/d                     # isolate y 
    dDec = X
    dRa = Y
    return dDec,dRa                         # in arc sec...

def chi_sq(data,model,sigma,ndim, jitt):
    sigma = np.sqrt(sigma**2+jitt**2)
    inv_sigma2 = 1.0/(sigma**2)
    N = np.size(data)
    return np.sum((data - model)**2 * inv_sigma2)/(N-ndim-1)

def axis_conv(axis,prec,scale):
    string = '{0:.'+str(prec)+'0f}'
    print(string)
    return ticker.FuncFormatter(lambda axis, pos: string.format(axis*scale))

def axis_conv2(axis,scale):
    return ticker.FuncFormatter(lambda axis, pos: '{0:.2f}'.format( axis-scale))

# Analysis options
target = 'HR7672'
ndim_rv1 = 7
ndim_rv2 = 8
ndim_tot = 11

# Conversions
secpyr= 3.15569e7
sperday = 86400.
msoltokg = 1.9891e30
deg2rad = np.pi/180
rad2deg = 180./np.pi
mtoAU = 1.0/149597870700.


# Julian date offset and distance
JD_off = 2440000        # days
d = 17.76               # parsecs
a = 18.3                # AU
inc = 97.3*deg2rad
w = 259*deg2rad
Om = 61*deg2rad
tp = 2014.6
e = 0.5
t = 2457915.5-JD_off
t = Time(t+JD_off,format='jd').decimalyear
p = 73.3

print(d,a,inc,w,Om,tp,e,t,p)

dec,ra = ast(t,a,p,e,inc,w,Om,tp,d)

print(a,b)

1/0
'''
#-----------#
50pc Best fit 
#-----------#
'''

print('Computing the quantiles.')
k_pc, w_pc, e_pc, p_pc, tp_pc, off_pc,doff_pc,jitt_pc, a_pc, Om_pc, inc_pc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), 
                                                                        zip(*np.percentile(samples, [16, 50, 84], axis=0)))

# convert list outputs to array
k_pc = np.asarray(k_pc)
w_pc = np.asarray(w_pc)
e_pc = np.asarray(e_pc)
p_pc = np.asarray(p_pc)
tp_pc = np.asarray(tp_pc)
off_pc = np.asarray(off_pc)
doff_pc = np.asarray(doff_pc)
jitt_pc = np.asarray(jitt_pc)
a_pc = np.asarray(a_pc)
Om_pc = np.asarray(Om_pc)
inc_pc = np.asarray(inc_pc)

#mt  = 1/((p_pc[0]*secpyr)**2 * G.value) *1/ (msoltokg)* a_pc[0]* (4*np.pi**2) # in meters

mt_pc = (a_pc[0]**3) / (p_pc[0])**2


1/0
# Calculate Chi squared for 50th percentile best fits

rv_chisq_pc1= radvel(t1,k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0])
rv_chisq_pc2= radvel_doff(t2,k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0],doff_pc[0])

chi2_pc_rv1=chi_sq(rv1,rv_chisq_pc1,rverr1,ndim_rv1,jitt_pc[0])
chi2_pc_rv2=chi_sq(rv2,rv_chisq_pc2,rverr2,ndim_rv2,jitt_pc[0])
chi2_pc_rv = (chi2_pc_rv1+chi2_pc_rv2)/2. # ??? 

print("""50th Percentile best fit result:
    Amp.     = {0[0]} +{0[1]} -{0[2]}
    omega    = {1[0]} +{1[1]} -{1[2]}
    Ecc.     = {2[0]} +{2[1]} -{2[2]}
    Period   = {3[0]} +{3[1]} -{3[2]}
    t_p      = {4[0]} +{4[1]} -{4[2]}
    G.Offset = {5[0]} +{5[1]} -{5[2]}
    D.Offset = {6[0]} +{6[1]} -{6[2]}
    Jitter   = {7[0]} +{7[1]} -{7[2]}
    SMA   = {8[0]} +{8[1]} -{8[2]} 
    Omega    = {9[0]} +{9[1]} -{9[2]} 
    Inclin   = {10[0]} +{10[1]} -{10[2]} 
""".format(k_pc, w_pc*rad2deg, e_pc, p_pc, tp_pc, off_pc, doff_pc,jitt_pc,  
    a_pc, Om_pc*rad2deg, inc_pc*rad2deg))

# Calculate radial velocity plotting points for 50th percentile
rv_pc = radvel(tl,k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0])-off_pc[0]
# rv_pc2 = radvel_doff(tl2,k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0],doff_pc[0])
# rv_pc = np.concatenate((rv_pc1,rv_pc2-doff_pc[0]))



tax = np.linspace(np.min(ta), np.min(ta)+1.1*p_pc[0], 1000) # astrometric time baseline for full orbit plotting
dec_pc,ra_pc = ast(tax,a_pc[0],p_pc[0],e_pc[0],inc_pc[0],w_pc[0],Om_pc[0],tp_pc[0],d)
dec_pc1,ra_pc1 = ast(tax,a_pc[0],p_pc[0],e_pc[0],0,w_pc[0],Om_pc[0],tp_pc[0],d)

extra_time = [2017.0,2020.0,2025.0]
extra_DEC, extra_RA = ast(extra_time,a_pc[0],p_pc[0],e_pc[0],inc_pc[0],w_pc[0],Om_pc[0],tp_pc[0],d)


#tp_pc[0]




# print(a_pc[0],p_pc[0],e_pc[0],inc_pc[0],w_pc[0],Om_pc[0],tp_pc[0],d)
# 1/0
# combine rv data for plotting and apply offset
rv_data_pc = np.concatenate((rv1-off_pc[0],rv2-doff_pc[0]-off_pc[0]))


# Calculate residuals and RMS for plotting
y_res_pc= radvel(t, k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0])-off_pc[0]
# y2_res_pc= radvel_doff(t2,k_pc[0],w_pc[0],e_pc[0],p_pc[0],tp_pc[0],off_pc[0],doff_pc[0])-doff_pc[0]
# y_res_pc = np.concatenate((y1_res_pc,y2_res_pc))

diff_pc = y_res_pc-rv_data_pc
rms_pc = np.sqrt(np.sum(diff_pc**2)/np.size(diff_pc))


'''
#---------------------#
50th Percentile Plots
#---------------------#
'''
### Plot RV with several samples overlapping

pl.figure(facecolor='white')
for k, w, e, p, tp, off, doff, jitt, a, Om, inc in samples[np.random.randint(len(samples), size=15)]:
    pl.plot(tl, radvel(tl, k, w, e, p, tp, off)-off, color="black", alpha=0.25,linewidth=0.5)
pl.plot(tl, rv_pc, color="firebrick", lw=1, alpha=0.9)
pl.errorbar(t, rv_data_pc, yerr=rverr, fmt=".k",ms=5)
pl.xlabel("Date (yrs)",fontsize=18)
pl.ylabel("Velocity (m/s)",fontsize=18)
# pl.title(target + ' RV best fit')
# pl.savefig(path + folder +run+ target + '_plot_manyrv_pc.pdf')



# Plot RV with residuals

fig2 = pl.figure(facecolor='white')
frame1=fig2.add_axes((.1,.3,.8,.6)) #xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
pl.errorbar(t, rv_data_pc, yerr=rverr, fmt=".k")
pl.plot(tl, rv_pc, color="firebrick", lw=1, alpha=0.9)
fig2.text(0.8, 0.35,'RMS = '+str("%.2f" %rms_pc)+' m/s',fontweight='bold', fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=fig2.gca().transAxes)
frame1.set_xticklabels([]) #Remove x-tic labels
pl.ylabel("RV (m/s)",fontweight='bold')
pl.title(target+ ' RV best fit',fontsize=14,fontweight='bold')
pl.grid(lw=1.2)
frame2=fig2.add_axes((.1,0.08,.8,.19)) # Residual plot
pl.errorbar(t,diff_pc,yerr = rverr,fmt='.k')
pl.plot((frame1.get_xlim()),(0,0),'-k',lw=1)
pl.grid(lw=1.2)
pl.xlabel("Time (yrs)",fontweight='bold')
pl.ylabel("Residuals (m/s)",fontweight='bold')
# pl.savefig(path + folder +run+ target + '_plot_rv_pc_res.pdf')
pl.close()

# Plot Astrometry
fig, ax = pl.subplots()
pl.plot(0,0,'*',color='yellow',markersize=15)
for k, w, e, p, tp, off, doff, jitt, a, Om, inc in samples[np.random.randint(len(samples), size=50)]:
    dectemp,ratemp = ast(tax,a,p,e,inc,w,Om,tp,d)
    pl.plot(ratemp,dectemp, color="darkgray", alpha=1,linewidth=0.25)
pl.plot(ra_pc, dec_pc, color="firebrick", lw=1, alpha=0.9)
# pl.plot(ra_pc1, dec_pc1, color="steelblue", lw=1, alpha=0.9)
pl.errorbar(dra, ddec, yerr=dec_err, xerr = ra_err, fmt=".k")
pl.plot(extra_RA, extra_DEC,".k")

ax.set_xticks(np.arange(-1.2,1.2,0.2))
ax.set_yticks(np.arange(-1,0.5,0.2))

pl.xlabel("$\Delta$RA [arcseconds]",fontsize=18)
pl.ylabel("$\Delta$DEC [arcseconds]",fontsize=18)
ax.annotate(str(extra_time[0]), fontsize=12,xy=(extra_RA[0],extra_DEC[0]), xycoords='data',
                xytext=(extra_RA[0]-0.15,extra_DEC[0]-0.05), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"))
ax.annotate(str(extra_time[1]),fontsize=12, xy=(extra_RA[1],extra_DEC[1]), xycoords='data',
                xytext=(extra_RA[1]-0.15,extra_DEC[1]-0.05), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"))
ax.annotate(str(extra_time[2]), fontsize=12, xy=(extra_RA[2],extra_DEC[2]), xycoords='data',
                xytext=(extra_RA[2]-0.15,extra_DEC[2]-0.05), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"))


# pl.title(target+' Astrometry best fit')
pl.gca().invert_xaxis()
pl.axis('equal')
pl.savefig(path + folder +run+ target + '_plot_many_ast_test2.eps')
pl.close()



sma = samples[:,8]
period = samples[:,3]
col = 'steelblue'
mass = sma**3/period**2
num_bins = 30

fig4, axes = pl.subplots(1, 1, facecolor = 'white')
n,bins,patches=axes.hist(mass, num_bins, normed=1, facecolor=col)
axes.set_ylabel('$\mathrm{Density}$',fontsize=16)
axes.set_xlabel('$\mathrm{Total\/Mass}$',fontsize=16)
fig4.savefig(path + folder +run+ target + "mass2.eps")
bincenters = 0.5*(bins[1:]+bins[:-1])
ind = np.argmax(n)
modemass = bincenters[ind]

mass_est = np.percentile(mass, [16, 50, 84])
print(mass_est[1],modemass,mass_est[2]-mass_est[1],mass_est[1]-mass_est[0])



1/0






'''
#-----------#
Posterior Plot
#-----------#
'''

# Create histograms and pick the mode value, specifically the bin center. 
print('Calculating the mode from histograms')

num_bins = 30
bestfit = np.zeros((1,ndim_tot))
col = 'steelblue'
fig3, axes = pl.subplots(3, 4, figsize=(28, 16), facecolor = 'white')
for x in range(0,4):
    if x == 1:
        n,bins,patches=axes[0][x].hist(samples[:,x]*rad2deg, num_bins, normed=1, facecolor=col)
    else:
        n,bins,patches=axes[0][x].hist(samples[:,x], num_bins, normed=1, facecolor=col)
    # axes[0][x].get_yaxis().set_ticks([])
    axes[0][x].set_xlabel(labels[x],fontsize=28)
    axes[0][x].set_ylabel('$\mathrm{Density}$',fontsize=28)

    axes[0][x].xaxis.set_major_locator(MaxNLocator(nbins = 4, prune=None))
    axes[0][x].yaxis.set_major_locator(MaxNLocator(nbins = 5, prune='lower'))

    # elif x == 2 or 3:
    #     axes[0][x].xaxis.set_major_locator(MaxNLocator(nbins = 4, prune=None))
    #     axes[0][x].yaxis.set_major_locator(MaxNLocator(nbins = 5, prune=None))
    bincenters = 0.5*(bins[1:]+bins[:-1])
    ind = np.argmax(n)
    bestfit[0][x] = bincenters[ind]
for x in range(4,8):
    n,bins,patches=axes[1][x-4].hist(samples[:,x], num_bins, normed=1, facecolor=col)
    axes[1][x-4].set_xlabel(labels[x],fontsize=28)
    # axes[1][x-4].get_yaxis().set_ticks([])
    axes[1][x-4].set_ylabel('$\mathrm{Density}$',fontsize=28)
    axes[1][x-4].xaxis.set_major_locator(MaxNLocator(nbins = 5))
    axes[1][x-4].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
    bincenters = 0.5*(bins[1:]+bins[:-1])
    ind = np.argmax(n)
    bestfit[0][x] = bincenters[ind]
for x in range(8,11):
    n,bins,patches=axes[2][x-8].hist(samples[:,x], num_bins, normed=1, facecolor=col)
    axes[2][x-8].set_xlabel(labels[x],fontsize=28)
    axes[2][x-8].set_ylabel('$\mathrm{Density}$',fontsize=28)
    axes[2][x-8].xaxis.set_major_locator(MaxNLocator(nbins = 5))
    axes[2][x-8].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
    bincenters = 0.5*(bins[1:]+bins[:-1])
    ind = np.argmax(n)
    bestfit[0][x] = bincenters[ind]
# axes[0][1].xaxis.set_major_formatter(axis_conv(axes[0][1].xaxis.get_major_ticks(),rad2deg))
axes[1][0].xaxis.set_major_formatter(axis_conv2(axes[1][0].xaxis.get_major_ticks(),1996))
# axes[1][0].get_xaxis().get_major_formatter().set_useOffset(False)
# axes[0][2].xaxis.set_major_formatter(axis_conv(axes[0][2].xaxis.get_major_ticks(),2,))



axes[2][1].xaxis.set_major_formatter(axis_conv(axes[2][1].xaxis.get_major_ticks(),0,rad2deg))
axes[2][2].xaxis.set_major_formatter(axis_conv(axes[2][2].xaxis.get_major_ticks(),0,rad2deg))
axes[2][3].axis('off')
pl.tight_layout()
# fig3.savefig(path + folder +run+ target + "_posterior3.eps")
pl.close()
1/0
# Create "mode best fit using cumsum percentile uncertanties"
'''
#-----------#
Mode best fit
#-----------#
'''

# k_mode = np.array([bestfit[0][0],k_pc[1],k_pc[2]])
# w_mode = np.array([bestfit[0][1],w_pc[1],w_pc[2]])
# e_mode = np.array([bestfit[0][2],e_pc[1],e_pc[2]])
# p_mode = np.array([bestfit[0][3],p_pc[1],p_pc[2]])
# tp_mode = np.array([bestfit[0][4],tp_pc[1],tp_pc[2]])
# off_mode = np.array([bestfit[0][5],off_pc[1],off_pc[2]])
# doff_mode = np.array([bestfit[0][6],doff_pc[1],doff_pc[2]])
# jitt_mode = np.array([bestfit[0][7],jitt_pc[1],jitt_pc[2]])
# a_mode = np.array([bestfit[0][8],a_pc[1],a_pc[2]])
# Om_mode = np.array([bestfit[0][9],Om_pc[1],Om_pc[2]])
# inc_mode = np.array([bestfit[0][10],inc_pc[1],inc_pc[2]])

# mt_mode = (a_mode[0]**3) / (p_mode[0])**2

# # Calculate Chi square for mode best fit

# rv_chisq_mode1 = radvel(t1,k_mode[0],w_mode[0],e_mode[0],p_mode[0],tp_mode[0],off_mode[0])
# rv_chisq_mode2 = radvel_doff(t2,k_mode[0],w_mode[0],e_mode[0],p_mode[0],tp_mode[0],off_mode[0],doff_mode[0])
# chisq_mode_rv1 = chi_sq(rv1,rv_chisq_mode1,rverr1,ndim_rv1,jitt_mode[0])
# chisq_mode_rv2 = chi_sq(rv2,rv_chisq_mode2,rverr2,ndim_rv2,jitt_mode[0])
# chisq_mode_rv  = (chisq_mode_rv1+chisq_mode_rv2)/2.


# # Print out mode results
# print("""Mode best fit result:
#     Amp.   = {0[0]} +{0[1]} -{0[2]}
#     omega  = {1[0]} +{1[1]} -{1[2]}
#     Ecc.   = {2[0]} +{2[1]} -{2[2]}
#     Period = {3[0]} +{3[1]} -{3[2]}
#     t_p    = {4[0]} +{4[1]} -{4[2]}
#     G.Offset = {5[0]} +{5[1]} -{5[2]}
#     D.Offset = {6[0]} +{6[1]} -{6[2]}
#     Jitter   = {7[0]} +{7[1]} -{7[2]}
#     SMA      = {8[0]} +{8[1]} -{8[2]} 
#     Omega    = {9[0]} +{9[1]} -{9[2]} 
#     Inclin   = {10[0]} +{10[1]} -{10[2]}
#     Mass     = {11} 
#     Chi2     = {12}
# """.format(k_mode, w_mode*rad2deg, e_mode, p_mode, tp_mode, off_mode, doff_mode,jitt_mode,  
#     a_mode, Om_mode*rad2deg,  inc_mode*rad2deg,mt_mode,chisq_mode_rv))


# # Calculate mode best fit rvs and astrometry for plotting
# rv_mode  = radvel(tl,k_mode[0],w_mode[0],e_mode[0],p_mode[0],tp_mode[0],off_mode[0])-off_mode[0]
# dec_mode,ra_mode = ast(tax,a_mode[0],p_mode[0],e_mode[0],inc_mode[0],w_mode[0],Om_mode[0],tp_mode[0],d)


# # combine rv data for plotting and apply offset
# rv_data_mode = np.concatenate((rv1-off_mode[0],rv2-doff_mode[0]-off_mode[0]))


# y_res_mode= radvel(t, k_mode[0],w_mode[0],e_mode[0],p_mode[0],tp_mode[0],off_mode[0])-off_mode[0]

# diff_mode = y_res_mode-rv_data_mode
# rms_mode = np.sqrt(np.sum(diff_mode**2)/np.size(diff_mode))


# '''
# #---------#
# Mode Plots
# #---------#
# '''


# Save data to file

print("Writing results to file...")


with open(path + folder +run+ target +'_analysis_results.txt','wb') as csvfile:
    writer = csv.writer(csvfile,delimiter='\t') 
    writer.writerow(['\n', startTime])
    writer.writerow(["target:"+ target])                 # writing a string + a single value
    writer.writerow(["Run: "+ run])                 # writing a string + a single value
    writer.writerow(["distance: "  "%.3f" % d]) 
    writer.writerow(["JD offset: "  "%.2f" % JD_off])   
    writer.writerow(['                       ','50th PC','+error','-error'])   
    writer.writerow(['Semi amplitude (K)      ',"%.6f" % k_pc[0], "%.6f" % k_pc[1],"%.6f" % k_pc[2]])  
    writer.writerow(['Arg. of periastron (w)  ',"%.6f" % w_pc[0], "%.6f" % w_pc[1],"%.6f" % w_pc[2]]) 
    writer.writerow(['Eccentricity (e)        ',"%.6f" % e_pc[0], "%.6f" % e_pc[1],"%.6f" % e_pc[2]]) 
    writer.writerow(['Period (p)              ',"%.6f" % p_pc[0], "%.6f" % p_pc[1],"%.6f" % p_pc[2]]) 
    writer.writerow(['Time of Periastron (Tp) ',"%.6f" % tp_pc[0], "%.6f" % tp_pc[1],"%.6f" % tp_pc[2]])
    writer.writerow(['Global offset           ',"%.6f" % off_pc[0], "%.6f" % off_pc[1],"%.6f" % off_pc[2]])
    writer.writerow(['Diff. offset            ',"%.6f" % doff_pc[0], "%.6f" % doff_pc[1],"%.6f" % doff_pc[2]])
    writer.writerow(['Jitter                  ',"%.6f" % jitt_pc[0], "%.6f" % jitt_pc[1],"%.6f" % jitt_pc[2]])
    writer.writerow(['semi-major axis (a)     ',"%.6f" % a_pc[0], "%.6f" % a_pc[1],"%.6f" % a_pc[2]])
    writer.writerow(['Long. of asc. node (Om) ',"%.6f" % Om_pc[0], "%.6f" % Om_pc[1],"%.6f" % Om_pc[2]]) 
    writer.writerow(['Inclination (i)         ',"%.6f" % inc_pc[0], "%.6f" % inc_pc[1],"%.6f" % inc_pc[2]])


print('Done.\n')
print('Runtime: \n',datetime.now() - startTime)
