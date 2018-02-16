
dprlx**2 *(np.cos(RA_J2000 + pm_ra* tl* Sec(DEC_J2000 + pm_dec* tl))* (-z + x* np.sin(DEC_J2000 + pm_dec* tl)) + y *np.sin(DEC_J2000 + pm_dec* tl)\
*np.sin(RA_J2000 + pm_ra *tl *Sec(DEC_J2000 + pm_dec* tl)))**2 + dRA**2 *prlx**2 *(y *np.cos(RA_J2000 + pm_ra* tl *Sec(DEC_J2000 + pm_dec* tl))\
*np.sin(DEC_J2000 + pm_dec* tl) + (z* - x *np.sin(DEC_J2000 + pm_dec* tl)) np.sin(RA_J2000 + pm_ra *tl *Sec(DEC_J2000 + pm_dec* tl)))**2\
+dpm_ra**2 *prlx**2* tl**2* (z *Sec(DEC_J2000 + pm_dec* tl) *np.sin(RA_J2000 +pm_ra*tl *Sec(DEC_J2000 + pm_dec*tl)) + \
(y* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl))) np.tan(DEC_J2000 + pm_dec*tl))**2\
+dDEC**2 prlx**2 (np.cos(DEC_J2000 +pm_dec*tl)* (x *np.cos(RA_J2000 + pm_ra*tl *Sec(DEC_J2000 + pm_dec*tl)) +y* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)))\
+pm_ra*tl* np.tan(DEC_J2000 +pm_dec*tl)* (z* Sec(DEC_J2000 + pm_dec*tl) np.sin(RA_J2000 + pm_ra*tl Sec(DEC_J2000 + pm_dec*tl)) + \
(y* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl))) *np.tan(DEC_J2000 + pm_dec*tl)))**2 +\
dpm_dec**2* prlx**2* tl**2* (np.cos(DEC_J2000 +pm_dec*tl)* (x* np.cos(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) +y *np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)))\
+pm_ra*tl *np.tan(DEC_J2000 + pm_dec*tl)* (z* Sec(DEC_J2000 + pm_dec*tl)* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl)) + (y* np.cos(RA_J2000 + pm_ra*tl\
*Sec(DEC_J2000 + pm_dec*tl)) - x* np.sin(RA_J2000 + pm_ra*tl* Sec(DEC_J2000 + pm_dec*tl))) *np.tan(DEC_J2000 + pm_dec*tl)))**2