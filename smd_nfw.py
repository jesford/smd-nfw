from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units

rbins = np.loadtxt('smd_in1.dat')
nbins = rbins.shape[0]
rbins = rbins * units.Mpc

data = np.loadtxt('smd_in2.dat')
nlens = data.shape[0]

rs = data[:,0] * units.Mpc
delta_c = data[:,1] #dimensionless
rho_crit = data[:,2] #assumed to be in Msun/pc^3 
sig_center = data[:,3] * units.Mpc

rho_crit = rho_crit * 10.**6 * units.Msun / (units.Mpc * units.pc * units.pc)

#if sig_center[0] == 0:

rbins_repeated = rbins.reshape(nbins,1).repeat(nlens,1)
rs_repeated = rs.reshape(nlens,1).repeat(nbins,1)
x = rbins_repeated.T/rs_repeated
x = x.value

#the 3 cases of dimensionless radius x
x_small = np.where(x < 1.-1.e-6)
x_big = np.where(x > 1.+1.e-6)
x_one = np.where(np.abs(x-1) <= 1.e-6)

#------------------------------------------------------------------------------
# calculate f

bigF = np.zeros_like(x)
f = np.zeros_like(x)

bigF[x_small] = (np.log((1./x[x_small]) + np.sqrt((1./(x[x_small]**2)) - 1.))
                / np.sqrt(1.-(x[x_small]**2)))
#check notes for whether acos == arccos?
bigF[x_big] = np.arccos(1./x[x_big]) / np.sqrt(x[x_big]**2 - 1.)

f = (1. - bigF) / (x**2 - 1.)
f[x_one] = 1./3.
if np.isnan(np.sum(f)) or np.isinf(np.sum(f)):
    print('\nERROR: f is not all real\n', f)


#------------------------------------------------------------------------------
# calculate g

firstpart = np.zeros_like(x)
secondpart = np.zeros_like(x)
g = np.zeros_like(x)

firstpart[x_small] = (( (4./x[x_small]**2) + (2./(x[x_small]**2 - 1.)) )
                      / np.sqrt(1. - x[x_small]**2))
firstpart[x_big] = (8./(x[x_big]**2 * np.sqrt(x[x_big]**2 - 1.)) +
                    4./((x[x_big]**2-1.)**1.5))

secondpart[x_small] = np.log((1. + np.sqrt((1. - x[x_small])/(1. + x[x_small])))/
                             (1. - np.sqrt((1. - x[x_small])/(1. + x[x_small]))))
secondpart[x_big] = np.arctan(np.sqrt((x[x_big] - 1.)/(1. + x[x_big])))

g = firstpart*secondpart + (4./(x**2))*np.log(x/2.) - (2./(x**2-1.))
g[x_one] = (10./3.) + 4.*np.log(0.5)
if np.isnan(np.sum(g)) or np.isinf(np.sum(g)):
    print('\nERROR: g is not all real\n', g)

#------------------------------------------------------------------------------
# calculate h

h = np.zeros_like(x)
h = (bigF + np.log(x/2.))/(x**2)
h[x_one] = 1. + np.log(0.5)
if np.isnan(np.sum(h)) or np.isinf(np.sum(h)):
    print('\nERROR: h is not all real\n', h)


#------------------------------------------------------------------------------
# calculate centered profiles

rs_dc_rcrit = rs*delta_c*rho_crit
rs_dc_rcrit_repeated = rs_dc_rcrit.value.reshape(nlens,1).repeat(nbins,1)

sigma_nfw = 2. * rs_dc_rcrit_repeated * f
mean_inside_sigma_nfw = 4. * rs_dc_rcrit_repeated * h
deltasigma_nfw = mean_inside_sigma_nfw - sigma_nfw

np.savetxt('sigma_PYTHON.dat', sigma_nfw, fmt='%15.8g')
np.savetxt('deltasigma_PYTHON.dat', deltasigma_nfw, fmt='%15.8g')
