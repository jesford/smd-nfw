from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units


def smd_nfw(rs, delta_c, rho_crit, sig_offset=None, rbins=None):

    if rbins is None:  #default r bins
        rmin, rmax = 0.1, 5. 
        rbins = np.logspace(np.log10(rmin), np.log10(rmax), num = 10)
        rbins = rbins * units.Mpc
    nbins = rbins.shape[0]


    nlens = rs.shape[0]
    if delta_c.shape[0] != nlens or rho_crit.shape[0] != nlens:
        raise ValueError("Input arrays (rs, delta_c, rho_crit) must have \
                          the same length (the number of clusters).")
    #rho_crit = rho_crit * 10.**6 * units.Msun / (units.Mpc * units.pc * units.pc)

    #if sig_offset is not None:
    #    sig_offset = sig_offset * units.Mpc

    rbins_repeated = rbins.reshape(nbins,1).repeat(nlens,1)
    rs_repeated = rs.reshape(nlens,1).repeat(nbins,1)
    x = rbins_repeated.T/rs_repeated
    x = x.value

    #the 3 cases of dimensionless radius x
    x_small = np.where(x < 1.-1.e-6)
    x_big = np.where(x > 1.+1.e-6)
    x_one = np.where(np.abs(x-1) <= 1.e-6)

    #--------------------------------------------------------------------------
    # calculate f

    bigF = np.zeros_like(x)
    f = np.zeros_like(x)

    bigF[x_small] = (np.log((1./x[x_small]) + np.sqrt((1./(x[x_small]**2))-1.))
                    / np.sqrt(1.-(x[x_small]**2)))

    bigF[x_big] = np.arccos(1./x[x_big]) / np.sqrt(x[x_big]**2 - 1.)

    f = (1. - bigF) / (x**2 - 1.)
    f[x_one] = 1./3.
    if np.isnan(np.sum(f)) or np.isinf(np.sum(f)):
        print('\nERROR: f is not all real\n', f)


    #--------------------------------------------------------------------------
    # calculate g

    firstpart = np.zeros_like(x)
    secondpart = np.zeros_like(x)
    g = np.zeros_like(x)

    firstpart[x_small] = (( (4./x[x_small]**2) + (2./(x[x_small]**2 - 1.)) )
                          / np.sqrt(1. - x[x_small]**2))
    firstpart[x_big] = (8./(x[x_big]**2 * np.sqrt(x[x_big]**2 - 1.)) +
                        4./((x[x_big]**2-1.)**1.5))

    secondpart[x_small] = (np.log((1. + np.sqrt((1. - x[x_small])/
                           (1. + x[x_small])))/(1. - np.sqrt((1. - x[x_small])
                            /(1. + x[x_small])))))
    secondpart[x_big] = np.arctan(np.sqrt((x[x_big] - 1.)/(1. + x[x_big])))

    g = firstpart*secondpart + (4./(x**2))*np.log(x/2.) - (2./(x**2-1.))
    g[x_one] = (10./3.) + 4.*np.log(0.5)
    if np.isnan(np.sum(g)) or np.isinf(np.sum(g)):
        print('\nERROR: g is not all real\n', g)

    #--------------------------------------------------------------------------
    # calculate h

    h = np.zeros_like(x)
    h = (bigF + np.log(x/2.))/(x**2)
    h[x_one] = 1. + np.log(0.5)
    if np.isnan(np.sum(h)) or np.isinf(np.sum(h)):
        print('\nERROR: h is not all real\n', h)

    #--------------------------------------------------------------------------
    # calculate & return centered profiles

    rs_dc_rcrit = rs*delta_c*rho_crit
    rs_dc_rcrit_repeated = rs_dc_rcrit.value.reshape(nlens,1).repeat(nbins,1)

    sigma_nfw = 2. * rs_dc_rcrit_repeated * f
    mean_inside_sigma_nfw = 4. * rs_dc_rcrit_repeated * h
    deltasigma_nfw = mean_inside_sigma_nfw - sigma_nfw

    return sigma_nfw, deltasigma_nfw
    
