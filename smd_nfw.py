from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units


#------------------------------------------------------------------------------

class SurfaceMassDensity(object):
    """Calculate NFW profiles for Sigma and Delta-Sigma."""
    def __init__(self, rs, delta_c, rho_crit, sig_offset=None, rbins=None):
        if rbins is None:
            rmin, rmax = 0.1, 5. 
            rbins = np.logspace(np.log10(rmin), np.log10(rmax), num = 50)
            rbins = rbins * units.Mpc
        else:
            rbins = rbins
        nbins = rbins.shape[0]
        
        nlens = rs.shape[0]
        if delta_c.shape[0] != nlens or rho_crit.shape[0] != nlens:
            raise ValueError("Input arrays (rs, delta_c, rho_crit) must have \
                              the same length (the number of clusters).")
        else:
            rs_dc_rcrit = rs * delta_c * rho_crit
            self._rs_dc_rcrit = rs_dc_rcrit.value.reshape(nlens,
                                                          1).repeat(nbins,1)
      
        if sig_offset is not None:
            if sig_offset.shape[0] != nlens:
                raise ValueError("sig_offset array must have length equal to \
                                  the number of clusters.")
        self._sigmoffset = sig_offset #None or array

        rbins_repeated = rbins.reshape(nbins,1).repeat(nlens,1)
        rs_repeated = rs.reshape(nlens,1).repeat(nbins,1)
        x = rbins_repeated.T/rs_repeated
        self._x = x.value

        #the 3 cases of dimensionless radius x
        self._x_small = np.where(self._x < 1.-1.e-6)
        self._x_big = np.where(self._x > 1.+1.e-6)
        self._x_one = np.where(np.abs(self._x-1) <= 1.e-6)

        
    def sigma_nfw(self):

        #calculate f
        bigF = np.zeros_like(self._x)
        f = np.zeros_like(self._x)

        bigF[self._x_small] = (np.log((1./self._x[self._x_small]) +
                               np.sqrt((1./(self._x[self._x_small]**2))-1.))
                               / np.sqrt(1.-(self._x[self._x_small]**2)))

        bigF[self._x_big] = (np.arccos(1./self._x[self._x_big])
                             / np.sqrt(self._x[self._x_big]**2 - 1.))

        f = (1. - bigF) / (self._x**2 - 1.)
        f[self._x_one] = 1./3.
        if np.isnan(np.sum(f)) or np.isinf(np.sum(f)):
            print('\nERROR: f is not all real\n', f)

        #calculate & return centered profiles
        sigma = 2. * self._rs_dc_rcrit * f

        return sigma
    

#------------------------------------------------------------------------------

    
    def deltasigma_nfw(self):

        # calculate g

        firstpart = np.zeros_like(self._x)
        secondpart = np.zeros_like(self._x)
        g = np.zeros_like(self._x)

        firstpart[self._x_small] = (( (4./self._x[self._x_small]**2) +
                                      (2./(self._x[self._x_small]**2 - 1.)) )
                                    / np.sqrt(1. - self._x[self._x_small]**2))
        
        firstpart[self._x_big] = (8./(self._x[self._x_big]**2 *
                                      np.sqrt(self._x[self._x_big]**2 - 1.)) +
                                        4./((self._x[self._x_big]**2-1.)**1.5))

        secondpart[self._x_small] = (np.log((1. + np.sqrt((1. -
                                                    self._x[self._x_small])/
                                    (1. + self._x[self._x_small])))/
                                    (1. - np.sqrt((1. - self._x[self._x_small])
                                    / (1. + self._x[self._x_small])))))
        
        secondpart[self._x_big] = np.arctan(np.sqrt((self._x[self._x_big] - 1.)
                                                / (1. + self._x[self._x_big])))

        g = firstpart*secondpart + ((4./(self._x**2))*np.log(self._x/2.) -
                                    (2./(self._x**2-1.)))
        g[self._x_one] = (10./3.) + 4.*np.log(0.5)
        if np.isnan(np.sum(g)) or np.isinf(np.sum(g)):
            print('\nERROR: g is not all real\n', g)

        # calculate & return centered profile
        deltasigma = self._rs_dc_rcrit * g
        
        return deltasigma






    #--------------------------------------------------------------------------
    # calculate h

    #h = np.zeros_like(x)
    #h = (bigF + np.log(x/2.))/(x**2)
    #h[x_one] = 1. + np.log(0.5)
    #if np.isnan(np.sum(h)) or np.isinf(np.sum(h)):
    #    print('\nERROR: h is not all real\n', h)

    #--------------------------------------------------------------------------


    #mean_inside_sigma_nfw = 4. * rs_dc_rcrit_repeated * h
    #deltasigma_nfw = mean_inside_sigma_nfw - sigma_nfw

    #alternate equivalent calculation (don't need sigma_nfw directly):
    #deltasigma_nfw = rs_dc_rcrit_repeated * g
    #np.testing.assert_allclose(deltasigma_nfw, rs_dc_rcrit_repeated * g, rtol=10**-3)
