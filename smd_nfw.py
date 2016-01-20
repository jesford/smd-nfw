from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units
from scipy.integrate import quad, dblquad

import utils

#------------------------------------------------------------------------------

class SurfaceMassDensity(object):
    """Calculate NFW profiles for Sigma and Delta-Sigma."""
    def __init__(self, rs, delta_c, rho_crit, sig_offset=None, rbins=None):
        
        if rbins is None:
            rmin, rmax = 0.1, 5. 
            rbins = np.logspace(np.log10(rmin), np.log10(rmax), num = 50)
            self._rbins = rbins * units.Mpc
        else:
            #check rbins input units & type
            self._rbins = utils.check_input(rbins, units.Mpc)
        
        #check rs input units
        self._rs = utils.check_input(rs, units.Mpc)

        #check rho_crit input units & type
        self._rho_crit = utils.check_input(rho_crit, units.Msun/units.Mpc/(units.pc**2))
              
        #check delta_c input units & type
        if hasattr(delta_c, 'unit'):
            raise ValueError("delta_c should be a dimensionless quantity.")
        self._delta_c = utils.check_array_or_list(delta_c)

        self._nbins = self._rbins.shape[0]
        self._nlens = self._rs.shape[0]
            
        if self._delta_c.shape[0] != self._nlens or self._rho_crit.shape[0] != self._nlens:
            raise ValueError("Input arrays (rs, delta_c, rho_crit) must have \
                              the same length (the number of clusters).")
        else:
            rs_dc_rcrit = self._rs * self._delta_c * self._rho_crit
            self._rs_dc_rcrit = rs_dc_rcrit.reshape(self._nlens,
                                                    1).repeat(self._nbins,1)

        if sig_offset is not None:
            self._sigmaoffset = utils.check_input(sig_offset, units.Mpc)
            if self._sigmaoffset.shape[0] != self._nlens:
                raise ValueError("sig_offset array must have length equal to \
                                  the number of clusters.")
        else:
            self._sigmaoffset = sig_offset #None

        #set self._x, self._x_big, self._x_small, self._x_one
        _set_dimensionless_radius(self)

            

    def sigma_nfw(self, offsets = None):
        """Returns NFW surface mass density profile (centered)."""

        def _centered_sigma(self, singlecluster = None):
            #perfectly centered cluster case
            
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
            if singlecluster is None:
                sigma = 2. * self._rs_dc_rcrit * f
            else:
                sigma = 2. * self._rs_dc_rcrit[singlecluster,0] * f

            return sigma

        
        def _offset_sigma(self):
            
            def dbl_integrand(theta, r_off, r_Mpc, ith_cluster):
                #theta, r_off are integration variables
                #r_Mpc is one of the radial measurement bins

                #temporarily make _x correspond to r_eq13
                r_eq13 = np.sqrt(r_Mpc**2 + r_off**2 -
                                 2.*r_Mpc*r_off*np.cos(theta))
                #r_eq13 = np.array([r_eq13])
                #print('\n r_EQ13: ', r_eq13, '\n')
                _set_dimensionless_radius(self, radii = r_eq13, singlecluster = ith_cluster)
                sigma = _centered_sigma(self, singlecluster = ith_cluster)
                #sigma = _sigma_singlecluster(self,
                #print('sigma\n', sigma, '\n')
                                
                inner_integral = sigma.value/(2.*np.pi)
                #print('inner_integral.shape', inner_integral.shape)
                #print('offsets', offsets)
                PofRoff = (r_off/(offsets.value[ith_cluster]**2) *
                           np.exp(-0.5 * (r_off/offsets.value[ith_cluster])**2))
                full_integrand = PofRoff * inner_integral
                return full_integrand

            sigma_sm = [] #TO DO: change to predefined np arrays?
            error_sm = []

            #for each cluster and physical measurement radius, calculate 
            # double integral to get sigma_smoothed & its error...
            # change eps* ? TO DO: add option to customize precision.
            rbins_saved = self._rbins.value

            rsdcrc_percluster = self._rs_dc_rcrit[:,0]

            #for every cluster and radius, do double integration
            for i in range(self._nlens):
                sigma_sm_ithcluster = []
                error_sm_ithcluster = []
                print('\nCalculating for cluster', i+1, 'of', self._nlens)
                
                for radius in rbins_saved:
                    #print('\nradius', radius)

                    #inner integral integrates theta 0 -> 2pi
                    #outer integral integrates r_off 0 -> Inf
                    I = dblquad(dbl_integrand, 0, np.inf,
                                lambda x: 0, lambda x: 2.*np.pi,
                                args = (radius,i), epsabs=0.1, epsrel=0.1)

                    sigma_sm_ithcluster.append(I[0])
                    error_sm_ithcluster.append(I[1])
                sigma_sm.append(sigma_sm_ithcluster)
                error_sm.append(error_sm_ithcluster)

            #reset _x to correspond to input rbins (default)
            _set_dimensionless_radius(self)

            sigma_sm = np.array(sigma_sm) * units.solMass / units.pc**2
            error_sm = np.array(error_sm) * units.solMass / units.pc**2
            
            return sigma_sm, error_sm


        if offsets is None:
            finalsigma = _centered_sigma(self)
            errsigma = None
        elif np.abs(offsets).sum() == 0:
            finalsigma = _centered_sigma(self)
            errsigma = None
        else:
            if type(offsets) == list:
                offsets = np.array(list)
            finalsigma, errsigma = _offset_sigma(self)
            
            
        return finalsigma #, errsigma ???

#------------------------------------------------------------------------------

    
    def deltasigma_nfw(self):
        """Returns NFW differential surface mass density profile (centered)."""
        #calculate g

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

        #calculate & return centered profile
        deltasigma = self._rs_dc_rcrit * g
        
        return deltasigma


    
def _set_dimensionless_radius(self, radii = None, singlecluster = None):
    if radii is None:
        radii = self._rbins
    #calculate x = radii / rs
    if singlecluster is None:
        radii_repeated = radii.reshape(self._nbins,1).repeat(self._nlens,1)
        rs_repeated = self._rs.reshape(self._nlens,1).repeat(self._nbins,1)
        x = radii_repeated.T/rs_repeated
    else:
        x = radii.reshape(1,1)/self._rs[singlecluster]

    #dimensionless radius
    self._x = x.value

    #set the 3 cases of dimensionless radius x
    self._x_small = np.where(self._x < 1.-1.e-6)
    self._x_big = np.where(self._x > 1.+1.e-6)
    self._x_one = np.where(np.abs(self._x-1) <= 1.e-6)



