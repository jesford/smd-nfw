from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units
from scipy.integrate import quad, dblquad


#------------------------------------------------------------------------------

class SurfaceMassDensity(object):
    """Calculate NFW profiles for Sigma and Delta-Sigma."""
    def __init__(self, rs, delta_c, rho_crit, sig_offset=None, rbins=None):

        #units_message = " must be in units of "
        
        if rbins is None:
            rmin, rmax = 0.1, 5. 
            rbins = np.logspace(np.log10(rmin), np.log10(rmax), num = 50)
            self._rbins = rbins * units.Mpc
        else:
            #check rbins input units
            self._rbins = _check_units(rbins, units.Mpc)
            #if hasattr(rbins, 'unit'):
            #    if rbins.unit != units.Mpc:
            #        raise ValueError("rbins" + units_message + "Mpc.")
            #    self._rbins = rbins
            #else:
            #    self._rbins = rbins * units.Mpc
            
            
        self._nbins = self._rbins.shape[0]

        #check rs input units
        self._rs = _check_units(rs, units.Mpc)
        #if hasattr(rs, 'unit'):
        #    if rs.unit != units.Mpc:
        #        raise ValueError("rs" + units_message + "Mpc.")
        #    self._rs = rs
        #else:
        #    self._rs = rs * units.Mpc

        #check rho_crit input units
        self._rho_crit = _check_units(rho_crit, units.Msun/units.Mpc/(units.pc**2))
        #if hasattr(rho_crit, 'unit'):
        #    if rho_crit.unit != units.Msun/units.Mpc/(units.pc**2):
        #        raise ValueError("rho_crit" + units_message + "Msun/Mpc/pc**2.")
        #else:
        #    rho_crit = rho_crit * units.Msun/units.Mpc/(units.pc**2)
        #self._rho_crit = rho_crit
              
        #check delta_c input units
        if hasattr(delta_c, 'unit'):
            raise ValueError("delta_c should be a dimensionless quantity.")
        self._delta_c = delta_c
                   
        self._nlens = rs.shape[0]
        if delta_c.shape[0] != self._nlens or rho_crit.shape[0] != self._nlens:
            raise ValueError("Input arrays (rs, delta_c, rho_crit) must have \
                              the same length (the number of clusters).")
                              
        else:
            rs_dc_rcrit = self._rs * self._delta_c * self._rho_crit
            self._rs_dc_rcrit = rs_dc_rcrit.reshape(self._nlens,
                                                    1).repeat(self._nbins,1)
            #self._rs = rs

                  
        if sig_offset is not None:
            if sig_offset.shape[0] != self._nlens:
                raise ValueError("sig_offset array must have length equal to \
                                  the number of clusters.")
        self._sigmaoffset = sig_offset #None or array

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
        #print('radii', radii)
        #print('self._rs[singlecluster]', self._rs[singlecluster])
        x = radii.reshape(1,1)/self._rs[singlecluster]

    #dimensionless radius
    self._x = x.value

    #set the 3 cases of dimensionless radius x
    self._x_small = np.where(self._x < 1.-1.e-6)
    self._x_big = np.where(self._x > 1.+1.e-6)
    self._x_one = np.where(np.abs(self._x-1) <= 1.e-6)



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


def _check_units(input, expected_units):
    if hasattr(input, 'unit'):
        if input.unit != expected_units:
            raise ValueError('Expected units of ' + str(expected_units))
    else:
        input = input * expected_units
    return input
