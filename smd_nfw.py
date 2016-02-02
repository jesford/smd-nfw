from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units
from scipy.integrate import quad, dblquad

import utils

#------------------------------------------------------------------------------

class SurfaceMassDensity(object):
    """Calculate NFW profiles for Sigma and Delta-Sigma.

    Parameters
    ----------
    rs : array_like
        Scale radii (in Mpc) for every halo. Should be 1D, optionally with
        astropy.units of Mpc.
    delta_c : array_like
        Characteristic overdensities for every halo. Should be 1D and have
        the same length as rs.
    rho_crit : array_like
        Critical energy density of the universe (in Msun/Mpc/pc^2) at every
        halo z. Should be 1D, optionally with astropy.units of
        Msun/Mpc/(pc**2), and have the same length as rs.
    offsets : array_like, optional
        Width of the Gaussian distribution of miscentering offsets (in
        Mpc), for every cluster halo. Should be 1D, optionally with
        astropy.units of Mpc, and have the same length as rs. (Note: it is
        common to use the same value for every halo, implying that they are
        drawn from the same offset distribution).
    rbins : array_like, optional
        Radial bins (in Mpc) at which profiles will be calculated. Should
        be 1D, optionally with astropy.units of Mpc.

    Methods
    ----------
    sigma_nfw(epsabs=0.1, epsrel=0.1)
        Calculate surface mass density Sigma.
    deltasigma_nfw()
        Calculate differential surface mass density DeltaSigma.

    See Also
    ----------
    ClusterEnsemble : parameters and profiles for a sample of clusters.
        This class provides an interface to SurfaceMassDensity, and tracks
        a DataFrame of parameters as well as nfw profiles for many clusters
        at once, only requiring the user to specify cluster z and richness,
        at a minimum.
    """
    def __init__(self, rs, delta_c, rho_crit, offsets=None, rbins=None):
        if rbins is None:
            rmin, rmax = 0.1, 5. 
            rbins = np.logspace(np.log10(rmin), np.log10(rmax), num = 50)
            self._rbins = rbins * units.Mpc
        else:
            #check rbins input units & type
            self._rbins = utils.check_units_and_type(rbins, units.Mpc)
        
        #check input units & types
        self._rs = utils.check_units_and_type(rs, units.Mpc)
        self._delta_c = utils.check_units_and_type(delta_c, None)
        self._rho_crit = utils.check_units_and_type(rho_crit,
                                        units.Msun/units.Mpc/(units.pc**2))

        self._nbins = self._rbins.shape[0]
        self._nlens = self._rs.shape[0]

        if offsets is not None:
            self._sigmaoffset = utils.check_units_and_type(offsets, units.Mpc)
            utils.check_input_size(self._sigmaoffset, self._nlens)
        else:
            self._sigmaoffset = offsets #None

        #check array sizes are compatible
        utils.check_input_size(self._rs, self._nlens)
        utils.check_input_size(self._delta_c, self._nlens)
        utils.check_input_size(self._rho_crit, self._nlens)
        utils.check_input_size(self._rbins, self._nbins)
        
        rs_dc_rcrit = self._rs * self._delta_c * self._rho_crit
        self._rs_dc_rcrit = rs_dc_rcrit.reshape(self._nlens,
                                                    1).repeat(self._nbins,1)
        
        #set self._x, self._x_big, self._x_small, self._x_one
        _set_dimensionless_radius(self)

#------------------------------------------------------------------------------

    def sigma_nfw(self, epsabs=0.1, epsrel=0.1):
        """Calculate NFW surface mass density profile.

        Generate the surface mass density profiles of each cluster halo,
        assuming a spherical NFW model. Optionally includes the effect of
        cluster miscentering offsets, if the parent object was initialized
        with offsets.

        Parameters
        ----------
        epsabs, epsrel : float, optional
            Absolute and relative tolerances of the double integration in
            the miscentering calculations (no effect for offsets=None).
            Defaults are both currently set to 0.1.

        Returns
        ----------
        Quantity
            Surface mass density profiles (ndarray, in astropy.units of
            Msun/pc/pc). Each row corresponds to a single cluster halo.
        """
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
                _set_dimensionless_radius(self, radii = r_eq13,
                                          singlecluster = ith_cluster)
                
                sigma = _centered_sigma(self, singlecluster = ith_cluster)
                inner_integral = sigma.value/(2.*np.pi)

                PofRoff = (r_off/(self._sigmaoffset.value[ith_cluster]**2) *
                           np.exp(-0.5 * (r_off/
                                    self._sigmaoffset.value[ith_cluster])**2))
                    
                full_integrand = PofRoff * inner_integral
                return full_integrand

            sigma_sm = [] #TO DO: change to predefined np arrays?
            error_sm = []

            #for each cluster and physical measurement radius, calculate 
            # double integral to get sigma_smoothed & its error...
            rbins_saved = self._rbins.value

            rsdcrc_percluster = self._rs_dc_rcrit[:,0]
            
            #for every cluster and radius, do double integration
            for i in range(self._nlens):
                sigma_sm_ithcluster = []
                error_sm_ithcluster = []
                #print('\nCalculating for cluster', i+1, 'of', self._nlens)
                
                for radius in rbins_saved:
                    #print('\nradius', radius)

                    #inner integral integrates theta: 0 -> 2pi
                    #outer integral integrates r_off: 0 -> Inf
                    I = dblquad(dbl_integrand, 0, np.inf,
                                lambda x: 0, lambda x: 2.*np.pi,
                                args = (radius,i))#, epsabs=epsabs, epsrel=epsrel)

                    sigma_sm_ithcluster.append(I[0])
                    error_sm_ithcluster.append(I[1])
                sigma_sm.append(sigma_sm_ithcluster)
                error_sm.append(error_sm_ithcluster)

            #reset _x to correspond to input rbins (default)
            _set_dimensionless_radius(self)

            sigma_sm = np.array(sigma_sm) * units.solMass / units.pc**2
            error_sm = np.array(error_sm) * units.solMass / units.pc**2
            
            return sigma_sm, error_sm


        if self._sigmaoffset is None:
            finalsigma = _centered_sigma(self)
            errsigma = None
        elif np.abs(self._sigmaoffset).sum() == 0:
            finalsigma = _centered_sigma(self)
            errsigma = None
        else:
            finalsigma, errsigma = _offset_sigma(self)
            
        print('\nThis is the slow dblquad python version!\n')
        return finalsigma #, errsigma ???

#------------------------------------------------------------------------------

    
    def deltasigma_nfw(self):
        """Calculate NFW differential surface mass density profile.

        Generate the surface mass density profiles of each cluster halo,
        assuming a spherical NFW model. Currently calculates centered
        profiles ONLY; DOES NOT have the miscentering implemented.

        Returns
        ----------
        Quantity
            Differential surface mass density profiles (ndarray, in
            astropy.units of Msun/pc/pc). Each row corresponds to a single
            cluster halo.
        """
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


#------------------------------------------------------------------------------


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



