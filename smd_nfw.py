"""Add the module level docstring..."""

from __future__ import absolute_import, division, print_function

import numpy as np
from astropy import units
from scipy.integrate import simps, romb, cumtrapz
#quad, dblquad

import utils


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
    sigma_nfw()
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

    def sigma_nfw(self):
        """Calculate NFW surface mass density profile.

        Generate the surface mass density profiles of each cluster halo,
        assuming a spherical NFW model. Optionally includes the effect of
        cluster miscentering offsets, if the parent object was initialized
        with offsets.

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
                print('\nERROR: f is not all real\n')#, f)

            #calculate & return centered profiles
            if singlecluster is None:
                if f.ndim == 2:
                    sigma = 2. * self._rs_dc_rcrit * f
                else:
                    rs_dc_rcrit_4D = self._rs_dc_rcrit.T.reshape(1,1, \
                                                    f.shape[2],f.shape[3])
                    sigma = 2. * rs_dc_rcrit_4D * f
                                    
            else:
                sigma = 2. * self._rs_dc_rcrit[singlecluster,0] * f

            #print('\nsigma.shape', sigma.shape)
            return sigma

        
        def _offset_sigma(self):

            #size of "x" arrays to integrate over
            numRoff = 300
            numTh = 100
            
            numRbins = self._nbins
            maxsig = self._sigmaoffset.value.max()
            roff_1D = np.linspace(0., 4.*maxsig, numRoff, endpoint=False) 
            theta_1D = np.linspace(0., 2.*np.pi, numTh, endpoint=False)
            rMpc_1D = self._rbins.value

            #reshape for broadcasting: (numTh,numRoff,numRbins)
            theta = theta_1D.reshape(numTh,1,1)
            roff = roff_1D.reshape(1,numRoff,1)
            rMpc = rMpc_1D.reshape(1,1,numRbins)
            
            r_eq13 = np.sqrt(rMpc ** 2 + roff ** 2 -
                             2. * rMpc * roff * np.cos(theta))
            
            #3D array r_eq13 -> 4D dimensionless radius (nlens)
            _set_dimensionless_radius(self, radii = r_eq13,
                                      integration = True)
                
            sigma = _centered_sigma(self)
            inner_integrand = sigma.value/(2.*np.pi)

            #integrate over theta axis: 
            sigma_of_RgivenRoff = simps(inner_integrand, x=theta_1D, axis=0)
            
            #theta is gone, now dimensions are: (numRoff,numRbins,nlens)
            sig_off_3D = self._sigmaoffset.value.reshape(1,1,self._nlens)
            roff_v2 = roff_1D.reshape(numRoff,1,1)
            PofRoff = (roff_v2/(sig_off_3D**2) *
                       np.exp(-0.5*(roff_v2 / sig_off_3D)**2))

            dbl_integrand = sigma_of_RgivenRoff * PofRoff
            
            #integrate over Roff axis (axis=0 after theta is gone):
            sigma_smoothed = simps(dbl_integrand, x=roff_1D, axis=0)
            
            #reset _x to correspond to input rbins (default)
            _set_dimensionless_radius(self)

            sigma_sm = np.array(sigma_smoothed.T) * units.solMass / units.pc**2
            
            return sigma_sm


        if self._sigmaoffset is None:
            finalsigma = _centered_sigma(self)
            errsigma = None
        elif np.abs(self._sigmaoffset).sum() == 0:
            finalsigma = _centered_sigma(self)
        else:
            finalsigma = _offset_sigma(self)
            
        print('\nThis is the new experimental python version!\n')
        return finalsigma

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


def _set_dimensionless_radius(self, radii = None, singlecluster = None,
                              integration = False):
    if radii is None:
        radii = self._rbins  #default radii
        
    #calculate x = radii / rs
    if integration == True:
        #radii is a 3D matrix (r_eq13 <-> numTh,numRoff,numRbins)
        d0, d1, d2 = radii.shape[0], radii.shape[1], radii.shape[2]

        #with each cluster's rs, we now have a 4D matrix:
        radii_4D = radii.reshape(d0, d1, d2, 1)
        rs_4D = self._rs.reshape(1, 1, 1, self._nlens)
        x = radii_4D/rs_4D
        
    elif singlecluster is None:
        #1D array of radii and clusters, reshape & broadcast
        radii_repeated = radii.reshape(1,self._nbins)
        rs_repeated = self._rs.reshape(self._nlens,1)
        x = radii_repeated/rs_repeated
    else:
        #1D array of radii and 1 cluster only
        x = radii.reshape(1,1)/self._rs[singlecluster]

    x = x.value
    
    if 0. in x:
        x[np.where(x == 0.)] = 1.e-10 #hack to avoid infs in sigma
        #print('Resetting x = zeros to 1e-10')
    
    #dimensionless radius
    self._x = x

    #set the 3 cases of dimensionless radius x
    self._x_small = np.where(self._x < 1.-1.e-6)
    self._x_big = np.where(self._x > 1.+1.e-6)
    self._x_one = np.where(np.abs(self._x-1) <= 1.e-6)
