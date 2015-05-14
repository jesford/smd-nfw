## Weak Lensing Profiles for NFW Halos

*Work in progress! I intend to convert the main c program (smd_nfw_omp.c) to python and/or cython and make it easier to use. But if the code above is helpful to you, please go ahead and use it, or, better yet, improve it and send me a pull request. Thanks!*

The c program **smd_nfw_omp.c** was used to calculate weak gravitational lensing magnification and shear profiles for the following published articles:

- "CFHTLenS: A Weak Lensing Shear Analysis of the 3D-Matched-Filter Galaxy Clusters," [J. Ford et al. (2015)](http://arxiv.org/abs/1409.3571)

- "Cluster Magnification & the Mass-Richness Relation in CFHTLenS," [J. Ford et al.](http://arxiv.org/abs/1310.2295)

This code calculates the surface mass density (Sigma) and differential surface mass density (Delta Sigma), which are probed by weak lensing magnification and shear, respectively. It can calculate these profiles for many halos at once, given their different masses and concentrations (parameterized using the scale radius r_s and concentration parameter delta_c), and their redshifts and a cosmological model (parameterized through rho_crit(z)). This is useful if you are interested in fitting a composite-halo model to an ensemble of different stacked galaxy clusters (as done in the above two papers), as opposed to simply fitting for a single average halo profile.

Optionally, this code will include the effects of cluster miscentering. Assuming that the cluster centroid offsets are described by a two-diminsional Gaussian centered on the identified center, the algorithm performs a convolution between the density profile expected for a perfectly centred halo, and this offset distribution. More details are in the above two papers. An excellent introduction to halo miscentering is given in [George et al. (2012)](http://arxiv.org/abs/1205.4262).


#### How to use this code

See the comments at the top of the file (smd_nfw_omp.c) for details on required input files, variables, output files, and compiling instructions.