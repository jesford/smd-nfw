## Weak Lensing Profiles for NFW Halos

*Note: this repository of C programs is being kept for historical reference. I highly recommend using the improved versions of these algorithms instead, which are available in the Python* [cluster-lensing](https://github.com/jesford/cluster-lensing) *project. The Python implementation has a much nicer API and allows you to easily build up a DataFrame of galaxy cluster cluster properties and NFW profiles, with and without miscentering effects, through the ClusterEnsemble object. See the* [project website](https://github.com/jesford/cluster-lensing).

The C program **smd_nfw_omp.c** was used to calculate weak gravitational lensing magnification and shear profiles for the following published articles:

- "CFHTLenS: A Weak Lensing Shear Analysis of the 3D-Matched-Filter Galaxy Clusters," [J. Ford et al. (2015)](http://arxiv.org/abs/1409.3571)

- "Cluster Magnification & the Mass-Richness Relation in CFHTLenS," [J. Ford et al. (2014)](http://arxiv.org/abs/1310.2295)

The second C program **smd_nfw.c** is identical, except that it does *not* use parallel processing, so may be easier to compile and run on some machines. See the comments at the top of the files for details on required input files, variables, output files, and compiling instructions.

These codes calculate the NFW surface mass density (Sigma) and differential surface mass density (DeltaSigma), which are probed by weak lensing magnification and shear, respectively. It can calculate these profiles for many halos at once, given their different masses and concentrations (parameterized using the scale radius r<sub>s</sub> and concentration parameter delta<sub>c</sub>), and their redshifts and a cosmological model (parameterized through rho<sub>crit</sub>(z)). This is useful if you are interested in fitting a composite-halo model to an ensemble of different stacked galaxy clusters (as done in the above two papers), as opposed to simply fitting for a single average halo profile.

Optionally, this code will include the effects of cluster miscentering. Assuming that the cluster centroid offsets are described by a two-diminsional Gaussian centered on the identified center, the algorithm performs a convolution between the density profile expected for a perfectly centred halo, and this offset distribution. More details are in the above two papers. An excellent introduction to halo miscentering is given in [George et al. (2012)](http://arxiv.org/abs/1205.4262).
