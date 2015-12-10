## Weak Lensing Profiles for NFW Halos

*Work in progress! I am converting the main c program (smd_nfw_omp.c) to python/cython and make it easier to use. See details below. But if the code above is helpful to you, please go ahead and use it, or, better yet, improve it and send me a pull request. Thanks!*

The c program **smd_nfw_omp.c** was used to calculate weak gravitational lensing magnification and shear profiles for the following published articles:

- "CFHTLenS: A Weak Lensing Shear Analysis of the 3D-Matched-Filter Galaxy Clusters," [J. Ford et al. (2015)](http://arxiv.org/abs/1409.3571)

- "Cluster Magnification & the Mass-Richness Relation in CFHTLenS," [J. Ford et al. (2014)](http://arxiv.org/abs/1310.2295)

The second c program **smd_nfw.c** is identical, except that it does *not* use parallel processing, so may be easier to compile and run on some machines. It is also currently being used to generate NFW profiles for the weak lensing project [cluster-lensing](https://github.com/jesford/cluster-lensing), which provides a much nicer python interface to this c program, through the ClusterEnsemble() object.

**smd_nfw.py** is a new pure python module, which can currently calculate NFW profiles for perfectly centered clusters only. In the [cluster-lensing](https://github.com/jesford/cluster-lensing) project, you can now use the "pure_python" option, to avoid dealing with the c code (hooray). *I'll add the miscentered NFW python calculation soon.*

These codes calculate the NFW surface mass density (Sigma) and differential surface mass density (DeltaSigma), which are probed by weak lensing magnification and shear, respectively. It can calculate these profiles for many halos at once, given their different masses and concentrations (parameterized using the scale radius r<sub>s</sub> and concentration parameter delta<sub>c</sub>), and their redshifts and a cosmological model (parameterized through rho<sub>crit</sub>(z)). This is useful if you are interested in fitting a composite-halo model to an ensemble of different stacked galaxy clusters (as done in the above two papers), as opposed to simply fitting for a single average halo profile.

Optionally, this code will include the effects of cluster miscentering. Assuming that the cluster centroid offsets are described by a two-diminsional Gaussian centered on the identified center, the algorithm performs a convolution between the density profile expected for a perfectly centred halo, and this offset distribution. More details are in the above two papers. An excellent introduction to halo miscentering is given in [George et al. (2012)](http://arxiv.org/abs/1205.4262).


#### How to use this code

See the comments at the top of the file (smd_nfw_omp.c) for details on required input files, variables, output files, and compiling instructions. **I recommended using the python class clusters.ClusterEnsemble() to interface to either of the smd_nfw.c or smd_nfw.py codes. See this repository:** [cluster-lensing](https://github.com/jesford/cluster-lensing)
