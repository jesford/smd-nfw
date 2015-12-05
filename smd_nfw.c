#include <omp.h>
#include <stdio.h>
#include <math.h>

int main()
{
  //-------------------------------------------------------------------------
  // SURFACE MASS DENSITY: SIGMA & DELTA-SIGMA
  //  
  // by Jes Ford... modified Dec 2015
  //
  // Function for calculating the surface mass density (Sigma) and 
  // differential surface mass density (DeltaSigma), according to the NFW 
  // profile. Optionally takes miscentering into account if the first value
  // of sig_center is non-0; if sig_center is zero then the code assumes no 
  // miscentering effects.
  //
  // This program DOES NOT use parallel processing.
  // (see smd_nfw_omp.c for faster parallel code with openMP).
  // 
  // Compile (works on my mac anyway) with: 
  //      "gcc -o smd_nfw smd_nfw.c -lm"
  //
  // This program requires 2 input files: 
  //
  //                      smd_in1.dat ... single column of {r}
  //                      smd_in2.dat ... 4 columns of {rs, delta_c, rho_crit, sig_center}
  //
  //  and gives 2 output files:
  //                      sigma.dat ........  #Rbins columns by #lenses rows [Msun/pc^2]
  //                      deltasigma.dat ...  #Rbins columns by #lenses rows [Msun/pc^2]
  //
  //
  // INPUT PARAMETERS FROM FILES:
  // r ............ radii of interest [Mpc]
  // rs ........... scale radius of cluster [Mpc]
  // delta_c ...... concentration parameter (~characteristic overdensity)
  // rho_crit ..... critical energy density [Msun/pc^3] 
  //                (depends on lens z and cosmology)
  // sig_center ... Gaussian spread of miscentering [Mpc]
  //
  // All input parameters are arrays of length = # lenses, except for
  // r, which has length = # radial bins.
  //
  //-------------------------------------------------------------------------
  // CHECK FILE LENGTHS AND READ IN ALL DATA

  //open files of input {R_Mpc, r_scale, delta_c, rho_crit, sig_center}
  FILE *d1,*d2;
  d1=fopen("smd_in1.dat", "r");
  d2=fopen("smd_in2.dat", "r");

  //find number of lines in each file and assign to number of objects
  long int lines=0;
  char check[100];
  while((fgets(check, sizeof check, d1)) != NULL) lines++;
  long int nbins=lines;
  //printf("\nNumber of R bins: %d\n", nbins);
  lines=0;
  while((fgets(check, sizeof check, d2)) != NULL) lines++;
  //printf("%s\n", check);
  long int nlens=lines;
  //printf("Number of Clusters: %d\n", nlens);

  //input data arrays
  //double data[nlens][3];
  double r[nbins];
  double rs[nlens], delta_c[nlens], rho_crit[nlens], sig_center[nlens];

  //close and reopen files
  fclose(d1);
  fclose(d2);
  d1=fopen("smd_in1.dat", "r");
  d2=fopen("smd_in2.dat", "r");

  //loop indices
  long int i, j, k, m, a, b; 

  //read in {R_Mpc, r_scale, delta_c, rho_crit, sig_center}
  for(a=0; a<nbins; a++) 
  {
    fscanf(d1, "%lf", &r[a]);
  }
  for(a=0; a<nlens; a++) 
  {
    fscanf(d2, "%lf %lf %lf %lf", &rs[a], &delta_c[a], &rho_crit[a], &sig_center[a]);
    rho_crit[a] = rho_crit[a]*pow(10.,6.); //convert [Msun/pc^3] -> [Msun/(Mpc * pc^2)]
  }

  //close files
  fclose(d1);
  fclose(d2);

  //declare output arrays
  double sigma_nfw[nlens][nbins];
  double deltasigma_nfw[nlens][nbins];
  double mean_inside_sigma_nfw[nlens][nbins];

  //for "approx=0" variable comparisons
  double verysmall = pow(10.,-8.);
  long int c1,c2,c3;


  if((sig_center[0]) == 0.)
  //------------------------------------------------------------------------
  // IF ASSUMING NO MISCENTERING...
  {
    //printf("\nAssuming Perfect Centers...\n");
    double f[nlens], g[nlens], h[nlens];
    double x, bigF, part_g, firstpart1, firstpart2;

    for(j=0; j<nbins; j++)
    {
      c1=0;
      c2=0;
      c3=0;
      for(i=0; i<nlens; i++)
      {
        x = r[j]/rs[i];  //dimensionless scalar (ratio)

        if((1.-x) >= verysmall) 
	{
          c1++;
	  bigF = log((1./x)+sqrt((1./(x*x))-1.))/sqrt(1.-(x*x));  //log should be ln
	  firstpart1 = ((4./(x*x))+(2./(x*x-1.)))/sqrt(1.-(x*x));
	  part_g = firstpart1*log((1.+sqrt((1.-x)/(1.+x)))/(1.-sqrt((1.-x)/(1.+x))));
          f[i] = (1.-bigF)/(x*x-1.);
          g[i] = part_g + (4./(x*x))*log(x/2.)-(2./(x*x-1.));
          h[i] = (bigF + log(x/2.))/(x*x);
        }
        else if((x-1.) >= verysmall)
	{
          c2++;
          bigF=acos(1./x)/sqrt(x*x-1.);
          firstpart2=(8./(x*x*sqrt(x*x-1.)))+(4./pow((x*x-1.),1.5));
          part_g=firstpart2*atan(sqrt((x-1.)/(1.+x)));
          f[i] = (1.-bigF)/(x*x-1.);
          g[i] = part_g + (4./(x*x))*log(x/2.)-(2./(x*x-1.));
          h[i] = (bigF + log(x/2.))/(x*x);
        }
        else
	{
          c3++;
          f[i] = 1./3.;
          g[i] = (10./3.)+4.*log(0.5);
          h[i] = 1. + log(0.5);
        }

        //if((i == 0)) printf("\nf[i]: %lf\n", f[i]);

        sigma_nfw[i][j] = 2.*rs[i]*delta_c[i]*rho_crit[i]*f[i];

        //options below are equivalent! However both formulations break
        //down for miscentering, so calculation is different for that block.

        //deltasigma_nfw[i][j] = rs[i]*delta_c[i]*rho_crit[i]*g[i];
        mean_inside_sigma_nfw[i][j] = 4.*rs[i]*delta_c[i]*rho_crit[i]*h[i];
        deltasigma_nfw[i][j] = mean_inside_sigma_nfw[i][j] - sigma_nfw[i][j];  
     

        //if((i == 0)) printf("\ni, j, sigma_nfw[i][j]: %ld, %ld, %lf\n", i,j,sigma_nfw[i][j]);
        //if((i == 0)) printf("\nrs[i], delta_c[i], rho_crit[i]: %lf, %lf, %lf\n", rs[i],delta_c[i],rho_crit[i]);
      }
    }
    //printf("\nc1,c2,c3: %d, %d, %d\n", c1,c2,c3);
    //printf("\nsigma_nfw[i][j]: %lf\n", sigma_nfw[i][j]);

  }
  else
  //------------------------------------------------------------------------
  //IF YOU ARE TAKING MISCENTERING IN ACCOUNT...
  {
    //printf("\nMiscentered Cluster Calculations...\n");

    double f, g, h;  //f,g,h are just scalars in this section
    double x, bigF, part_g, firstpart1, firstpart2;

    long int numRp=20;   //precision of integration over r < min(Rbins)
    long int numRc=300;  //precision of integration over R_centoff
    long int numTh=100;  //precision of integration over theta

    double maxsig=0.;  //max miscentering sigma_centoff
    for(i=0; i<nlens; i++) 
      if((sig_center[i]) > maxsig) maxsig=sig_center[i];

    //printf("\nmaxsig[k]: %lf\n", maxsig);

    //Rp is a composite of the numRp linear bins interior to  
    // min(r), and midpoints of r bins (the actual measurement bins)
    double Rp[nbins+numRp-1];
    double deltaRp[nbins+numRp-1];
    for(j=0; j<(nbins+numRp-1); j++)
    {
      if(j<numRp) 
      {
        deltaRp[j]=r[0]/numRp;
        Rp[j]=(j+0.5)*deltaRp[j];
      }
      else if(j>=numRp) 
      {
        deltaRp[j]=r[(j-numRp+1)]-r[(j-numRp)];
        Rp[j]=(r[(j-numRp)]+r[(j-numRp+1)])/2.;
      }

      //printf("\nRp[j]: %lf\n", Rp[j]);
      //printf("\ndeltaRp[j]: %lf\n", deltaRp[j]);
    }

    //R_c = R_centoff array spanning range of possible miscenterings 
    // (PofRc exponential drops to essentially zero, when Rc ~ 4*sig_center)
    double Rc[numRc];
    double PofRc[numRc][nlens]; //P(R_centoff) probability of miscentering offsets
    for(k=0; k<numRc; k++)
    {
      Rc[k] = 4.*k*maxsig/numRc;
      //printf("\nRc[k]: %lf\n", Rc[k]);

      for(i=0; i<nlens; i++)
        PofRc[k][i]=(Rc[k]/(sig_center[i]*sig_center[i]))*exp(-0.5*pow((Rc[k]/sig_center[i]),2.));
    } 

    double theta[numTh];
    for(m=0; m<numTh; m++) 
      theta[m] = 2.*M_PI*(m+1.)/numTh;

    double dRc,dtheta;
    dRc = (Rc[numRc-1]-Rc[0])/(numRc-1.);        //differential spacing
    dtheta=(theta[numTh-1]-theta[0])/(numTh-1.); //differential angle

    //divide into nt=8 threads and parallel process...
    long int nt=8; //number of threads  ????????????????
    long int tid;  //thread id number: 0 to (nt-1)

    double sigmaof_rgivenRc=0.;

    double sigma_smoothed[nlens][nbins][nt];
    double sigma_smoothed_Rp[nlens][nbins+numRp-1][nt];
    double mean_inside_sigma_smoothed[nlens][nbins];

    //Force all array values to zero
    //(NECESSARY since they are ADDED TO below)
    for(i=0; i<nlens; i++)
    {
      for(j=0; j<nbins; j++)
      {
        mean_inside_sigma_smoothed[i][j]=0.;
        for(tid=0; tid<nt; tid++) sigma_smoothed[i][j][tid]=0.;
      }
      for(j=0; j<(nbins+numRp-1); j++)
        for(tid=0; tid<nt; tid++) sigma_smoothed_Rp[i][j][tid]=0.;
    }


    //printf("\nBeginning PARALLEL processing...\n");

    //#pragma omp parallel for private(tid, i, j, k, m, x, bigF, sigmaof_rgivenRc) num_threads(nt) schedule(dynamic)
    for(i=0; i<nlens; i++)  //lens loop, threads share execution of this
    {
      //tid=omp_get_thread_num();
      tid=0;

        //----------------------------------------
        // SIGMA CALCULATION

	for(j=0; j<nbins; j++)  //R bins loop
	{

          //if(i==0) printf("\nj = %d", j);
          for(k=0; k<numRc; k++)  //R_centoff loop
	  {
            sigmaof_rgivenRc=0.;

            for(m=0; m<numTh; m++)  //theta loop
	    {
              // x = r_offset/r_scale [note: instead of sqrt(abs(...)) doing 4th-root of square]
              // THIS x is for sigma (different than below)
              x = pow(pow((r[j]*r[j] + (Rc[k]*Rc[k]) - 2.*r[j]*Rc[k]*cos(theta[m])),2.),0.25)/rs[i];


              if((1.-x) >= verysmall) 
              {
	        bigF = log((1./x)+sqrt((1./(x*x))-1.))/sqrt(1.-(x*x));  //log should be ln
                f = (1.-bigF)/(x*x-1.);
              }
              else if((x-1.) >= verysmall)
	      {
                bigF=acos(1./x)/sqrt(x*x-1.);
                f = (1.-bigF)/(x*x-1.);
              }
              else f = 1./3.;

              //CONVOLUTION: INTEGRAL OVER THETA
              //EQ 7 (Johnston et al. 2007)
              sigmaof_rgivenRc += (2.*rs[i]*delta_c[i])*(rho_crit[i]*f*(dtheta/2./M_PI));
              //#pragma omp flush

	    } //end m loop

            //INTEGRAL OVER R_centoff
            //EQ 9 (Johnston et al. 2007, with his Rs -> my Rc)
            sigma_smoothed[i][j][tid] += sigmaof_rgivenRc*PofRc[k][i]*dRc;

            //if((i == 349) && (k == numRc-1) && (j == 0)) printf("\nsigma_smoothed[i][j][tid]: %lf\n", sigma_smoothed[i][j][tid]);

          } //end k loop
        } //end j loop

        //----------------------------------------
        // MEAN-SIGMA(<r) CALCULATION
        for(j=0; j<(nbins+numRp-1); j++)  //Rp (extended R) bins loop
        {
          //if(i==0) printf("\njp = %d", jp);
          for(k=0; k<numRc; k++)  //R_centoff loop
	  {
            sigmaof_rgivenRc=0.;

            for(m=0; m<numTh; m++)  //theta loop
	    {
              //#pragma omp flush
              // x = r_offset/r_scale [note: instead of sqrt(abs(...)) doing 4th-root of square]
              // THIS x is for calculating deltasigma from mean(sigma(<r))
              x = pow(pow((Rp[j]*Rp[j] + (Rc[k]*Rc[k]) - 2.*Rp[j]*Rc[k]*cos(theta[m])),2.),0.25)/rs[i]; 

              if((1.-x) >= verysmall) 
              {
	        bigF = log((1./x)+sqrt((1./(x*x))-1.))/sqrt(1.-(x*x));  //log should be ln
                f = (1.-bigF)/(x*x-1.);
              }
              else if((x-1.) >= verysmall)
	      {
                bigF=acos(1./x)/sqrt(x*x-1.);
                f = (1.-bigF)/(x*x-1.);
              }
              else f = 1./3.;


              //CONVOLUTION: INTEGRAL OVER THETA
              //EQ 7 (Johnston et al. 2007)
              sigmaof_rgivenRc += (2.*rs[i]*delta_c[i])*(rho_crit[i]*f*(dtheta/2./M_PI));

	    } //end m loop

            //INTEGRAL OVER R_centoff
            //EQ 9 (Johnston et al. 2007, with his Rs -> my Rc)
            sigma_smoothed_Rp[i][j][tid] += sigmaof_rgivenRc*PofRc[k][i]*dRc;


          } //end k loop
        } //end j loop

        //----------------------------------------

	//printf("\ni,tid,sigma_smoothed[0][0][tid]: %ld %ld %lf \n", i, tid, sigma_smoothed[0][0][tid]);

    } //end i loop
    //printf("\nFinished PARALLEL processing.\n");

    //for(tid=0; tid<nt; tid++) printf("\nsigma_smoothed (i=349,j=0): %lf", sigma_smoothed[349][0][tid]);
    //for(tid=0; tid<nt; tid++) printf("\nsigma_smoothed_Rp (i=349,j=0): %lf", sigma_smoothed_Rp[349][0][tid]);


    //COMBINE THREADS
    for(i=0; i<nlens; i++)
    {
      for(j=0; j<nbins; j++)
      {
        sigma_nfw[i][j] = 0.;
	deltasigma_nfw[i][j] = 0.;
        mean_inside_sigma_nfw[i][j] = 0.;

        for(tid=0; tid<nt; tid++)
	{
          //FINAL output (smoothed sigma)
          sigma_nfw[i][j] += sigma_smoothed[i][j][tid];

          //INTEGRAL OVER INSIDE (r<R)
          //EQ 8 (George et al. 2012)
          for(k=0; k<(j+numRp); k++)
            mean_inside_sigma_smoothed[i][j] += (sigma_smoothed_Rp[i][k][tid]*Rp[k]*deltaRp[k])*(2./(r[j]*r[j]));

	}

        //FINAL output (smoothed deltasigma)
        deltasigma_nfw[i][j] = mean_inside_sigma_smoothed[i][j] - sigma_nfw[i][j]; 

      } //end j loop
    } //end i loop

  }  //END MISCENTERING OPTION

  //---------------------------------------------------------------------
  // PRINT OUTPUT

  //create and open files for writing
  FILE *out1, *out2;
  out1=fopen("sigma.dat", "w");
  out2=fopen("deltasigma.dat", "w");

  //print output data
  for(a=0; a<nlens; a++) 
  {
    for(b=0; b<nbins; b++) 
    {
      fprintf(out1,"%lf  ", sigma_nfw[a][b]);
      fprintf(out2,"%lf  ", deltasigma_nfw[a][b]);
    }
    fprintf(out1,"\n");
    fprintf(out2,"\n");
  }


  //close files
  fclose(out1);
  fclose(out2);

  //the end
  //printf("\nsmd_nfw_omp FINISHED.\n");
  return 0;
}
