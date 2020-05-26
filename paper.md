---
title: 'fastSF: A parallel code for computing the structure functions of turbulence'

tags:
  - C++
  - structure functions
  - turbulence
  - fluid dynamics

authors:
  - name: Shubhadeep Sadhukhan
    orcid: 0000-0002-7278-5041
    affiliation: 1
  - name: Shashwat Bhattacharya
    orcid: 0000-0001-7462-7680
    affiliation: 2
  - name: Mahendra K. Verma
    orcid: 0000-0002-3380-4561
    affiliation: 1
  

affiliations:
 - name: Department of Physics, Indian Institute of Technology Kanpur, Kanpur 208016, India
   index: 1
 - name: Department of Mechanical Engineering, Indian Institute of Technology Kanpur 208016, India
   index: 2

date: 18 February 2020

bibliography: paper.bib

---

# Summary

Turbulence is a complex phenomenon in fluid dynamics involving nonlinear interactions between multiple scales. Structure function is a popular diagnostics tool to study the statistical properties of turbulent flows [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]. Some of the earlier works comprising of such analysis are those of @Gotoh:PF2002, @Kaneda:PF2003, and @Ishihara:ARFM2009 for three-dimensional (3D) hydrodynamic turbulence; @Yeung:PF2005 and @Ray:NJP2008 for passive scalar turbulence; @Biferale:NJP2004 for two-dimensional (2D) hydrodynamic turbulence; and @Kunnen:PRE2008, @Kaczorowski:JFM2013, and @Bhattacharya:PF2019 for turbulent thermal convection. Structure functions are two-point statistical quantities; thus, an accurate computation of these quantities requires averaging over many points. However, incorporation of a large number of points makes the computations very expensive and challenging. Therefore, we require an efficient parallel code for accurate computation of structure functions. In this paper, we describe the design and validation of the results of ``fastSF``, a parallel code to compute the structure functions for a given velocity or scalar field. 

 ``fastSF``, written in C++, is a fast and efficient code that uses vectorization for computing the structure functions. The code employs MPI (Message Passing Interface) parallelization with equal load distribution. The user has a choice on the type (scalar or vector) and the dimensions of the fields to be read by the code, and the range of the orders of the structure functions to be computed. The code writes the computed structure functions to `hdf5` files that can be further processed by the user.
 
 The code uses the following libraries:
 
1.   blitz++ (version 1.0.2)
2.   hdf5 (version 1.8.20)
3.   h5si (version 1.1.1)
4.   yaml-cpp (version 0.3.0)
5.   mpich (version 3.3.2)
6.   cmake (version 3.16.4) 


In the next section, we will briefly discuss the performance and scaling of ``fastSF`` in a Cray XC40 system.



# Performance and scaling of `fastSF`

Typical structure function computations in literature involve calculation of the velocity or scalar difference using loops over two points. These computations require six nested `for` loops for 3D fields that makes the computations very expensive for large grids. `fastSF` is scalable over many processors due to vectorization and equal load distribution. In our code, we employ vectorization and loops over only one point, thus requiring three loops instead of six for 3D fields. The new algorithm enhances the performance approximately 20 times over the earlier schemes due to vectorization. Further, we ensure equal load distribution among processors to enhance scalabality of the code over many processors. Please refer to the code and the documentation for details.

We demonstrate the scaling of `fastSF` for the third-order longitudinal structure function for an idealized velocity field on a $128^3$ grid.  For our computation we employ a maximum of 1024 processors. We take the velocity field as
$$\mathbf{u} = 
\begin{bmatrix} 
x \\ y \\z
\end{bmatrix}.$$
We perform four runs on a Cray XC40 system (Shaheen II of KAUST) for this problem using 16, 64, 256, and 1024 processors. In Fig. \ref{Scaling}, we plot the inverse of time taken in seconds versus the number of processors. The best fit curve for these data points yields
$$T^{-1} \sim p^{0.986 \pm 0.002},$$
Thus, the data-points follow $T^{-1} \sim p$ curve to a good approximation. Thus, we conclude that our code exhibits strong scaling. 

![Scaling of `fastSF` for the computation of third-order longitudinal velocity structure function using 16, 64, 256, and 1024 processors of Shaheen II. All the runs were conducted on a $128^3$ grid.  We observe a linear scaling. \label{Scaling}](SF_scaling.png)

 


# Conclusions

This paper provides a brief description ``fastSF``, an efficient parallel C++ code that computes structure functions for given velocity and scalar fields. This code is scalable over many processors. An earlier version of the code was used by @Bhattacharya:PF2019 for analyzing the structure functions of turbulent convection. We are currently using this code to investigate the structure functions of two-dimensional turbulence with large-scale forcing. We believe that ``fastSF`` will be useful to turbulence community.  


# Acknowledgements

We thank R. Samuel, A. Chatterjee, S. Chatterjee, and M. Sharma for useful discussions during the development of ``fastSF``. Our computations were performed on Shaheen II at KAUST supercomputing laboratory, Saudi Arabia, under the project k1416. 

---

# References


