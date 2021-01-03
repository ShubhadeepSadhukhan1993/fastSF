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

Turbulence is a complex phenomenon in fluid dynamics involving nonlinear interactions between multiple scales. Structure functions are popular diagnostics in the study of statistical properties properties of turbulent flows [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]. Some of the earlier works comprising of such analysis are those of @Gotoh:PF2002, @Kaneda:PF2003, and @Ishihara:ARFM2009 for three-dimensional (3D) hydrodynamic turbulence; @Yeung:PF2005 and @Ray:NJP2008 for passive scalar turbulence; @Biferale:NJP2004 for two-dimensional (2D) hydrodynamic turbulence; and @Kunnen:PRE2008, @Kaczorowski:JFM2013, and @Bhattacharya:PF2019 for turbulent thermal convection. Structure functions are two-point statistical quantities; thus, an accurate computation of these quantities requires averaging over many points. However, incorporation of a large number of points makes the computations very expensive and challenging. Therefore, we require an efficient parallel code for accurate computation of structure functions. In this paper, we describe the design and validation of the results of ``fastSF``, a parallel code to compute the structure functions for a given velocity or scalar field. 

``fastSF`` is a C++ application for computing the structure functions of scalar and vector fields on Cartesian grids of a 2D or 3D periodic box, stored as HDF5 files. The code employs MPI (Message Passing Interface) parallelization with equal load distribution and vectorization for efficiency on SIMD architectures. The user can select the range of the orders of the structure functions to be computed and the computed structure functions are written to HDF5 files that can be further processed by the user.

We are not aware of any other open soure or commercial packages for computing structure functions; prior studies have relied on in-house software that was never publicly released.  As an open source package, `fastSF` provides a standard high-performance implementation and thus facilitates wider use of structure functions.

``fastSF`` uses MPI [@Pacheco:book:PP] for parallelism, HDF5 [@HDF5_web] via H5SI [@H5SI_web] for reading gridded field data and writing structure functions, as well as blitz++ [@Blitz_web] for vectorized computation and yaml-cpp [@YAML_web] for reading control parameters.
 In the next section, we will briefly explain the velocity and scalar structure functions in turbulent flows.



# Velocity and scalar structure functions

We denote the velocity and scalar fields using $\boldsymbol{u}$ and $\theta$  respectively. The velocity difference between any two points $\boldsymbol{r}$ and $\boldsymbol{r}+\boldsymbol{l}$ is $\delta \boldsymbol{u} = \boldsymbol{u(r}+ \boldsymbol{l)}-\boldsymbol{u(r)}$. The difference in the parallel components of the velocity field along $\boldsymbol{l}$ is $\delta u_\parallel=\delta \boldsymbol{u}\cdot \hat{\boldsymbol{l}}$.  The corresponding difference in the perpendicular component is $\delta u_\perp= |\delta \boldsymbol{u} - \delta u_\parallel \hat{\boldsymbol{l}}|$. Assuming statistical homogeneity, we define the longitudinal velocity structure functions of order $q$ as
$$ S_q^{u_\parallel}(\boldsymbol{l}) = \langle (\delta u_\parallel)^q \rangle = \langle [\{\boldsymbol{u(r+l)}-\boldsymbol{u(r)}\}\cdot \hat{\boldsymbol{l}}]^q \rangle, \quad \quad (1)$$ 
and the transverse velocity structure functions of order 
$q$ as 
$$ S_q^{u_\perp}(\boldsymbol{l}) = \langle (\delta u_\perp)^q \rangle = \langle |\delta \boldsymbol{u} - \delta u_\parallel \hat{\boldsymbol{l}}|^q \rangle. \quad \quad (2)$$ 
Here, $\langle \cdot \rangle$ denotes spatial averaging. Similarly, we can define the scalar structure functions for the scalar field as 
$$ S_q^\theta(\boldsymbol{l}) = \langle (\delta \theta)^q\rangle = \langle [\theta (\boldsymbol{r+l}) - \theta(\boldsymbol{r})]^q \rangle. \quad \quad(3)$$

For isotropic turbulence (in addition to being homogeneous), the structure functions become functions of $l$, where $l=|\boldsymbol{l}|$. The second-order velocity structure function $S_q^{u_{\parallel}}(l)$ provides an estimate for the energy in the eddies of size $l$ or less [@Davidson:book:Turbulence]. 

![For 3D homogeneous isotropic turbulence: plots of the negative of normalized third, fifth and seventh-order longitudinal velocity structure functions vs. $l$. The negative of the normalized third-order structure function is close to $4/5$ (dashed line) in the inertial range. \label{SF_Hydro}](docs/figs/SF_hydro.png)

For 3D incompressible hydrodynamic turbulence with homegeneity and isotropy, the third-order longitudinal velocity structure function in the inertial range (scales lying between the large-scale forcing regime and the small-scale dissipation regime) is given by [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]
$$S_3^{u_\parallel}(l) = -\frac{4}{5} \epsilon l \sim -l, \quad \quad (4)$$
where $\epsilon$ is the viscous dissipation rate. 
For an arbitrary order $q$, @She:PRL1994 proposed that the longitudinal structure functions scale as $S_q^{u_\parallel} (l) \sim \zeta_q$, where the exponent $\zeta_q$ is given by
$$\zeta_q = \frac{q}{9} + 2\left(1 - \left( \frac{2}{3} \right)^{q/3} \right). \quad \quad (5)$$


Figure~\ref{SF_Hydro} exhibits the plots of the negative of the normalized 3rd, 5th, and 7th-order longitudinal velocity structure functions computed using the simulation data of 3D hydrodynamic turbulence [@Sadhukhan:PRF2019]. The structure functions are normalized by $(\epsilon l)^{\zeta_q}$, where $\zeta_q$ is given by Eq. (5).  In the inertial range (0.2 < l < 0.7), the normalized third-order longitudinal velocity structure function is fairly close to $4/5$ (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of l, thus exhibiting consistency with She-Leveque's model.




In the next section, we provide a brief description of the code.

# Design of the Code
In this section, we present a sketch of the structure function computation for the velocity structure functions.  We employ vectorization and loops over $\boldsymbol{l}$, thus requiring three loops for 3D fields and two loops for 2D fields. In the following, we provide the algorithm for structure function computation for a 2D velocity field.

**Pseudo-code**

*Data*: Velocity field $\boldsymbol{u}$ in domain $(L_x, L_z)$; number of processors $P$.

*Procedure*:
 
* Divide $\boldsymbol{l}$'s among various processors. The process of data division among the processors  has been described later in this section. 
 
* For every processor:
     
    * for $\boldsymbol{l}= (l_x,l_z)$ assigned to the processor:
        
        * Compute $\delta \boldsymbol{u}(l_x,l_z)$ by taking the difference between two points with the same indices in pink and green subdomains as shown in Fig. \ref{Schematic}. This feature enables vectorized subtraction operation.
        
        * $\delta u_{\parallel}(l_x,l_z) = \delta \boldsymbol{u} \cdot \hat{\boldsymbol{l}}$ (Vectorized). 
        
        * $\delta u_{\perp}(l_x,l_z) = |\delta \boldsymbol{u} - \delta u_{\parallel} \hat{\boldsymbol{l}}$| (Vectorized). 
        
        * for order $q$:
        
            * $S_q^{u_{\parallel}}(l_x,l_z) =$ Average of $\delta u_{\parallel}^q$ (Vectorized).
            
            * $S_q^{u_{\perp}}(l_x,l_z) =$ Average of $\delta u_{\perp}^q$ (Vectorized).
            
            * Send the values of $S_q^{u_{\parallel}}(l_x,l_z)$, $S_q^{u_{\perp}}(l_x,l_z)$, $q$, $l_x$, and $l_z$ to the root process.
            
* The root process stores $S_q^{u_{\parallel}} (l_x, l_z)$ and $S_q^{u_{\perp}} (l_x, l_z)$.
            
* Stop

![The velocity difference $\delta \boldsymbol{u}(\boldsymbol{l})$ is computed by taking the difference between two points with the same indices in the pink and the green subdomains. For example, $\boldsymbol{u}(\boldsymbol{l}) - \boldsymbol{u}(0,0) = \boldsymbol{u}_B - \boldsymbol{u}_A$, where $B$ and $A$ are the origins of the green and the pink subdomains. This feature enables vecotrization of the computation. \label{Schematic}](docs/figs/Schematic.png)

Since $S_q^u(\boldsymbol{l})$ is important for intermediate scales (inertial range) only, we vary $\boldsymbol{l}$ upto half the domain size, that is, upto ($L_x/2, L_z/2$), to save computational cost. The $\boldsymbol{l}$'s are divided among MPI processors along $x$ and $z$ directions. Each MPI processor computes the structure functions for the points assigned to it and has access to the entire input data. 
After computing the structure function for a given $\boldsymbol{l}$, each processor communicates the result to the root process, which stores the $S_q^{u_\parallel}(\boldsymbol{l})$ and $S_q^{u_\perp}(\boldsymbol{l})$ arrays.

It is clear from Fig. \ref{Schematic} that the sizes of the pink or green subdomains are $(L_x-l_x)(L_z-l_z)$, which are function of $\boldsymbol{l}$'s.  This function decreases with increasing $\boldsymbol{l}$ leading to larger computational costs for small $l$ and less cost of larger $l$.   Hence, a straightforward division of the domain among the processors along $x$ and $z$ directions will lead to a load imbalance.   Therefore, we assign both large and small $\boldsymbol{l}$'s to each processor to achieve equal load distribution. We illustrate the above idea  using the following example.

Consider a one-dimensional domain of size $L=15$, for which the possible $l$'s are
$$l=\{0, 1, 2, 3 ... 15\}.$$ 
We need to compute the structure functions for $l$ ranging from 0 to 7. We divide the task among four processors, with two $l$'s assigned to each processor. The following distribution of $l$'s ensures equal load distribution:
$$\mbox{Processor 0: } \quad l=\{0,7\}, \quad \sum(L-l)=(15-0)+(15-7) = 23,$$
$$\mbox{Processor 1: } \quad l=\{1, 6\}, \quad \sum(L-l)=(15-1)+(15-6) = 23,$$
$$\mbox{Processor 2: } \quad l=\{2,5\}, \quad \sum(L-l)=(15-2)+(15-5) = 23,$$
$$\mbox{Processor 3: } \quad l=\{3, 4\}, \quad \sum(L-l)=(15-3)+(15-4) = 23.$$
Similarly, if two processors are used, then the following distribution results in load balance. 
$$\mbox{Processor 0: } \quad l=\{0, 7, 2, 5\},$$
$$\mbox{Processor 1: } \quad l=\{1, 6, 3, 4\}.$$
 This idea of load distribution has been implemented in our program and has been extended to higher dimensions. 

Note that for 2D, $l_x>0$, but $l_z$ can take both positive and negative values. However, for isotropic turbulence, the structure functions for $+l_z$ and $-l_z$ are statistically equal. Therefore, in our computations, we restrict  to $l_x>0$, $l_z>0$. For anisotropic turbulence, not discussed here, the structure functions will depend on $(l_x,l_z)$ rather than $l$; our code will be extended to such systems in future. 

For 3D turbulence, the structure functions will depend on $(l_x,l_y,l_z)$. We divide the tasks among processors over $l_x$ and $l_y$ as done  above for 2D turbulence. The aforementioned algorithm can be easily extended to the 3D case. We employ a similar method for the computation of scalar structure functions as well.
 
In the next section, we discuss the scaling of our code.

# Scaling of `fastSF`

`fastSF` is scalable over many processors due to vectorization and equal load distribution. We demonstrate the scaling of `fastSF` for the third-order longitudinal structure function for an idealized velocity field on a $128^3$ grid.  For our computation we employ a maximum of 1024 cores. We take the velocity field as
$$\boldsymbol{u} = 
\begin{bmatrix} 
x \\ y \\z
\end{bmatrix}.$$
We perform four runs on a Cray XC40 system (Shaheen II of KAUST) for this problem using a total of 16, 64, 256, and 1024 cores. We used 16 cores per node for each run. In Fig. \ref{Scaling}, we plot the inverse of time taken in seconds versus the number of cores. The best fit curve for these data points yields
$$T^{-1} \sim p^{0.986 \pm 0.002},$$
Thus, the data-points follow $T^{-1} \sim p$ curve to a good approximation. Hence, we conclude that our code exhibits good scaling. 

![Scaling of `fastSF` for the computation of third-order longitudinal velocity structure function using 16, 64, 256, and 1024 processors of Shaheen II. All the runs were conducted on a $128^3$ grid.  We observe a linear scaling. \label{Scaling}](docs/figs/SF_scaling.png)


# Conclusions

This paper provides a brief description of ``fastSF``, an efficient parallel C++ code that computes structure functions for given velocity and scalar fields. This code is shown to be scalable over many processors. An earlier version of the code was used by @Bhattacharya:PF2019 for analyzing the structure functions of turbulent convection.  We believe that ``fastSF`` will be useful to turbulence community as it facilitates wider use of structure functions.  


# Acknowledgements

We thank Roshan Samuel, Anando Chatterjee, Soumyadeep Chatterjee, and Manohar Sharma for helpful discussions during the development of ``fastSF``. We are grateful to Jed Brown, Ilja Honkonen, and Chris Green for a careful review of our work and their useful suggestions. Our computations were performed on Shaheen II at KAUST supercomputing laboratory, Saudi Arabia, under the project k1416. 

---

# References


