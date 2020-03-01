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
 
1.   blitz++ (version xx)
2.  hdf5
3. hfsi
4. yaml
5. mpich
6. cmake


In the next section, we will briefly define the velocity and the scalar structure functions in turbulence.

# Velocity and scalar structure functions

We denote the velocity and scalar fields using $\mathbf{u}$ and $\theta$  respectively. The velocity difference between any two points $\mathbf{r}$ and $\mathbf{r+l}$ is $\delta \mathbf{u} = \mathbf{u(r+l)}-\mathbf{u(r)}$. The difference in the parallel components of the velocity field along $\mathbf{l}$ is $\delta u_\parallel=\delta \mathbf{u}\cdot \hat{\mathbf{l}}$.  The corresponding difference in the perpendicular component is $\delta u_\perp= |\delta \mathbf{u} - \delta u_\parallel \hat{\mathbf{l}}|$. Assuming statistical homogeneity, we define the longitudinal velocity structure functions of order $q$ as
$$ S_q^{u_\parallel}(\mathbf{l}) = \langle \delta u_\parallel^q \rangle = \langle [\{\mathbf{u(r+l)}-\mathbf{u(r)}\}\cdot \hat{\mathbf{l}}]^q \rangle, \quad \quad (1)$$ 
and the transverse velocity structure functions of order 
$q$ as 
$$ S_q^{u_\perp}(\mathbf{l}) = \langle \delta u_\perp^q \rangle = \langle |\delta \mathbf{u} - \delta u_\parallel \hat{\mathbf{l}}|^q \rangle. \quad \quad (2)$$ 
Here, $\langle \cdot \rangle$ denotes ensemble averaging. Similarly, we can define the scalar structure functions for the scalar field as 
$$ S_q^\theta(\mathbf{l}) = \langle \delta \theta^q\rangle = \langle [\theta (\mathbf{r+l}) - \theta(\mathbf{r})]^q \rangle. \quad \quad(3)$$

For sotropic turbulence (in addition to being homogeneous), the structure functions become functions of $l$, where $l=|\mathbf{l}|$. The second-order velocity structure function $S_q^{u_{\parallel}}(l)$ provides an estimate for the energy in the eddies of size $l$ or less [@Davidson:book:Turbulence]. 

In the next section, we provide a brief description of the code.

# Design of the Code
First we present a sketch of the structure function computaion for the velocity structure functions.  Typical structure function computations [Eqs (1-3)] in literature involve calculation of the velocity difference using loops over $\mathbf{r}$ and $\mathbf{l}$. These computations require six nested `for` loops for 3D fields that makes the computations very expensive for large grids. In our code, we employ vectorization and loops over only $\mathbf{l}$, thus requiring three loops instead of six for 3D fields. The new algorithm enhances the performance approximately 20 times over the earlier schemes due to vectorization.  In the following, we provide the algorithm for structure function computation for a 2D velocity field.

**Pseudo-code**

*Data*: Velocity field $\mathbf{u}$ in domain $(L_x, L_z)$; number of processors $P$.

*Procedure*:
 
* Divide $\mathbf{l}$'s among various processors. The process of data division among the processors  has been described later in this section. 
 
* For every processor:
     
    * for $\mathbf{l}= (l_x,l_z)$ assigned to the processor:
        
        * Compute $\delta \mathbf{u}(l_x,l_z)$ by taking the difference between two points with the same indices in pink and green subdomains as shown in Fig. \ref{Schematic}. This feature enables vectorized subtraction operation.
        
        * $\delta u_{\parallel}(l_x,l_z) = \delta \mathbf{u} \cdot \hat{\mathbf{l}}$ (Vectorized). 
        
        * $\delta u_{\perp}(l_x,l_z) = |\delta \mathbf{u} - \delta u_{\parallel} \hat{\mathbf{l}}$| (Vectorized). 
        
        * for order $q$:
        
            * $S_q^{u_{\parallel}}(l_x,l_z) =$ Average of $\delta u_{\parallel}^q$ (Vectorized).
            
            * $S_q^{u_{\perp}}(l_x,l_z) =$ Average of $\delta u_{\perp}^q$ (Vectorized).
            
            * Send the values of $S_q^{u_{\parallel}}(l_x,l_z)$, $S_q^{u_{\perp}}(l_x,l_z)$, $q$, $l_x$, and $l_z$ to the master processor. 
            
* The master processor stores $S_q^{u_{\parallel}} (l_x, l_z)$ and $S_q^{u_{\perp}} (l_x, l_z)$.
            
* Stop

![The velocity difference $\delta \mathbf{u}(\mathbf{l})$ is computed by taking the difference between two points with the same indices in the pink and the green subdomains. For example, $\mathbf{u}(\mathbf{l}) - \mathbf{u}(0,0) = \mathbf{u}_B - \mathbf{u}_A$, where $B$ and $A$ are the origins of the green and the pink subdomains. This feature enables vecotrization of the computation. \label{Schematic}](Schematic.png)

Since $S_q^u(\mathbf{l})$ is important for intermediate scales (inertial range) only, we vary $\mathbf{l}$ upto half the domain size, that is, upto ($L_x/2, L_z/2$), to save computational cost. The $\mathbf{l}$'s are divided among MPI processors along $x$ and $z$ directions. Each MPI processor computes the structure functions for the points assigned to it and has access to the entire input data. Thus, we save considerable time that would otherwise be spent on communication between the processors during the calculation of the velocity or the scalar difference. After computing the structure function for a given $\mathbf{l}$, each processor communicates the result to the master processor, which stores the $S_q^{u_\parallel}(\mathbf{l})$, $S_q^{u_\perp}(\mathbf{l})$ and $S_q^{\theta}(\mathbf{l})$ arrays.

It is clear from Fig. \ref{Schematic} that the sizes of the pink or green subdomains are $(L_x-l_x)(L_z-l_z)$, which is function of $\mathbf{l}$'s.  This function decreases with increasing $\mathbf{l}$ leading to larger computational costs for small $l$ and less cost of larger $l$.   Hence, a straightforward division of the domain among the processors along $x$ and $z$ directions will lead to a load imbalance.   Therefore, we assign both large and small $\mathbf{l}$'s to each processor to achieve equal load distribution. We illustrate the above idea  using the following example.

Consider a one-dimensional domain of size $L=15$, for which the possible $l$'s are
$$l=\{0, 1, 2, 3 ... 15\}.$$ 
We need to compute the structure functions for $l$ ranging from 0 to 7. We divide the task among four processors, with 2 $l$'s assigned to each processor. The following distribution of $l$'s ensures equal load distribution:
$$\mbox{Processor 0: } \quad l=\{0,7\}, \quad \sum(L-l)=(15-0)+(15-7) = 23,$$
$$\mbox{Processor 1: } \quad l=\{1, 6\}, \quad \sum(L-l)=(15-1)+(15-6) = 23,$$
$$\mbox{Processor 2: } \quad l=\{2,5\}, \quad \sum(L-l)=(15-2)+(15-5) = 23,$$
$$\mbox{Processor 3: } \quad l=\{3, 4\}, \quad \sum(L-l)=(15-3)+(15-4) = 23.$$
Similarly, if two processors are used, then the following distribution results in load balance. 
$$\mbox{Processor 0: } \quad l=\{0, 7, 2, 5\},$$
$$\mbox{Processor 1: } \quad l=\{1, 6, 3, 4\}.$$
 This idea of load distribution has been implemented in our program and has been extended for higher dimensions. 

Note that for 2D, $l_x>0$, but $l_z$ can take both positive and negative values. However, for isotropic turbulence, the structure functions for $+l_z$ and $-l_z$ are statistically equal. Therefore, in our computations, we restrict  to $l_x>0$, $l_z>0$. For anisotropic turbulence, not discussed here, the structure functions will depend on $(l_x,l_z)$ rather than $l$; our code will be extended to such systems in future. 

For 3D turbulence, the structure functions will depend on $(l_x,l_y,l_z)$. We divide the tasks among processors over $l_x$ and $l_y$ as done  above for 2D turbulence. The aforementioned algorithm can be easily extended to the 3D case. We employ a similar method for the computation of scalar structure functions as well.
 
In the next section, we discuss the scaling of our code.

# Scaling of `fastSF`

`fastSF` is scalable over many processors due to vectorization and equal load distribution. We demonstrate the scaling of `fastSF` for the third-order longitudinal structure function for an idealized velocity field on a $128^3$ grid.  For our computation we employ a maximum of 1024 processors. We take the velocity field as
$$\mathbf{u} = 
\begin{bmatrix} 
x \\ y \\z
\end{bmatrix}.$$
We perform four runs on a Cray XC40 system (Shaheen II of KAUST) for this problem using 16, 64, 256, and 1024 processors. In Fig. \ref{Scaling}, we plot the inverse of time taken in seconds versus the number of processors. The best fit curve for these data points yields
$$T^{-1} \sim p^{0.986 \pm 0.002},$$
Thus, the data-points follow $T^{-1} \sim p$ curve to a good approximation. Thus, we conclude that our code exhibits strong scaling. 

![Scaling of `fastSF` for the computation of third-order longitudinal velocity structure function using 16, 64, 256, and 1024 processors of Shaheen II. All the runs were conducted on a $128^3$ grid.  We observe a linear scaling. \label{Scaling}](SF_scaling.png)

 
# Validation

We validate `fastSF` by comparing the numerical results with analytical results for idealized $\mathbf{u}$ and $\theta$ fields as well as with the predictions of K41 [@Kolmogorov:Dissipation; @Kolmogorov:Structure].

### Problem 1

We consider the following 2D velocity and scalar fields:
$$\mathbf{u} = 
\begin{bmatrix} 
x \\ z
\end{bmatrix}, \quad \theta = x+z.
$$
For the above fields, it can be analytically shown that the longitudinal and the transverse velocity structure functions and the scalar structure functions are
$$S_q^{u_\parallel} = (l_x^2 + l_z^2)^{q/2} = l^q,$$
$$S_q^{u_\perp} = 0,$$
$$S_q^\theta = (l_x+l_z)^q.$$
We run ``fastSF`` to compute the velocity and scalar structure functions for the above fields. The resolution of the fields and the domain size are $32^2$ and $1 \times 1$ respectively. We plot the second and the third-order longitudinal velocity structure functions versus $l$ in Fig. \ref{SFTest}. Clearly, $S_2^{u_\parallel}(l)$ and $S_3^{u_\parallel}(l)$ equal $l^2$ and $l^3$ respectively, consistent with the analytical results. Figure \ref{SFScalar} exhibits the density plots of the computed second-order scalar structure function $S_2^{\theta}(\mathbf{l})$ along with $(l_x + l_z)^2$. The two plots are very similar, thus showing the ``fastSF`` computes the scalar structure function correctly.

![For the velocity field defined in Problem 1: plots of the second and third-order longitudinal structure functions vs. $l$. The second and third-order structure functions equal $l^2$ and $l^3$ respectively.\label{SFTest}](SF_test.png)

The above problem is used as a test case for the the code. The user is required to execute the shell script `fastSF/runTest.sh` to run the test case. On doing so, the code generates the velocity and the scalar fields as per the above relation. After computing the structure functions, the code computes the percentage difference between the theoretical and the computed values of the structure functions. If the error does not exceed $1\times 10^{-10}$, the code is deemed to have passed.

![For the scalar field defined in Problem 1: (a) Density plot of the second-order scalar structure function as function of the displacement vector. (b) Density plot of $(l_x+l_z)^2$, which is the analytical value of the second-order scalar structure function. The two density plots are very similar.\label{SFScalar}](SF_scalar.png)


### Problem 2

Here, we consider the classical results of Kolmogorov [@Kolmogorov:Dissipation; @Kolmogorov:Structure] for 3D incompressible hydrodynamic turbulence with homegeneity and isotropy. In such flows, for the inertial range, which comprises of scales lying between the large-scale forcing regime and the small-scale dissipation regime, the third-order longitudinal velocity structure function is given by
$$S_3^{u_\parallel}(l) = -\frac{4}{5} \epsilon l,$$
where $\epsilon$ is the viscous dissipation rate [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book; @Verma:book:ET]. For an arbitrary order $q$, @She:PRL1994 proposed that the longitudinal structure functions scale as $S_3^{u_\parallel}(l) \sim l^{\zeta_q}$, where the exponent $\zeta_q$ is given by 
$$ \zeta_q = \frac{q}{9} + 2 \left ( 1 - \left ( \frac{2}{3} \right )^{q/3} \right ).$$ 


We compute the longitudinal velocity structure functions of $q=3,5,7$ using the simulation data of 3D hydrodynamic turbulence with Reynolds number (Re) of 5700. The simulation was performed using TARANG [@Verma:Pramana2013tarang; @Chatterjee:JPDC2018] on a $512^3$ grid with the domain size of ($2\pi \times 2\pi \times 2\pi$). For more details on the simulation, refer to @Sadhukhan:PRF2019. We run ``fastSF`` on Shaheen II to compute the structure functions, employing 4096 MPI processes. 

![For 3D homogeneous isotropic turbulence (Problem 2): plots of the negative of normalized third, fifth and seventh-order structure functions vs. $l$. The negative of the normalized third-order structure function is close to $4/5$ (dashed line) in the inertial range. \label{Hydro}](SF_hydro.png)

We normalize the third, fifth, and seventh-order longitudinal velocity structure functions with $(\epsilon l)^{\zeta_q}$, where $\zeta_q$ is given by She-Leveque's relation. We plot the negative of these quantities versus $l$ in Fig. \ref{Hydro}. 
The figure clearly shows that in the inertial range ($0.2<l<0.8$), the normalized third-order longitudinal velocity structure function is fairly close to $4/5$ (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of $l$, thus exhibiting consistency with She-Leveque's model. Note that we expect more accurate results for higher resolution simulations [@Verma:Pramana2013tarang].

The results obtained from Problems 1 and 2 thus validate ``fastSF``. 

# Conclusions

This paper describes the design and the validations of ``fastSF``, an efficient parallel C++ code that computes structure functions for given velocity and scalar fields. We validate ``fastSF`` using two test cases. In the first case, we compute the structure functions using hypothetical velocity and scalar fields, and find them to be consistent with analytical results. In the second case, we compute the velocity structure functions using the fields obtained from the simulations of 3D hydrodynamic turbulence, and show consistency with Kolmogorov's K41 result. We believe that ``fastSF`` will be useful to turbulence community.  


# Acknowledgements

We thank R. Samuel, A. Chatterjee, S. Chatterjee, and M. Sharma for useful discussions during the development of ``fastSF``. Our computations were performed on Shaheen II at KAUST supercomputing laboratory, Saudi Arabia, under the project k1416. 

---

# References


