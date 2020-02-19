---
title: 'Kolmogorov41: A hybrid parallel code for computing velocity and scalar structure functions'

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

Turbulence is a complex phenomenon in fluid dynamics involving nonlinear interactions between multiple scales. Structure function is a popular diagnostics tool to study the statistical properites of turbulent flows [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book]. Some of the earlier works comprising of such analysis are those of @Gotoh:PF2002, @Kaneda:PF2003, and @Ishihara:ARFM2009 for three-dimensional homogeneous isotropic turbulence; @Yeung:PF2005 and @Ray:NJP2008 for passive scalar turbulence; @Biferale:NJP2004 for two-dimensional turbulence; and @Kunnen:PRE2008, @Kaczorowski:JFM2013, and @Bhattacharya:PF2019 for turbulent thermal convection. Structure functions are two-point statistical quantities; thus, an accurate computation of these quantities requires averaging over many points. However, incorporation of a large number of points makes the computations very expensive and challenging. Therefore, we require an efficient parallel code for accurate computation of structure functions. In this paper, we describe the design and validation of the results of ``fastSF``, a hybrid parallel code to compute velocity and scalar structure functions. 

 ``fastSF``, written in C++, employs a combination of distributed (MPI) memory parallelization [@Pacheco:book:PP] to compute the structure functions for a given velocity or scalar field. In this code, we minimize the communication between the processors by sharing the input field data in all the MPI processes. Thus, we save considerable time that would otherwise be spent on communication. The user has a choice on the type (scalar or vector) and the dimensions of the fields to be read by the code, and the range of the orders of the structure functions to be computed. The code writes the computed structure functions to hdf5 files that can be further processed by the user. We remark that our code has been used by  for analysing the longitudinal velocity structure functions of turbulent thermal convection.

In the next section, we will briefly define the velocity and the scalar structure functions in turbulence.

# Velocity and scalar structure functions

Let $\mathbf{u}$ and $\theta$ be velocity and scalar fields respectively. For any two points $\mathbf{r}$ and $\mathbf{r+l}$, we define the velocity differential to be $\delta \mathbf{u} = \mathbf{u(r+l)}-\mathbf{u(r)}$. Further, we denote $\delta u_\parallel=\delta \mathbf{u}\cdot \hat{\mathbf{l}}$ as the component of the velocity differential along the vector $\mathbf{l}$, and $\delta u_\perp= |\delta \mathbf{u} - \delta u_\parallel \hat{\mathbf{l}}|$ as the component of the velocity differential perpendicular to $\mathbf{l}$. Assuming statistical homogeneity, we define the longitudinal velocity structure functions of order $q$ as
$$ S_q^{u_\parallel}(\mathbf{l}) = \langle \delta u_\parallel^q \rangle,$$ 
and the transverse velocity structure functions of order 
$q$ as 
$$ S_q^{u_\perp}(\mathbf{l}) = \langle \delta u_\perp^q \rangle. $$ 

Here, $\langle \cdot \rangle$ denotes ensemble averaging. We can also define the scalar differential to be $\delta \theta = \theta (\mathbf{r+l}) - \theta(\mathbf{r})$, and the scalar structure functions (assuming statistical homogeneity) as 
$$ S_q^\theta(\mathbf{l}) = \langle \delta \theta^q\rangle. $$
If the turbulence is isotropic in addition to being homogeneous, the structure functions become functions of $l$, where $l=|\mathbf{l}|$. The second-order velocity structure function provides an estimate of energy in the eddies of size $l$ or less [@Davidson:book:Turbulence]. 

In the next section, we provide a brief description of the code.

# Design of the Code
``fastSF``  is written in C++. The code uses the Blitz++ library for vectorized array operations. The velocity structure functions ($S_q^{u_\parallel}(\mathbf{l})$ and $S_q^{u_\perp}(\mathbf{l})$) and the scalar structure functions ($S_q^{\theta}(\mathbf{l})$) are computed and stored in a Blitz array. ``fastSF`` contains separate functions for computing the scalar and the velocity structure functions. Also, separate functions are called for two-dimensional or three-dimensional input fields. However, the basic design of all the functions is the same, as described below.

Fig. 1 exhibits a brief schematic diagram for the methodology to compute the structure functions. The code navigates through different points in the domain, with each point corresponding to the position vector $\mathbf{l}$. For every $\mathbf{l}$, two subdomains of equal dimensions are considered. The origin of the first subdomain (pink) coincides with that of the main domain. The origin of the second subdomain (green) lies on the position $\mathbf{l}$. The velocity (or scalar) differential field array is then computed by subtracting the velocity (or scalar) field in all the points in the first subdomain from the corresponding points in the second subdomain. The code computes $\delta U_{\parallel}^q$ and $\delta U_{\perp}^q$ arrays for every velocity differential, and then calculates the mean of these arrays to obtain the velocity structure functions. Similarly, the code computes the mean of the scalar differential array to obtain the scalar structure functions. Note that we employ vectorisation in computing the aforementioned arrays as this makes the code efficient and fast. 

Note that an alternative straighter method is to navigate a pair of points through the entire domain and compute the structure functions for every distance between the two points. However, this requires six nested loops for three-dimensional domains (four nested loops for two dimensional domains) because it is not possible to employ vectorisation using this technique. We remark that the strategy used in our code makes it approximately twenty times faster than what it would have been if six nested `for` loops were used for computing the structure functions.

We employ Message Passing Interface (MPI) to parallelise the code.  The points ($\mathbf{l}$) of the domain are distributed among the processors; thus, each processor computes the structure functions for the points assigned to it. The input fields are accessible to all the processors to minimize the need of communication. The function `compute_index_list` distributes the load equally among the processors. The master processor stores the structure function array; it receives the structure functions computed for different $\mathbf{l}$ from all the processors and assigns them to the appropriate index of the structure function array.  

The current version of ‘fastSF’ computes the structure functions correctly only for homogeneous and isotropic turbulence. To save time and memory, ‘fastSF’ computes and stores the structure functions for only the positive values of $l_x$, $l_y$, and $l_z$. Note that this still gives the correct values for isotropic turbulence because in such case, the structure functions depend only on the magnitude of $\mathbf{l}$ and not its orientation.

```
void SF_scalar_2D(Array<double,2> T)
 {
     if (rank_mpi==0) {
         cout<<"\nComputing S(lx, lz) using 2D scalar field data..\n";
     }
     
    int starPt = 0;
    int endPt = Nx*Nz/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Nz);
    Array<double,2> dT;
    
    for (int ix=starPt; ix<endPt; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int z=index_list(ix, 1, rank_mpi);
       			
        dT.resize(Nx-x,Nz-z);
        int count=(Nx-x)*(Nz-z);
        //Calculating the scalar differential
        dT(Range::all(),Range::all())=T(Range(x,Nx-1),Range(z,Nz-1))
                                      -T(Range(0,Nx-x-1),Range(0,Nz-z-1));
        		
        for (int p=0; p<=q2-q1; p++){
            double St = sum(pow(dT(Range::all(),Range::all()),q1+p))/(count);
            Array<int, 1> X, Z, p_arr;
            Array<double, 1> St_arr;
            
            if (rank_mpi==0) {
                X.resize(P);
                Z.resize(P);
                p_arr.resize(P);
                St_arr.resize(P);
            }
        
            MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&St, 1, MPI_DOUBLE_PRECISION, St_arr.data(), 
                       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
            MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (rank_mpi==0) {
                for (int i=0; i<P; i++) {
                    SF_Grid2D_scalar(X(i), Z(i), p_arr(i)) = St_arr(i);
                } 
            }
        }
    }
    if (rank_mpi==0) {
        SF_Grid2D_scalar(0,0,Range::all())=0;
    }
 }
 ```

In the next section, we discuss the validation of our code.
 
# Results

We validate `Kolmogorov41` by using it to compute the structure functions for specific cases, and then comparing our results with those obtained analytically.
We validate `Kolmogorov41` by comparing the numerical results with analytical results for idealized $\mathbf{u}$ and $\theta$ fields as well as with the predictions of K41 [@Kolmogorov:Dissipation; @Kolmogorov:Structure].

### Problem 1

We consider the following two-dimensional velocity and scalar fields:
$$\mathbf{u} = 
\begin{bmatrix} 
x \\ z
\end{bmatrix}, \quad \theta = x+z.
$$
For the above fields, it can be analytically shown that the longitudinal and the transverse velocity structure functions and the scalar structure functions are
$$S_q^{u_\parallel} = (l_x^2 + l_z^2)^{q/2} = l^q,$$
$$S_q^{u_\perp} = 0,$$
$$S_q^\theta = (l_x+l_z)^q.$$
We run ``Kolmogorov41`` to compute the velocity and scalar structure functions for the above fields. The resolution of the fields and the domain size are $32^2$ and $1 \times 1$ respectively. We plot the second and the third-order longitudinal velocity structure functions versus $l$ in Fig. \ref{SFTest}. Clearly, $S_2^{u_\parallel}(l)$ and $S_3^{u_\parallel}(l)$ equal $l^2$ and $l^3$ respectively, consistent with the analytical results. Figure \ref{SFScalar} exhibits the density plots of the computed second-order scalar structure function $S_2^{\theta}(\mathbf{l})$ along with $(l_x + l_z)^2$. The two plots are very similar, thus showing the ``Kolmogorov41`` computes the scalar structure function correctly.

![For the velocity field defined in Problem 1: plots of the second and third-order longitudinal structure functions vs. $l$. The second and third-order structure functions equal $l^2$ and $l^3$ respectively.\label{SFTest}](SF_test.png)

The above problem is used as a test case for the the code. The user is required to execute the shell script `Kolmogorov41-master/runTest.sh` to run the test case. On doing so, the code generates the velocity and the scalar fields as per the above relation. After computing the structure functions, the code computes the percentage difference between the theoretical and the computed values of the structure functions. If the error does not exceed $1\times 10^{-10}$, the code is deemed to have passed.

![For the scalar field defined in Problem 1: (a) Density plot of the second-order scalar structure function as function of the displacement vector. (b) Density plot of $(l_x+l_z)^2$, which is the analytical value of the second-order scalar structure function. The two density plots match identically.\label{SFScalar}](SF_scalar.png)


### Problem 2

Here, we consider the classical results of Kolmogorov [@Kolmogorov:Dissipation; @Kolmogorov:Structure] for three-dimensional incompressible homogeneous isotropic turbulence. In such flows, for the inertial range, which comprises of scales lying between the large-scale forcing regime and the small-scale dissipation regime, the third-order longitudinal velocity structure function is given by
$$S_3^{u_\parallel}(l) = -\frac{4}{5} \epsilon l,$$
where $\epsilon$ is the viscous dissipation rate [@Kolmogorov:Dissipation; @Kolmogorov:Structure; Frisch:book; Verma:book:ET]. For an arbitrary order $q$, @She:PRL1994 proposed that the longitudinal structure functions scale as $S_3^{u_\parallel}(l) \sim l^{\zeta_q}$, where the exponent $\zeta_q$ is given by 
$$ \zeta_q = \frac{q}{9} + 2 \left ( 1 - \left ( \frac{2}{3} \right )^{q/3} \right ).$$ 


We compute the longitudinal velocity structure functions of $q=3,5,7$ using the simulation data of three-dimensional homogeneous isotropic turbulence with Reynolds number (Re) of 5700. The simulation was performed using TARANG [@Verma:Pramana2013tarang; @Chatterjee:JPDC2018] on a $512^3$ grid with the domain size of ($2\pi \times 2\pi \times 2\pi$). For more details on the simulation, refer to @Sadhukhan:PRF2019. We run ``Kolmogorov41`` on a Cray XC40 system (Shaheen II of KAUST) to compute the structure functions, employing 512 MPI and 32 OpenMP processes. The code took $5 \times 10^4$ seconds to complete the computations under the aforementioned parallelization configuration.

![For 3D homogeneous isotropic turbulence (Problem 2): plots of the negative of normalized third, fifth and seventh-order structure functions vs. $l$. The negative of the normalized third-order structure function is close to $4/5$ (dashed line) in the inertial range. \label{Hydro}](SF_hydro.png)

We normalize the third, fifth, and seventh-order longitudinal velocity structure functions with $(\epsilon l)^{\zeta_q}$, where $\zeta_q$ is given by She-Leveque's relation. We plot the negative of these quantities versus $l$ in Fig. \ref{Hydro}. 
The figure clearly shows that in the inertial range ($0.2<l<0.8$), the normalized third-order longitudinal velocity structure function is fairly close to $4/5$ (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of $l$, thus exhibiting consistency with She-Leveque's model. Note that we expect more accurate results for higher resolution simulations [@Verma:Pramana2013tarang].

The results obtained from Problems 1 and 2 thus validate ``Kolmogorov41``. 

# Conclusions

This paper describes the design and the validations of ``Kolmogorov41``, a hybrid parallel C++ code that computes structure functions for given velocity and scalar fields. We validate ``Kolmogorov41`` using two test cases. In the first case, we compute the structure functions using hypothetical velocity and scalar fields, and find them to be consistent with analytical results. In the second case, we compute the velocity structure functions using the fields obtained from the simulations of three-dimensional homogeneous and isotropic turbulence, and show consistency with Kolmogorov's K41 result. We believe that ``Kolmogorov41`` will be useful to turbulence community.  


# Acknowledgements

We thank R. Samuel, A. Chatterjee, S. Chatterjee, and M. Sharma for useful discussions during the development of ``Kolmogorov41``. Our computations were performed on Shaheen II at KAUST supercomputing laboratory, Saudi Arabia, under the project k1416. 

---

# References

