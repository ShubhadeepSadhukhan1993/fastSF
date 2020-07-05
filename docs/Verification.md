# Testing and Validation

We validate `fastSF` by comparing the numerical results with analytical results for idealized \f$\mathbf{u}\f$ and \f$\theta\f$ fields as well as with the predictions of K41 [@Kolmogorov:Dissipation; @Kolmogorov:Structure].

### Problem 1

We consider the following 2D velocity and scalar fields:
\f$\mathbf{u} = 
\begin{bmatrix} 
x \\ z
\end{bmatrix}, \quad \theta = x+z.
\f$
For the above fields, it can be analytically shown that the longitudinal and the transverse velocity structure functions and the scalar structure functions are
\f$S_q^{u_\parallel} = (l_x^2 + l_z^2)^{q/2} = l^q,\f$
\f$S_q^{u_\perp} = 0,\f$
\f$S_q^\theta = (l_x+l_z)^q.\f$
We run ``fastSF`` to compute the velocity and scalar structure functions for the above fields. The resolution of the fields and the domain size are \f$32^2\f$ and \f$1 \times 1\f$ respectively. We plot the second and the third-order longitudinal velocity structure functions versus \f$l\f$ in Fig. \ref{SFTest}. Clearly, \f$S_2^{u_\parallel}(l)\f$ and \f$S_3^{u_\parallel}(l)\f$ equal \f$l^2\f$ and \f$l^3\f$ respectively, consistent with the analytical results. Figure \ref{SFScalar} exhibits the density plots of the computed second-order scalar structure function \f$S_2^{\theta}(\mathbf{l})\f$ along with \f$(l_x + l_z)^2\f$. The two plots are very similar, thus showing the ``fastSF`` computes the scalar structure function correctly.

![For the velocity field defined in Problem 1: plots of the second and third-order longitudinal structure functions vs. \f$l\f$. The second and third-order structure functions equal \f$l^2\f$ and \f$l^3\f$ respectively.\label{SFTest}](SF_test.png)

The above problem is used as a test case for the the code. The user is required to execute the shell script `fastSF/runTest.sh` to run the test case. On doing so, the code generates the velocity and the scalar fields as per the above relation. After computing the structure functions, the code computes the percentage difference between the theoretical and the computed values of the structure functions. If the error does not exceed \f$1\times 10^{-10}\f$, the code is deemed to have passed.

![For the scalar field defined in Problem 1: (a) Density plot of the second-order scalar structure function as function of the displacement vector. (b) Density plot of \f$(l_x+l_z)^2\f$, which is the analytical value of the second-order scalar structure function. The two density plots are very similar.\label{SFScalar}](SF_scalar.png)


### Problem 2

Here, we consider the classical results of Kolmogorov [@Kolmogorov:Dissipation; @Kolmogorov:Structure] for 3D incompressible hydrodynamic turbulence with homegeneity and isotropy. In such flows, for the inertial range, which comprises of scales lying between the large-scale forcing regime and the small-scale dissipation regime, the third-order longitudinal velocity structure function is given by
\f$S_3^{u_\parallel}(l) = -\frac{4}{5} \epsilon l,\f$
where \f$\epsilon\f$ is the viscous dissipation rate [@Kolmogorov:Dissipation; @Kolmogorov:Structure; @Frisch:book; @Verma:book:ET]. For an arbitrary order \f$q\f$, @She:PRL1994 proposed that the longitudinal structure functions scale as \f$S_3^{u_\parallel}(l) \sim l^{\zeta_q}\f$, where the exponent \f$\zeta_q\f$ is given by 
\f$ \zeta_q = \frac{q}{9} + 2 \left ( 1 - \left ( \frac{2}{3} \right )^{q/3} \right ).\f$ 


We compute the longitudinal velocity structure functions of \f$q=3,5,7\f$ using the simulation data of 3D hydrodynamic turbulence with Reynolds number (Re) of 5700. The simulation was performed using TARANG [@Verma:Pramana2013tarang; @Chatterjee:JPDC2018] on a \f$512^3\f$ grid with the domain size of (\f$2\pi \times 2\pi \times 2\pi\f$). For more details on the simulation, refer to @Sadhukhan:PRF2019. We run ``fastSF`` on Shaheen II to compute the structure functions, employing 4096 MPI processes. 

![For 3D homogeneous isotropic turbulence (Problem 2): plots of the negative of normalized third, fifth and seventh-order structure functions vs. \f$l\f$. The negative of the normalized third-order structure function is close to \f$4/5\f$ (dashed line) in the inertial range. \label{Hydro}](SF_hydro.png)

We normalize the third, fifth, and seventh-order longitudinal velocity structure functions with \f$(\epsilon l)^{\zeta_q}\f$, where \f$\zeta_q\f$ is given by She-Leveque's relation. We plot the negative of these quantities versus \f$l\f$ in Fig. \ref{Hydro}. 
The figure clearly shows that in the inertial range (\f$0.2<l<0.8\f$), the normalized third-order longitudinal velocity structure function is fairly close to \f$4/5\f$ (represented by dashed line), consistent with Kolmogorov's theory. Moreover, the normalized fifth and seventh-order structure functions show a plateau for the same range of \f$l\f$, thus exhibiting consistency with She-Leveque's model. Note that we expect more accurate results for higher resolution simulations [@Verma:Pramana2013tarang].

The results obtained from Problems 1 and 2 thus validate ``fastSF``. 
