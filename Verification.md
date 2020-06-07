# Testing and Validation

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