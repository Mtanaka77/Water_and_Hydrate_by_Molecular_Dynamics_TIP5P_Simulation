## Molecular Dynamics Simulation of Water and Ice by TIP5P Code ##

Microwave heating of water and ice is studied by molecular dynamics simulations of the 5-point TIP5P code. 
This is the second generation code to the first one of 3-point SPC/E shake/rattle algorithm (JCP, 2007). 
There are five-water molecules to calculate positive two hydrogens H1, H2 and negative two hydrogens L1, L2 
called dummy sites. The fifth oxygen site is to correlate with adjacent molecules with the Lennard-Jones potential Psi(r)= eps_A/r^12 -eps_B /r^6 (Ref. 1). 

Simulations of molecular dynamics of water and ice are shown at the end of README.md and/or 
by PDF documents on this repository. These are "Water_TIP5P_Simulation.pdf" for numerical coding 
and "Molecular_Dynamics_Simulation_of_Water_and_Ice_by_TIP5P_code.pdf" for physics simulations.   

The ice state of freezing by microwaves, which is our theory discovery in J.Chem.Phys. 2007 (Ref. 2), 
remains basically the same due to the structure of six-membered ice ! 
Methane hydrate simulations are to be executed by the TIP5P code, like with the SPC/E code (Ref. 3).


### Procedure of Water Molecules by 5-Points Method ###

A. Five sites are oxygen (O), hydrogen 1 and 2 (H), and hydrogen virtual L sites (L). 
Their charges are 0, 0.241e, 0.241e, -0.241e and -0.241e, respectively. 
The L1 and L2 are called the dummy sites which have null masses.

B. Separate {\bf R}_{j}, {\bf V}_{j} and {\bf r}_{i} for water with i=1-N, j=1-N/5 molecules, and 
{\bf s}_{i}= (x_{i},y_{i},z_{i}) means for the three sites. The separation is done at the 
starting step only, once determined at t=0, they are constant in time.

C. The half time step is first executed for a predictor step, and the full time step is made for
the three sites of k=1-3, and the L sites are also calculated by algebraic operation.

D. Before the end of one step, the forces are calculated at {\bf r}_{i}= {\bf R}_{j} +A^(-1){\bf s}
with the three sites. The L sites are also calculated by algebraic operation.

E. After correction of quaternions, go to the beginning of the cycle.  The leap-frog method is used for the plasmas and waters.


### Each Step of Molecular Dynamics Code ###

Each step of molecular dynamics simulation consistes of tranlation (step 1), rotation (steps 2-4), 
and adding to the fields (steps 5-8). The step 0 is made only initially. 

0. Read the positions (x,y,z)_{i=1,N}, and quaternions (e0,e1,e2,e3)_{j=1,N/5} from the file by 'read(30) e0,e1,e2,e3'.

1. Summation of five sites of water and to make advance in time, 
'd{\bf V}_{j}/dt=\sum_{k=1,5} {\bf F}_{k}/m_{j}, 
d{\bf R}_{j}/dt={\bf V}_{j}' for each of the translation motion.

2. For the rotation motion 'd{\bf L}_{j}/dt=\sum_{i} (yr_{i}F_{i}^z-zr_{i}F_{i}^y, 
zr_{i}F_{i}^x-xr_{i}F_{i}^z, xr_{i}F_{i}^y-yr_{i}F_{i}^x)', where 
xr_{i}=x_{i}-XC, yr_{i}=y_{i}-YC, zr_{i}=z_{i}-ZC and (XC,YC,ZC) are the gravity center
of the j-th molecule. F_x, F_y, F_z stand for the x,y,z direction of forces. 
The summation over each molecule is made over the five sites. 

3. 'omega_{j}=(A_{alpha,1)L_{x}+A_{alpha,2)L_{y}+A_{alpha,3)L_{z})/Im_{j,alpha}', 
for A_{alpha,1}, A_{alpha,2}, A_{alpha,3} and inertia moments Im_{j,alpha} with 
the directions alpha=x,y,z.

4. 'd{\bf q}_{j}/dt =(1/2)Q(e0,e1,e2,e3)(omega_{j,x),omega_{j,y),omega_{j,z),0), 
d{\bf q}_{j}/dt of Q and omega's have four components found in Goldstein's book.

5. Get a new rotation matrix A_{alpha,beta}(e0,e1,e2,e3) in p.205 of Goldstein's book.

6. x_{i}= X_{j} +(A_{11}xr_{i}+A_{21}yr_{i}+A_{31}zr_{i}, 
   y_{i}= Y_{j} +(A_{12}xr_{i}+A_{22}yr_{i}+A_{32}zr_{i},
   z_{i}= Z_{j} +(A_{13}xr_{i}+A_{23}yr_{i}+A_{33}zr_{i},
where the three components are (xr,yr,zr)_{i}=(A_{11}*(x_{i}-XC)+A_{12}*(y_{i}-YC)+A_{13}*(z_{i}-ZC), 
A_{21}*(x_{i}-XC)+A_{22}*(y_{i}-YC)+A_{23}*(z_{i}-ZC),
A_{31}*(x_{i}-XC)+A_{32}*(y_{i}-YC)+A_{33}*(z_{i}-ZC)), 
and the position {\bf R}_{j}. 
The dummy sites are calculated by algebraic vector operation.

7. Forces by Coulombic interactions and Lennard-Jones potentials are calculated using five sites.
This is the most time consuming part of the TIP5P code.

8. The correction and normalization by quaternions are made in every 10-step interval. 
Then, go to the next time step as step 1.

Note that the choice of a time step is important. For TIP5P case, it may be dt=0.025, else 
the code is inaccurate or/and goes overflow.

The equations of A_{ij} and e0(i),e1(i),e2(i),e3(i) are written in the PDF file, "Water_TIP5P_Simulation.pdf".
Checks of equations of the TIP5P code are also shown.

### The Lennard-Jones Potential ###

With the Coulombic interactions, the Lennard-Jones 12-6 potential is adopted,
A=3.85x10^(-8) erg Ang^(12), B=4.36x10^(-11) erg Ang^(6) for TIP5P-Ewald sum's case, 
and A=4.17x10^(-8) erg Ang^(12), B=4.24x10^(-11) erg Ang^(6) for TIP4P.

Other parameters are: r(OH)=0.9572 Ang, Delta(HOH)=104.52 Ang. r(OM)=0.15 Ang is used 
for the TIP4P case, where the equipartition line of the virtual M sites is on the plane 
that equally separates the HOH angle. 

### To Start a Run ###

To start a simulation of water cluster with the TIP5P code, the adjacent 4x4 hydrogen pairs 
are summed electrostatically, and oxygen pairs are coupled by TIP5P Lennard-Jones potentials.

1. To get an initial state, we make the size of at least a 6x6x6 water cluster for numerical stability.
Short-range Coulombic and long-range forces are best separated for interactions, where the short-range 
forces are made to be spatially dumped.

2. Around a given temperature, a dryrun is executed at least for 5 periods that is 10^(-9) seconds.
The long dryrun of a simulation is very important !!

3. Then, we may be ready to apply the electric field E_x= E_0*sin(\omega *time) in the x-direction to excite the 
electric dipole interactions of water. For the moment, we give the electric field 10 GHz where 
the electric field E_0 and electric dipole p_0 are of the order of 5x10^(-3) eV.
The E_0 value can be as small as 5x10^(-4) eV for the 6x6x6 water cluster.

### To Obtain the Initial Equilibrium for 298 K ###

To give random noises, we will use salt ions of Na(+) and Cl(-) initially as the dryrun. 
The 6x6x6 water clusters use 64 Na(+) and 64 Cl(-) ions and a run time is t=1700, which is 
enough for that amount of water. We can see randomized water clusters. 
Afterwards, the salt ions are gradually removed, and the dryrun for water is continued for 5 periods
of t=50,000.
The dryrun is shown with color pictures in the PDF document, "Water_TIP5P_Simulation.pdf".


### Figures of Simulations by TIP5P Code ###

Molecular dynamics simulation of water and ice are shown using the TIP5P code, 
and the current file is "Water_TIP5P_Simulation.pdf" in this repository.

Figures in the color PDF file are water and ice states that include the energy of molecules, 
3-dimensional scatter plots, pair distribtion functions of oxygen and hidrogen atoms,
and cosine distribution functions in the x-direction.

Generally speaking, water is heated by microwaves where its heating efficiency is highest at 0 Celsius.
Below 273 K, however, ice is frozen and not heated by the electric field.

### Necessary Files ###

Files for simulation 

!*   1. @p3mtip5p07a.f03 : simulation code                       

!*   2. param_tip5p_D07a.h : parameter file 1, be constants      

!*   3. TIP507_config.start0 : parameter file 2, kstart=0 or 2   

!*      and/or TIP507_config.start1 : kstart=1 or 3              

Post-processing files                                

!*   * @lplotip507.f03 - the final history of energy             

!*   * @dipol_seqtip507.f03 - dipole Ex field                    

!*   * @iceplotip507.f03  - for 3D plot for x,y,z                

!*   * @wat_radtip507.f03 - radial distribution functions        
                                                         
References of Numerical Technique

!* 1) M.Tanaka, J.Comput.Phys., vol. 79, 206 (1988).   

!* 2) M.Tanaka, J.Comput.Phys., vol.107, 124 (1993).   

!* 3) M.Tanaka, Comput.Phys.Comm., vol.87, 117 (1995). 

!*  4) M.Tanaka, Comput.Phys.Comm., vol.241, 56 (2019). 


### References ### 

1. Classical Mechanics, H. Goldstein, C. Poolee, J. Safko, 3rd Edition, Pearson Education Inc., England (2003); 
古典力学，吉岡書店 (2006).
 
2. M. Tanaka and M. Sato, Microwave heating of water, ice and saline solution: Molecular dynamics study, J.Chem.Phys., 126, 034509 1-9 (2007).

3. M. Tanaka, M. Sato, S. Nakatani, Microwave heating and collapse of methane hydrate by molecular dynamics simulations, 
arXiv:1909.01024, Cornell University, USA (2019).

