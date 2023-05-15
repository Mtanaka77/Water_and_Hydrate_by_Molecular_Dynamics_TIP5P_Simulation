## Molecular_Dynamics_of_Water_by_TIP5P_Code ##

Procedures of water molecules dynamics simulation are shown for the 5-points molecules. 
This approach is done with five-water molecules with two hydrogens and two L1, L2 hydrogens of dummy sites. 
The fifth site of an oxygen site is used with Lennard-Jones potential Psi(r)= eps_A/r^12 -eps_B /r^6 (Ref. 1). 

The ice state of freezing due to microwaves, our theory discovery in J.Chem.Phys. 2007 (Ref. 2), remains the same, due to
the structure of six-membered ice ! 

### Procedure of Water Molecules by 5-Points Method ###

A. Five sites are oxygen (O), hydrogen 1 and 2 (H), and hydrogen virtual L sites (L). 
Their charges are 0, 0.241e, 0.241e, -0.241e and -0.241e, respectively. The L1 and L2 are called the dummy sites which have null masses.

B. Separate {\bf R}_{i}, {\bf V}_{i} and {\bf r}_{i,k} for water with i=1-N molecules, and 
{\bf s}_{i,k}= (x_{i,k},y_{i,k},z_{i,k}) means for the five sites k=1-5. The separation is done at the starting step only, once determined at t=0, they are constant in time.

C. The half time step is first executed for a predictor step, and the full time step is made for
advance of time. 

D. Before the end of one step, the forces are calculated at {\bf r}_{i,k}= {\bf R}_{i} +A^(-1){\bf s}_{k}
with the three sites of k=1-3, and the L sites are also calculated by algebraic operation.

E. After correction of quatenions, go to the beginning of the clcle.  The leap-frog method is used for the plasmas and waters.


### Each Step of Molecular Dynamics ###

Each step of molecular dynamics is tranlation (1), rotation (2-4), and adding to the fields (5-8).

0. Read quaternions e0-e3 from the file by 'read(30) e0,e1,e2,e3'

1. Sum up five sites and advance, 'd{\bf V}_{i}/dt=\sum_{k=1,5} {\bf F}_{i,k}/m_{i}, 
d{\bf R}_{i}/dt={\bf V}_{i}' for each of the translation motion.

2. For the rotation motion 'd{\bf L}_{i}/dt=\sum_{k} (y_{i,k}F_{i,k}^z-z_{i,k}F_{i,k}^y,
z_{i,k}F_{i,k}^x-x_{i,k}F_{i,k}^z,x_{i,k}F_{i,k}^y-y_{i,k}F_{i,k}^x)'
The summations are made over the five sites of k=1,5. 

3. 'omega_{i,alpha}=(A_{alpha,1)L_{x}+A_{alpha,2)L_{y}+A_{alpha,3)L_{z})/Im_{i,alpha}', 
the angular frequency for speceis A_{alpha,beta} and inertia moments Im_{i,alpha} 
for the directions alpha=x,y,z.

4. 'd{\bf q}_{i}/dt =(1/2)Q(e0,e1,e2,e3)(omega_{i,x),omega_{i,y),omega_{i,z),0), 
d{\bf q}_{i}/dt of Q and omega's has four components found in Goldstein's book.

5. Get a new rotation matrix A_{ij}(e0,e1,e2,e3) in p.205 of Goldstein's book.

6. x_{i}= X_{i} +(A_{11}x_{i}+A_{21}y_{i}+A_{31}z_{i}, 
   y_{i}= Y_{i} +(A_{11}x_{i}+A_{21}y_{i}+A_{31}z_{i},
   z_{i}= Z_{i} +(A_{11}x_{i}+A_{21}y_{i}+A_{31}z_{i},
at the three sites {\bf r}_{i,k} and the position {\bf R}_{i}. The dummy sites are
determined by algebraic vector operation.

7. Forces by Coulombic and Lennard-Jonnes potentials are calculated using five sites.

8. Correction and normalization by quaternions are made at every 10-step interval. Then, go to the next time step as (1).

Note that the choice of a time step is important. For TIP5P case, it may be dt=0.025, else the code is inaccurate or goes overflow.


### The Lennard-Jones Potential ###

With the Coulombic interactions, the Lennard-Jones 12-6 potential is adopted,
A=3.85x10^(-8) erg Ang^(12), B=4.36x10^(-11) erg Ang^(6) for TIP5P-Ewald sum's case, 
and A=4.17x10^(-8) erg Ang^(12), B=4.24x10^(-11) erg Ang^(6) for TIP4P.

Other parameters are: r(OH)=0.9572 Ang, Delta(HOH)=104.52 Ang. r(OM)=0.15 Ang is used 
for the TIP4P case, where the equipartition line of the virtual M sites is on the plane 
that equally separates the HOH angle. 

The A_{ij} and e0,e1,e2,e3(i) equations are written in PDF, "Water_TIP5P_Simulation.pdf".

### To Start a Run ###

To start a simulation of water cluster with the TIP5P code, the adjascent 4x4 hydrogen pairs are summed  
electrostatically, and oxygen pairs are coupled by TIP5P Lennard-Jones potentials.

1. To get an initial state, we make the size of at least a 6x6x6 water cluster for numerical stability.
Short-range and long-range Coulombic forces are best separated for interactions, anf the short-range 
forces are made to be spatially dumped.

2. Around a given temperature, a dryrun is executed at least for 5 periods that is 10^(-9) seconds.
The long dryrun is very important !!

3. Then, we apply the electric field E_x= E_0*sin(omega*time) in the x-direction to excite the 
electric dipole interactions of water. For the moment, we give the electric field 10 GHz where 
the electric field E_0 and electric dipole p_0 are of the order of 5x10^(-3) eV.

### To Obtain the Initial Equilibrium for 298 K ###

We will use salt ions of Na(+) and Cl(-) initially as the dryrun to give random noises. 
The 6x6x6 water clusters have 64 Na(+) and 64 Cl(-) ions, and a run time is t=1700.
We can see random water clusters. 
Afterwards, the salt ions are gradually removed, and the dryrun is continued for 5 periods.

The dryrun is shown with color pictures in PDF, "Water_TIP5P_Simulation.pdf".
The code including salt ions could be slower than the water code even if the salt ions are not present.

### References ### 

1. Classical Mechanics, H. Goldstein, C. Poolee, J. Safko, 3rd Edition, Pearson Education Inc., England, 2003; 
古典力学，吉岡書店，2006.
![image](https://github.com/Mtanaka77/Molecular_Dynimics_of_Water_by_TI5P/assets/111667711/aacddf22-0d92-4f3c-ae38-e1115261fd58)

2. M.Tanaka and M.Sato, Microwave heating of water, ice and saline solution: Molecular dynamics study, J.Chem.Phys., 126, 034509 1-9 (2007).
