# EPITIME

## EPidemic Integral models TIMe profile Explorer

This is a computational framework for the simulation of two classes of 
integral epidemic models: an age of infection model and an information
dependent behavioural model. 

The framework combines structure preserving Non-Standard Finite Difference 
(NSFD) discretizations with modular implementations in MATLAB and Python, 
together with routines for parameter handling, input validation, performance
assessment, and graphical interaction. 

The proposed methods preserve key qualitative properties of the continuous 
problems, including positivity, boundedness, invariant regions, and correct
long term behaviour, independently of the time step.  

The software is supported by numerical schemes for both model classes and
their main analytical properties, including first order convergence. 
See [[1](https://github.com/ghitan/EPITIME#references)] for more detailed 
information on software architecture and its use.

Examples of use are also provided in the package, illustrating experiments 
on a few meaningful problems: 
asymptotic behaviour [[2](https://github.com/ghitan/EPITIME#references)], 
inverse reconstruction of an infectivity kernel from COVID 19 incidence data 
[[3](https://github.com/ghitan/EPITIME#references)],
and behavioural dynamics under different memory kernels
[[4](https://github.com/ghitan/EPITIME#references)]. 

Overall, EPITIME provides a reliable and accessible computational environment
for the numerical study of renewal epidemic models.

### Authors 

* Bruno Buonomo, Eleonora Messina, Claudia Panico<br/>
  Dept. of Mathematics and Applications "Renato Caccioppoli"<br/>
  University "Federico II", Naples, Italy
* Mario Pezzella<br/>
  Institute for Applied Mathematics “Mauro Picone”<br/>
  National Research Council, Naples, Italy
* Gaetano Zanghirati<br/>
  Dept. of Mathematics and Computer Science<br/>
  University of Ferrara, Ferrara, Italy

**Version:** 1.0<br/>
**Release date:** April 2026

### References

1. ___EPITIME: A Computational Framework for Integral Epidemic Models 
   with Structure-Preserving Discretizations___<br/>
   B. Buonomo, E. Messina, C. Panico, M. Pezzella, G. Zanghirati<br/>
   (2026) [arXiv:2605.00067v1](https://arxiv.org/abs/2605.00067v1),
   submitted.
2. ___A non-standard numerical scheme for an age-of-infection epidemic model___<br/>
   E. Messina, M. Pezzella, A. Vecchio<br/>
   Journal of Computational Dynamics, 2022, 9(2): 239–252.<br/>
   [DOI: 10.3934/jcd.2021029](https://www.aimsciences.org/article/doi/10.3934/jcd.2021029)
3. ___Nonlocal finite difference discretization of a class of renewal equation
   models for epidemics___<br/>
   E. Messina, M. Pezzella, A. Vecchio<br/>
   Mathematical Biosciences and Engineering, 2023, 20(7): 11656–11675.<br/>
   [DOI: 10.3934/mbe.2023518](https://www.aimspress.com/article/id/6458c560ba35de3de31bb8be)
4. ___An integral renewal equation approach to behavioural epidemic
   models with information index___<br/>
   B. Buonomo, E. Messina, C. Panico, A. Vecchio<br/>
   Journal of Mathematical Biology, 2025, 90(8).<br/>
   [DOI: 10.1007/s00285-024-02172-y](https://link.springer.com/article/10.1007/s00285-024-02172-y)

### Citation

If you use this software, please cite the paper [[1](https://github.com/ghitan/EPITIME#references)]
or, better, please check there if a final version has been published.

### License

EPITIME is distributed under the terms of the GNU GPL v. 3 license 
(see the attached LICENSE.txt file).

### Software content

The current version of the EPITIME package is organized in the 
following directory structure:

* Matlab
	* Behavioural_Model
		* `NSFD_behavioural.m`
		* `plotData.m`
		* `Test1_Paper_JMB.m`
		* `Test2_Paper_JMB.m`
	* Age_of_Infection_Model
		* `NSFD_AoI.m`
		* `NSFD_AoI_LIVE.mlx`
		* `EPITIME_SimulationTool.mlapp`
		* `AoI_Simulation_Tool.mlx`
		* `Simple_example.m`
		* `Test1_Paper_JCD.m`
		* `Test2_Paper_JCD.m`
		* `Calibration_output.mat`
		* `Italy-daily-trend.csv`
	* General_Framework
		* `NSFD_Renewal.m`
		* `Test1_Paper_MBE.m`
		* `Test2_Paper_MBE.m`
		* `Test3_Paper_MBE.m`
* Python
	* Behavioural_Model
		* `NSFD_behavioural.py`
		* `plotData.py`
		* `Test1_Paper_JMB.py`
		* `Test2_Paper_JMB.py`
	* Age_of_Infection_Model
		* `NSFD_AoI.py`
		* `AoI_Simulation_Tool.ipynb`
		* `Test_NSFD.py`
		* `Test1_Paper_JCD.py`
		* `Test2_Paper_JCD.py`
 
