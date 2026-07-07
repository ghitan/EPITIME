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
See [[1](#epitime-paper)] for more detailed 
information on software architecture and its use.

Examples of use are also provided in the package, illustrating experiments 
on a few meaningful problems: 
asymptotic behaviour [[2](#jcd-paper)], 
inverse reconstruction of an infectivity kernel from COVID 19 incidence data 
[[3](#mbe-paper)],
and behavioural dynamics under different memory kernels
[[4](#jmb-paper)]. 

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
**Release date:** April-July 2026

### References

1. <a name="epitime-paper"></a>___EPITIME: A Computational Framework for Integral Epidemic Models 
   with Structure-Preserving Discretizations___<br/>{#epitime-paper}
   B. Buonomo, E. Messina, C. Panico, M. Pezzella, G. Zanghirati<br/>
   (2026) [arXiv:2605.00067v1](https://arxiv.org/abs/2605.00067v1),
   submitted.
2. <a name="jcd-paper"></a>___A non-standard numerical scheme for an age-of-infection epidemic model___<br/>
   E. Messina, M. Pezzella, A. Vecchio<br/>
   Journal of Computational Dynamics, 2022, 9(2): 239–252.<br/>
   [DOI: 10.3934/jcd.2021029](https://www.aimsciences.org/article/doi/10.3934/jcd.2021029)
3. <a name="mbe-paper"></a>___Nonlocal finite difference discretization of a class of renewal equation
   models for epidemics___<br/>
   E. Messina, M. Pezzella, A. Vecchio<br/>
   Mathematical Biosciences and Engineering, 2023, 20(7): 11656–11675.<br/>
   [DOI: 10.3934/mbe.2023518](https://www.aimspress.com/article/id/6458c560ba35de3de31bb8be)
4. <a name="jmb-paper"></a>___An integral renewal equation approach to behavioural epidemic
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

* **Matlab**<br/>
  This folder contains the Matlab source code of the EPITIME software,
  as well as test scripts to reproduce paper results.
	* **Age_of_Infection_Model**<br/>
	  The souces in this folder implement and test the age-of-infection epidemic model    
		* `NSFD_AoI.m`<br/>
		  Main routine implementing in Matlab the unconditionally positive, monotonicity-preserving, 
		  non standard finite
		  difference (NSFD) scheme for the solution of the age-of-infection (AoI)
		  epidemic model [[2](#jcd-paper)].
		* `NSFD_AoI_LIVE.mlx`<br/>
		  Matlab live script to directly test the AoI model. Some descriptions are merged
		  with the source code lines, to help the user follow the computations and
		  better understand the results.
		* `EPITIME_SimulationTool.mlapp`<br/>
		  Complete graphics user interface (GUI) in Matlab to interactively work with the
		  AoI epidemic model. It also exports the results to the Matlab main workspace.
		* `AoI_Simulation_Tool.mlx`<br/>
		  Matlab live script counterpart of the GUI to test the AoI model. Descriptions and
		  recalls are merged with the code, to help the user better understand
		  what the code is doing.
		* `Simple_example.m`<br/>
		  Prototype script showing how to set the data, call the main `NSFD_AoI.m` routine 
		  and visualize the experimental results. 
		* `Test1_Paper_JCD.m`<br/>
		  Matlab test script to reproduce the experimental order of convergence of the AoI model
		  (see [[2](#jcd-paper), section 5]).
		* `AoI_Final_Size.m`<br/>
		  Matlab test script to investigate the long-time behaviour of the AoT model.
		* `AoI_Peak_Estimation.m`<br/>
		  Matlab scripot to investigate the epidemic peak forecasted by the given AoI model.
		* `AoI_Performance_Benchmark.m`<br/>
		  Test script for a computational complexity analysis of the `NSFD_AoI()` Matlab routine.
		* `Calibration_output.mat`<br/>
		  Matlab binary data file for the parameters calibration of the COVID-19 
		  infectivity kernel estimation (see [[1](#epitime-paper), section 5.2])
		* `Italy-daily-trend.csv`<br/>
		  The official 5-years COVID-19 infection data in Italy, as reported by the
		  ''*Il Sole 24 Ore*'', used to infer the infectivity kernel (see [[1](#epitime-paper), section 5.2])
	* **Behavioural_Model**<br/>
      The souces in this folder implement and test the behavioural epidemic model    
		* `NSFD_behavioural.m`<br/>
		  Weighted non standard finite difference (NSFD) scheme for the solution of the behavioural
		  integro-differential epidemic model with demographic turnover and information index
		  [[4](#jmb-paper)].
		* `plotData.m`<br/>
		  Customized plotting routine to show the experimental results of section 6 in  
		  [[4](#jmb-paper)].
		* `Behavioural_Trapezoidal_IF.m`<br/>		
		  Script to test the exponential memory kernel (see [[4](#jmb-paper), section 6.2]).
		* `Behavioural_Unimodal_IF.m`<br/>
		  Script to test the unimodal memory kernel (see [[4](#jmb-paper), section 6.2]).
	* **General_Framework**<br/>
	  Sources in this folder refers to the renewal framework.
		* `NSFD_Renewal.m`<br/>
		  Main routine for the unconditionally positive, non standard finite difference
		  (NSFD) scheme for the solution of a comprehensive renewal framework, which includes
		  different variants of the age-of-infection (AoI) model
		  (see [[3](#mbe-paper), section 2]).
		* `Test1_Paper_MBE.m`<br/>
		  Test script to reproduce the results of AoI model with a low-regularity kernel
		  (see [[3](#mbe-paper), section 5]).
		* `Test2_Paper_MBE.m`<br/>
		  Test script to reproduce the AoI modelling with symptomatic and asymptomatic infection
		  (see [[3](#mbe-paper), section 5]).
		* `Test3_Paper_MBE.m`<br/>
		  Test script to reproduce the results of the AoI simulation with heterogeneous mixing
		  (see [[3](#mbe-paper), section 5]).
* **Python**<br/>
  This folder contains the Pyython sources of the EPITIME software, together
  with some examples and test scripts.
	* Age_of_Infection_Model<br/>
	  Python sources folder for the age-of-infection model.
		* `NSFD_AoI.py`<br/>
		  Main function implementing in Python the unconditionally positive, monotonicity-preserving, 
		  non standard finite
		  difference (NSFD) scheme for the solution of the age of infection (AoI)
		  epidemic model [[2](#jcd-paper)].
		* `AoI_Simulation_Tool.ipynb`<br/>
		  Jupiter notebook to test the AoI model, where descriptions and
		  recalls are merged with the code, to help the user better understand
		  what the code is doing.
		* `Test_NSFD.py`<br/>
		  Prototype script to run the age-of-infection model esperiments with 
		  the NSFD scheme.
		* `Test1_Paper_JCD.py`<br/>
	      Python code to reproduce the experimental order of convergence of the AoI model
	      (see [[2](#jcd-paper), section 5]).
		* `AoI_Final_Size.py`<br/>
		  Python code to investigate the long-time behaviour of the AoT model.
		* `AoI_Peak_Estimation.py`<br/>
		  Python code to investigate the epidemic peak forecasted by the given AoI model.
		* `AoI_Performance_Benchmark.py`<br/>
		  Code for a computational complexity analysis of the `NSFD_AoI()` Python routine.
	* Behavioural_Model<br/>
	  Python sources folder for the age-of-infection behavioural model.	
		* `NSFD_behavioural.py`<br/>
		  Python implementation of the weighted non standard finite difference (NSFD) scheme 
		  for the solution of the behavioural integro-differential epidemic model with 
		  demographic turnover and information index
		  [[4](#jmb-paper)].
		* `plotData.py`<br/>
		  Custom plotting function to show the experimental results of section 6 in
		  [[4](#jmb-paper)].
		* `Behavioural_Trapezoidal_IF.py`<br/>
		  Python code to test the exponential memory kernel (see [[4](#jmb-paper), section 6.2]).
		* `Behavioural_Unimodal_IF.py`<br/>
		  Python code to test the unimodal memory kernel (see [[4](#jmb-paper), section 6.2]).
		* `Behavioural_Peaks_Estimation.py`<br/>
		  Python code to investigate the epidemic peaks forecasted by the given model.
		* `Behavioural_Performance_Benchmark.py`<br/>
		  Code for a computational complexity analysis of the `NSFD_behavioural()` Python routine.

### Acknowledgemnts

We thank the anonymous referees for suggesting to add the following test scripts to the repository
(both .m and .py):
AoI\_Peak\_Estimation, AoI\_Performance\_Benchmark, Behavioural\_Peaks\_Estimation, and
Behavioural\_Performance\_Benchmark.
 
