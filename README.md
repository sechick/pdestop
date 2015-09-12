# PDE solver

This site has Matlab code for approximating the solution to certain optimal stopping problems. Examples of such problems are found in the following papers:

Chick SE and Gans N, 2009, Economic Analysis of Simulation Selection Problems, Management Science, 55(3): 421--437, and Electronic Companion.

Chick SE and Frazier P, 2012, Sequential Sampling with Economics of Selection Procedures, Management Science, 58(3), 550--569.

Both of these papers approximate the solution to optimal sequential sampling problems when there are sampling costs to infer the unknown mean of an alternative through statistical sampling, either when there are discounting penalties (C&G) or no discounting penalties (C&F) associated with delays in selecting an alternative.
Both of these papers assume that there is offline learning, and that a decision is made at the end of the sampling to select the alternative, to obtain a one time reward equal to the sampling mean, or to get a 0 reward.
As such, these papers provide results for certain optimal stopping / stoppable bandit problems.

The solutions are approximated using a reverse time diffusion approximation, in the spirit of Chernoff and coauthors, and adapts techniques of Chernoff and Petkau and others, to compute the solution to the resulting continuous time optimal stopping problem with a finite difference technique with grid sizes which vary.
This code includes a couple of improvements in accuracy which have been found to be useful in the computations.

# Workflow to get files for C&G and C&F papers

The code can be used to produce all graphs with PDE computations in the above papers. It does not (yet) produce the numerical tables based on Monte Carlo simulations.

1. Make a local copy of the repo.
2. Launch Matlab and change your directory to the root directory of this code, where PDELocalSetPaths.m is stored.
3. Run PDELocalSetPaths.m to perform some initializations for the path variables.
4. Edit SolvePlotCFCG.m, and cut/paste chunks of code into the Matlab browser.
   - One first generates two sets of files with the solutions for standardized versions of the sequential sampling problems, CF<n>.mat for undiscounted problems, and CG<n>.mat for problems with positive discount rates. Each of these two sets of files can take up approximately 20-30Mb of hard drive space (and are stored in Matfile\). These need be produced only once. They take about 3 min each to generate with a Fujitsu T902 tablet with quad core and 2.6MHz chip.
   - Then, one can load in the solutions with the function PDESolnLoad, in order to approximate the solutions for problems with your parameter values, as functions of the stored solutions for the standardized versions of the sequential sampling problems.

SolvePlotCFCG.m gives such examples. The files DoCGPlots() and DoCFPlots() can be used to produce graphs like those in the papers listed above.



For future work:
1. Online learning for both discounted and undiscounted rewards.
2. Delayed samples (as in the repo github:sechick\htadelay) solved with the Chernoff scaling rather than the `normal' scaling in the number of samples.
3. Monte carlo sampling tools.

The code is in development stages and is provided as is for academic use. 

(c) 2015, S Chick, All rights reserved.

