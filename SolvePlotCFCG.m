% SolvePlotCFCG.m:
% A macro with the code required to:
%   1. Set up the path structure in matlab which allows the general PDE
%   solution code to be called.
%   2. Create files with the solutions to the free boundary PDE problems in
%   Chick & Frazier and in Chick & Gans, for offline learning, and an
%   extension which accounts for online learning (requires c=1 for C&F and
%   requires additional conditions for C&G, to be able to use the online
%   learning solution).
%   3. Loads the files with the solutions into memory, and generates a
%   number of diagnostic plots for those solutions.
%   4. Generates the plots for the C&F paper and for the C&G paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     SETUP:To be called every time.
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDELocalSetPaths(); % set up paths to PDE package
PDELocalInit;       % set up some variables (all starting with PDE - avoid variables with such names

onflag = true;
onflag = false;
THoriz = -1;    % set to below 0 for infinite horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     CREATE FILES WITH PDE SOLUTIONS: Only needs to be called once to create discounted and undiscounted reward base files.
%%     Once the files are created, they need to be recreated - they can be loaded with next step
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the files for C&F and C&G with base case being offline learning,
% call with onflag = true if you would also like to create the files which
% have the online learning equivalents (not fully tested).
% WARNING: The Online Learning for C&F requires a FINITE horizon and
% requires sampling costs per sample to be 1.
% NOTE: This only needs to be done ONCE!!!! (Can redo if you wish to change
% the parameters of the computations so that it is higher resolution or
% lower resolution)
[rval] = PDECreateSolnFiles(PDEmatfilebase, onflag, THoriz);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     LOAD FILES WITH PDE SOLUTIONS: Each time the CF and/or CG files are required, call the following
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following needs to be done every time the files are required.
% Load the files just created into memory. If only one specific file is
% required, it can be loaded individually. See the code in the function for
% how to do so, or view the 
[cgSoln, cfSoln, cgOn, cfOn] = PDELoadSolnFiles(PDEmatfilebase, onflag);
%[rval, cfSoln] = PDESolnLoad(['MatFiles\CF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a lot of plots for those files: display the solutions to the
% standardized PDE for these combinations
if ~exist('fignum','var'), fignum = 20; end;
fignum = PDEPlotSolnFigFiles(fignum, [cgSoln, cfSoln, cgOn, cfOn] ) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the plots from Chick & Frazier (2012)
% Tables are not yet supported, only the PDE solutions, not the Monte
% Carlo results at present.
[ rval, fignum, ~ ] = DoCFPlots( fignum, cfSoln );      % generate plots from Chick & Frazier (2012)
%[ rval, fignum, cfSoln ] = DoCFPlots( fignum, [PDEmatfilebase PDEnodiscbase] );
%[ rval, fignum, cfSoln ] = DoCFPlots( fignum ); % uses default location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the plots from Chick & Gans (2009)
% Tables are not yet supported, only the PDE solutions, not the Monte
% Carlo results at present.
[ rval, fignum, ~ ] = DoCGPlots( fignum, cgSoln );      % generate plots from Chick & Frazier (2012)
