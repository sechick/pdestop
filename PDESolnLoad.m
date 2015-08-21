function [rval, PDESolnStructure] = PDESolnLoad(B1filestring,filestart,fileend)
% PURPOSE: 
%   Load in the data from files in order to have estimates of the
%   standardized optimal expected dynamic reward (OEDR) for the Bandit 
%   papers.  These files can then be used in sequential sampling
%   optimization routines, such as those of the following papers:
%       Chick and Frazier, Management Science
%       Chick and Gans, Management Science, Economics of Sequential
%       Sampling...
% ASSUMPTIONS:
%   Data should have been saved in a sequence of .mat files whose name
%   scheme is the value of the variable B1filestring<n>.mat, where n=0
%   contains the general information for the PDE solution, and n=1,2,...
%   containts additional information related to the solution.
%   Such files can be created with the routine PDESolnCompute.m
% INPUTS: 
%   B1filestring : base name of the files (default: 'PDEData')
%   filestart, fileend: (optional) indices for file to load. While there would usually
%       be a file B1filestring<0>.mat which has fields StartFileVal and 
%       EndFileVal, and those are used to load B1filestring<StartFileVal>
%       to B1filestring<EndFileVal> by default, those values can be overridden by
%       putting in optional values for filestart, fileend. 
% OUTPUTS: 
%   rval: true if loaded ok, false if there was an error in loading the data
%   PDESolnStructure: a structure containing information about the PDE
%       problem data, as well as its solution, including 
%       - header information about problem structure, 
%       - matrices of data for PDE, and of stopping boundaries
%       - computed data which indicates how to reach into the matrices.

% initialize output values and set up housekeeping
rval = 0;               % assume failure to load unless loading is completed
PDESolnStructure = [];  % default to empty solution structure
routinename = 'PDESolnLoad';

% error checking on file name
if (nargin<1) || (length(B1filestring)==0)
    B1filestring = 'PDEData';
    warning(sprintf('%s: using default PDE solution file name %s\n',routinename,B1filestring));
end
    
% attempt to read file with global info about PDE solution
ijk = 0;    
mymat = strcat(strcat(B1filestring,int2str(ijk)),'.mat');
[fid, msg] = fopen(mymat,'r');       % check if file exists
if fid == -1
    error('%s: unable to read PDE solution header %s\n',routinename,mymat);
    return;
else
    fclose(fid);            % file exists, close it up 
	S = load(mymat);        % reopen as a matlab file, read in info

    % check range of files to load: using subset of range if requested by
    % user
    if nargin >= 2
        S.StartFileVal = max(S.StartFileVal, filestart);
    end
    if nargin >= 3
        S.EndFileVal = min(S.EndFileVal, fileend);
    end
    PDESolnStructure.Header = S;    % save this info as a header, use the info in the header to guide the loading of the remaining stuff
    
    % now, try to read in the rest of the file
    % Assumes header has fields S.StartFileVal and S.EndFileVal
    for ijk = S.StartFileVal:S.EndFileVal
        mymat = strcat(strcat(B1filestring,int2str(ijk)),'.mat');
        [fid, msg] = fopen(mymat,'r');       % check if file exists
        if fid == -1
            error('%s: unable to read PDE data file %s\n',routinename,mymat);
            return;
        else
            fclose(fid);            % file exists, close it up 
            PDESolnStructure.Data(ijk) = load(mymat);        % reopen as a matlab file, read in info
        end
    end
    
    % Finally, try to do some computations to build up the file boundary
    % and create indices into the solution data structures above
    startedloading = false;
    S = PDESolnStructure.Header;
    Data = PDESolnStructure.Data;
    for ijk = S.StartFileVal:S.EndFileVal
    	if (startedloading) 
            lastsval = accumsvec(length(accumsvec));
            [a b] = min(lastsval > Data(ijk).svec);              % find indx where new file contains strictly larger s values
            accumsvec=[accumsvec Data(ijk).svec(b:length(Data(ijk).svec))];            % append the s values and boundaries for those larger s values
%            accumwvec=[accumwvec Data(ijk).wvec(b:length(Data(ijk).svec))];            % append the s values and boundaries for those larger s values
            accumupper=[accumupper Data(ijk).upvec(b:length(Data(ijk).svec))];
            accumlower=[accumlower Data(ijk).downvec(b:length(Data(ijk).svec))];
            fileindx=[fileindx ijk*ones(size(Data(ijk).svec(b:length(Data(ijk).svec))))];
        else
            accumsvec=Data(ijk).svec;
%            accumwvec=Data(ijk).wvec;
            accumupper=Data(ijk).upvec;
            accumlower=Data(ijk).downvec;
            fileindx=ijk*ones(size(Data(ijk).svec));
            startedloading=true;
        end
        
    end
    Computed.accumsvec = accumsvec;
    Computed.fileindx = fileindx;       % index into which file structure contains the rigth element...
%    Computed.accumwvec = accumwvec;
    Computed.slower = min(accumsvec);
    Computed.supper = max(accumsvec);
    Computed.accumupper = accumupper;
    Computed.accumlower = accumlower;
    PDESolnStructure.Computed = Computed;

    % if we get to here, then the files have been correctly loaded
    rval = 1;
    
end
