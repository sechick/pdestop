function [ figout ] = PDEPlotSolnFigFiles( fignum, PDEstructvec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(PDEstructvec)
    if PDEstructvec(i).Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
        if ~exist('fignum','var'), fignum = 20; end;
        fignum = UtilPlotDiagnostics(fignum, PDEstructvec(i));           % generate a bunch of plots which can be used as diagnostics of computations
    end
end

figout = fignum;

end

