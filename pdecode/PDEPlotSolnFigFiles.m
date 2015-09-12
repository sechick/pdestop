function [ figout ] = PDEPlotSolnFigFiles( fignum, PDEstructvec )
%PDEPlotSolnFigFiles: Plots a set of diagnostic plots for the files / blocks
% computed for each PDE solution in the vector of PDE solutions
% (PDEscructvec). In addition to a 'header' file for the PDE solution,
% there are a number of additional files per solution, one per block in the
% Chernoff & Petkau solution technique. In each successive block, the time
% scale and space scale so as to allow the solution to cover several orders
% of magnitude in the number of samples tested.

    for i=1:length(PDEstructvec)
        if PDEstructvec(i).Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
            if ~exist('fignum','var'), fignum = 20; end;
            fignum = UtilPlotDiagnostics(fignum, PDEstructvec(i));           % generate a bunch of plots which can be used as diagnostics of computations
        end
    end

    figout = fignum;

end

