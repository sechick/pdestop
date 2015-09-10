function [ figout ] = PDEComparePDEboundApproxbound( fignum, pdeSoln )
%PDEComparePDEboundApproxbound: Summary of this function goes here
%   starting at figure fignum, plot the pde solution's approximation to the
%   upper boundary, found in pdeSoln structure, as computed by the finite
%   difference method, as well as with the functional approximation to the
%   boundary, as passed in the handle to pdeSoln.... (normally this would
%   be the C&G or C&F approximation to the upper boundary...)

myfontsize=16;
mysmallfontsize=14;
points = 144*3; %spacing between labels on contours - made so that only one label appears per line
fracheight = 0.9;   % take up 90% of screen
square = true;      % make plot format to be 'square' if true

scale = pdeSoln.Header.PDEscale;
fName = pdeSoln.Header.fName;
figdir = pdeSoln.Header.PDEparam.figdir;

% start to reconstruct the values for plotting....
sbvec = pdeSoln.Computed.accumsvec;
up = pdeSoln.Computed.accumupper;
lo = pdeSoln.Computed.accumlower;
findex = pdeSoln.Computed.fileindx;

    if isa(pdeSoln.Header.PDEparam.approxmethod,  'function_handle')
    %    dt = 1/scale.gamma/s0 - 1/scale.gamma/(s0+ds)
    %    dmu=dw/beta
        mysmallfontsize = 14;
        myfontsize = 16;
        fignum=fignum+1;figure(fignum);
        tstcurv=pdeSoln.Header.PDEparam.approxmethod(sbvec);
        plot(sbvec,up,'--',sbvec,tstcurv,'-.');set(gca,'FontSize',mysmallfontsize);
        legend('From PDE','Quick Approx.','Location','SouthEast')
        xlabel('Reverse time scale,  s=1/(\gamma t)','FontSize',myfontsize,'FontName','Times'); ylabel('Rescaled mean, w','FontSize',myfontsize,'FontName','Times')
        %title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')

        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        mytitle = strcat(fName,'Fig-BoundaryWS-');
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

        fignum=fignum+1;figure(fignum);
        loglog(1/scale.gamma./sbvec,up/scale.beta,'--',1/scale.gamma./sbvec,tstcurv/scale.beta,'-.');set(gca,'FontSize',mysmallfontsize);
        legend('From PDE','Quick Approx.','Location','NorthEast')
        xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        %title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        mytitle = strcat(fName,'Fig-BoundaryYT-');
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );
    end
    
figout = fignum;

end

