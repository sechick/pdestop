function [rval] = UtilPlotDiagnostics( cfSoln )
% UtilPlotDiagnostics: Plot a bunch of plots for the solution of the PDE
% found in cfSoln. Useful for visualizing the data structures, and for
% diagnostics for the data sets in question. plots are saved in the
% directory named in 'figdir', or to 'Figure\ subdirectory if the directory
% is not passed. plots will not be saved if figdir is passed as [].
% Figdir should end with a '\' if it is not empty string.
%

    myfontsize=16;
    mysmallfontsize=14;
    points = 144*3; %spacing between labels on contours - made so that only one label appears per line

    % do a bit of input parameter checking
    figdir = cfSoln.Header.PDEparam.figdir;
    figsave = true;

    % create the directory for the figures if it does not exist already and the figures are
    % to be saved
    if ~isdir(figdir) & figsave
        mkdir(figdir);
    end

    % print data to console with information about file
    cfSoln.Header
    PDEscale=cfSoln.Header.PDEscale
    PDEparam=cfSoln.Header.PDEparam

    firsts = cfSoln.Header.firsts;
    lasts = cfSoln.Header.lasts;
    lasthelds = cfSoln.Header.lasthelds;
    fName = cfSoln.Header.fName;
    
    alpha = cfSoln.Header.PDEscale.alpha;
    beta = cfSoln.Header.PDEscale.beta;
    gamma = cfSoln.Header.PDEscale.gamma;
    approxmeth = cfSoln.Header.PDEparam.approxmethod;
    
    
    % for each computed block of data, plot a bunch of diagnostics
    for ijk = cfSoln.Header.StartFileVal:cfSoln.Header.EndFileVal
        block = ijk
 
        svec = cfSoln.Data(ijk).svec;
        wvec = cfSoln.Data(ijk).wvec;
        Bwsmatrix = cfSoln.Data(ijk).Bwsmatrix;
        ENwsmatrix = cfSoln.Data(ijk).ENwsmatrix;
        EPCSwsmatrix = cfSoln.Data(ijk).EPCSwsmatrix;
        upvec = cfSoln.Data(ijk).upvec;
        downvec = cfSoln.Data(ijk).downvec;
        up1 = cfSoln.Data(ijk).up1;
        down1 = cfSoln.Data(ijk).down1;

        dw = wvec(2)-wvec(1);

        maxup1 = max(up1) 
        maxupvec = max(upvec) 
        maxwvec = max(wvec)

        
        % First, plot a contour plot of the benefit (above 0) of continuing, plus a
        % boundary of the upper and lower continuation set
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(svec,wvec,Bwsmatrix); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        plot(svec,upvec,'--',svec,downvec,'--');set(gca,'FontSize',mysmallfontsize);
        if ijk>1
            plot([ lasts(ijk-1) lasts(ijk-1) ] ,[min(wvec) max(wvec)],'-.r');
        end
        plot([ firsts(ijk) firsts(ijk) ] ,[min(wvec) max(wvec)],'-.g');
        plot([ lasthelds(ijk) lasthelds(ijk) ] ,[min(wvec) max(wvec)],'-.k');
        tmp=axis;tmp(4)=1.2*max(max(up1),max(upvec));tmp(3)=1.2*min(min(down1),min(downvec));axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Stdized E[value of continuing over stopping | (w_s, s)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContWS',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(1/gamma./svec,wvec/beta,Bwsmatrix/beta); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        plot(1/gamma./svec,upvec/beta,'--',1/gamma./svec,downvec/beta,'--');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(max(up1)/beta,max(upvec)/beta);tmp(3)=1.2*min(min(down1),min(downvec))/beta;axis(tmp);
        xlabel('Effective number of samples, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        title('E[value of continuing over stopping | (y_t/n_t, n_t)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContYT',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end	

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        plot(svec,up1,'-.',svec,upvec,'--');set(gca,'FontSize',mysmallfontsize);
        hold on
        if isa(approxmeth, 'function_handle')
            theoryvec=approxmeth(svec,PDEscale,PDEparam);   
            plot(svec,theoryvec,'-');
            legend('not adjusted','bias adjusted','theory bound');
        else
            legend('not adjusted','bias adjusted');
        end
        plot(svec,down1,'-.',svec,downvec,'--');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(max(up1),max(upvec));tmp(3)=1.2*min(min(down1),min(downvec));axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Stopping boundaries in (w_s,s) scale','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContCPbias',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end	

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        plot(svec,-downvec,'-.',svec,upvec,'--');set(gca,'FontSize',mysmallfontsize);
        hold on
        if isa(approxmeth, 'function_handle')
            theoryvec=approxmeth(svec,PDEscale,PDEparam);   
            plot(svec,theoryvec,'-');
            legend('- downvec','upvec','theory bound');
        else
            legend('- downvec','upvec');
        end
        tmp=axis;tmp(4)=1.2*max(max(-downvec),max(upvec))+dw;tmp(3)=0;axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Compare upper bound with -lower bound, (w_s, s) coord','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigDiffUpDownWS',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end	

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
%        [C, h]=contour(1/gamma./svec,wvec/beta,ENwsmatrix/gamma); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        plot(wvec(:)/beta,ENwsmatrix(:,end)/gamma,'-');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(min(downvec(end))/beta,0,'x',max(upvec(end))/beta,0,'x');set(gca,'FontSize',mysmallfontsize);
        plot(min(cfSoln.Data(ijk).down1(end))/beta,0,'o',(cfSoln.Data(ijk).up1(end))/beta,0,'o');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(1)=(cfSoln.Data(ijk).downvec(end)-dw)/beta;tmp(2)=(cfSoln.Data(ijk).upvec(end)+dw)/beta;axis(tmp);
%        tmp=axis;tmp(4)=1.2*max(ENCin)/gamma;tmp(3)=-1;axis(tmp);
%        xlabel('Effective number of samples, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('E[num samples | s]','FontSize',myfontsize,'FontName','Times')
        title('E[num samples] in (y_t/n_t,n_t) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigENYT',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end	

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
%        [C, h]=contour(1/gamma./svec,wvec/beta,1-EPCSwsmatrix); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        plot(wvec(:),1-EPCSwsmatrix(:,end),'-',wvec(:),1 - normcdf(abs(wvec(:))/sqrt(cfSoln.Header.lasts(ijk)),0,1),'-.');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(downvec(end),0,'x',upvec(end),0,'x');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=max(dw,max(1-EPCSwsmatrix(:,end)));tmp(3)=0;axis(tmp);
%        xlabel('Effective number of samples, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('1-E[PCS]','FontSize',myfontsize,'FontName','Times')
        title('1-E[PCS given (w, s)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigEPCSWS',int2str(ijk),'.eps');
        if figsave print('-deps',mytitle); end	

        WidthContinInw = upvec(end) - downvec(end)
        
        % plot value function as computed for same value of s on
        % consecutive plots
        if ijk < cfSoln.Header.EndFileVal
            svaltotest=cfSoln.Header.lasts(ijk)
            vala = interp2(cfSoln.Data(ijk).svec,cfSoln.Data(ijk).wvec,cfSoln.Data(ijk).Bwsmatrix,svaltotest,cfSoln.Data(ijk).wvec)/beta; 
            valb = interp2(cfSoln.Data(ijk+1).svec,cfSoln.Data(ijk+1).wvec,cfSoln.Data(ijk+1).Bwsmatrix,svaltotest,cfSoln.Data(ijk).wvec)/beta; 
            if ~exist('fignum','var'), fignum = 20; end;
            fignum=fignum+1;figure(fignum);
            plot(cfSoln.Data(ijk).wvec,vala,'-r',cfSoln.Data(ijk).wvec,valb,'-.k',cfSoln.Data(ijk).up1(end),0,'x',cfSoln.Data(ijk).down1(end),0,'x')
            title(sprintf('Solution in two blocks at s=%f',svaltotest));
            legend('block j','block j+1');
            tmp=axis;tmp(1)=min(cfSoln.Data(ijk).down1(end))-dw;tmp(2)=max(cfSoln.Data(ijk).up1(end))+dw;axis(tmp);
            mytitle = strcat(figdir,fName,'RippleCheckAbs',int2str(ijk),'.eps');
            if figsave print('-deps',mytitle); end	

            if ~exist('fignum','var'), fignum = 20; end;
            fignum=fignum+1;figure(fignum);
            plot(cfSoln.Data(ijk).wvec,(vala-valb)./vala,'-',cfSoln.Data(ijk).up1(end),0,'x',cfSoln.Data(ijk).down1(end),0,'x')
            maxrelerr=max(abs(vala-valb)./vala)
            maxabserr=max(abs(vala-valb))
            title(sprintf('Solution in two blocks at s=%f',svaltotest));
            legend('relative error: block j - block j+1');
            tmp=axis;tmp(1)=min(cfSoln.Data(ijk).down1(end))-dw;tmp(2)=max(cfSoln.Data(ijk).up1(end))+dw;axis(tmp);
            mytitle = strcat(figdir,fName,'RippleCheckRel',int2str(ijk),'.eps');
            if figsave print('-deps',mytitle); end	
        end
    end

    accumsvec = cfSoln.Computed.accumsvec;
    accumlower = cfSoln.Computed.accumlower;
    accumupper = cfSoln.Computed.accumupper;
    
    % summary statistics computed from the file
    cfSoln.Computed
    if ~exist('fignum','var'), fignum = 20; end;
    fignum=fignum+1;figure(fignum);
    hold off
    plot(accumsvec,-accumlower,'-.',accumsvec,accumupper,'--');set(gca,'FontSize',mysmallfontsize);
    hold on
    if isa(approxmeth, 'function_handle')
        theoryvec=approxmeth(accumsvec,PDEscale,PDEparam);   
        plot(accumsvec,theoryvec,'-');
        legend('- downvec','upvec','theory bound');
    else
        legend('- downvec','upvec');
    end
    tmp=axis;tmp(4)=1.2*max(max(10*dw,max(-downvec)),max(upvec));tmp(3)=0;axis(tmp);
    xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
    title('Compare upper bound with -lower bound, (w_s, s) coord','FontSize',myfontsize,'FontName','Times')
    mytitle = strcat(figdir,fName,'FigDiffUpDownWS',int2str(ijk),'.eps');
    if figsave print('-deps',mytitle); end	

end