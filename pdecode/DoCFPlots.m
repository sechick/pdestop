function [ rval, figout, pdeSolnStruct ] = DoCFPlots( fignum, pdeSoln )

% look in default location for base name of input files for discounted
% reward case, otherwise use the argument passed to this function.
%rval = 0;               % assume failure to load unless loading is completed
[ST,~] = dbstack;
routinename = ST.name;  % get function name

PDELocalInit;
if nargin < 1
    basename = [PDEmatfilebase PDEnodiscbase];   % give location
    [rval, pdeSolnStruct] = PDESolnLoad(basename);
elseif isstruct(pdeSoln)
    pdeSolnStruct = pdeSoln;
    basename = strcat(pdeSoln.Header.PDEparam.matdir,pdeSoln.Header.PDEparam.BaseFileName);
    rval = 1;
else
    basename = pdeSoln;
    [rval, pdeSolnStruct] = PDESolnLoad(basename);
end

if ~rval
    warning(routinename, ': unable to open ', basename);
	return; 
end

myfontsize=16;
mysmallfontsize=14;
points = 144*3; %spacing between labels on contours - made so that only one label appears per line
fracheight = 0.9;   % take up 90% of screen
square = true;      % make plot format to be 'square' if true

%pdeSolnStruct.Header;
PDEscale2=pdeSolnStruct.Header.PDEscale;
PDEparam2=pdeSolnStruct.Header.PDEparam;

% do a bit of input parameter checking
figdir = PDEparam2.figdir;
figsave = true;

% create the directory for the figures if it does not exist already and the figures are
% to be saved
if ~isdir(figdir) && figsave
    mkdir(figdir);
end

%firsts = pdeSolnStruct.Header.firsts;
%lasts = pdeSolnStruct.Header.lasts;
%lasthelds = pdeSolnStruct.Header.lasthelds;
fName = pdeSolnStruct.Header.fName;

% for figure 3
PDEscale2.sigma = 1e5;
PDEscale2.c = 1;
PDEscale2.discrate = 0;
PDEscale2.P = 0;
PDEscaleout = PDEScaleStandardize(PDEscale2);
PDEscale2 = PDEscaleout; % this convuluted code here to verify some pass by reference versus pass by value for matlab structures

%alpha = PDEscale2.alpha;
beta = PDEscale2.beta;
gamma = PDEscale2.gamma;
%approxmeth = PDEparam2.approxmethod;


%% Generate plots similar to those in Chick & Frazier, 2012
% First, figures in spirit of Fig 2 and Fig 3
    for ijk = pdeSolnStruct.Header.StartFileVal:pdeSolnStruct.Header.EndFileVal
        svec = pdeSolnStruct.Data(ijk).svec;
        wvec = pdeSolnStruct.Data(ijk).wvec;

        upvec = pdeSolnStruct.Data(ijk).upvec;
        downvec = pdeSolnStruct.Data(ijk).downvec;
        %up1 = pdeSolnStruct.Data(ijk).up1;
%        down1 = pdeSolnStruct.Data(ijk).down1;

%        Vvec = stopfunc(wvec(:),sval,PDEscale2,pdeSoln.Header.PDEparam); % compute value of stopping immediately
        Bwsmatrix = pdeSolnStruct.Data(ijk).Bwsmatrix;
%        ENwsmatrix = pdeSolnStruct.Data(ijk).ENwsmatrix;
%        EPCSwsmatrix = pdeSolnStruct.Data(ijk).EPCSwsmatrix;
        Vvec = zeros(size(Bwsmatrix)); 
        for j=1:length(svec)
            vveccol  = PDEGetVals(pdeSolnStruct,wvec,svec(j)) ;
            Vvec(:,j) = vveccol;
        end
        
%        dw = wvec(2)-wvec(1);

        % Plots similar to figure 2
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(svec,wvec,Bwsmatrix); 
        clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        plot(svec,upvec,'--','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        plot(svec,downvec,'--','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.1*max(upvec);  tmp(3)=-tmp(4); axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); 
        ylabel('Rescaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('E[value of sequential sampling]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-Like2-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );
        
        % Plots similar to figure 3

        maxbnd = 1.1*max(upvec)/beta;
        minV = 10^(floor(log10(2*maxbnd)))/4;
        V = [ minV/5000 minV/1000 minV/200 minV/50 minV/10 minV/4 minV/2 minV:minV:(2*maxbnd) ];

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(1./svec/gamma,wvec/beta,Bwsmatrix/beta,V); 
        clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        plot(1./svec/gamma,upvec/beta,'--','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        plot(1./svec/gamma,downvec/beta,'--','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.1*max(upvec)/beta;  tmp(3)=-tmp(4); axis(tmp);
        xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); 
        ylabel('Posterior mean, y_t / n_t','FontSize',myfontsize,'FontName','Times')
        title('E[value of sequential sampling]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-Like3-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 5
    % Try to get one-step estimate of boundary
    scale = PDEscale2; %pdeSoln.Header.PDEscale;

    % start to reconstruct the values for plotting....
    sbvec = pdeSoln.Computed.accumsvec;
    up = pdeSoln.Computed.accumupper;
    %lo = pdeSoln.Computed.accumlower;
    %findex = pdeSoln.Computed.fileindx;

    betastepvec=[1 4 10 40 100 400 1000];
    numbetas=length(betastepvec);
    betamatrix=zeros(numbetas,length(up));
    mysmalleps = 1e-10;

    for j=1:numbetas; % find the optimal stopping boundary for the KG_\beta rule, for all \beta in betastepvec
        ybndup1step = up/scale.beta; %initialize vector
        tvec = 1/scale.gamma./sbvec;
        ytinit = ybndup1step(length(ybndup1step))*tvec(length(ybndup1step));
        betamatrix(j,:) = ytinit;
        nrepslookahead=betastepvec(j);
        disp(nrepslookahead);
        for i=length(ybndup1step):-1:1
            [ytinit, ~, ~]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),scale.sigma,tvec(i),nrepslookahead,nrepslookahead*scale.c);
            ybndup1step(i)=ytinit;
            ytinit=ytinit*0.99;
        end
        betamatrix(j,:)=max(mysmalleps,ybndup1step);
    end
    kgstar=max(betamatrix);

    if ~exist('fignum','var'), fignum = 20; end;
    fignum=fignum+1;figure(fignum);
    %loglog(1/scale.gamma./sbvec,up/scale.beta,'--',1/scale.gamma./sbvec,tstcurv/scale.beta,'-.',tvec,ybndup1step./tvec,':o',tvec,ybndupNstep./tvec,'-x')
    %legend('From PDE','Quick Approx.','One-step lookahead','N-step lookahead','Location','NorthEast')
    tstcurv=pdeSoln.Header.PDEparam.approxmethod(sbvec);
    hold off;
    loglog(1/scale.gamma./sbvec,tstcurv/scale.beta,'-.r','LineWidth',1.5);
    hold on
    loglog(1/scale.gamma./sbvec,up/scale.beta,'-k','LineWidth',1.5);
    loglog(1/scale.gamma./sbvec,kgstar./tvec,'--','LineWidth',1.5);
    loglog(tvec', betamatrix ./ (ones(length(betastepvec),1)*tvec),'LineWidth', 1 );
    legendvec = { 'Quick approx','PDE','KG_* (approx)' };
    for j=1:length(betastepvec)
        legendvec = [legendvec, strcat('KG_\beta: \beta=',int2str(betastepvec(j)))];
    end
    
    xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); 
    ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
    set(gca,'FontSize',mysmallfontsize);
    legend(legendvec,'Location','SouthWest');
    
    PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
    tmp=axis; tmp(1) = 0.9*min(1/scale.gamma./sbvec); tmp(2) = 1.1*max(1/scale.gamma./sbvec);
    tmp(3)=1; tmp(4) = 5*max(up)/scale.beta; 
    axis(tmp);
    mytitle = strcat(fName,'Fig-Like5-',int2str(ijk));
    PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    % extras
    pdeSolnStructtmp = pdeSolnStruct;
    pdeSolnStructtmp.Header.PDEscale = PDEscale2;
    fignum = PDEComparePDEboundApproxbound( fignum, pdeSolnStructtmp );

    disp('not yet implemented: C&F 2012 tables 1, 2, 3');

%%%%%%%%%%%%%%%%%%%%%%%%
    figout = fignum;
end