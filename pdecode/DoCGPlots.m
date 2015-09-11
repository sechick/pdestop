function [ rval, figout, pdeSolnStruct ] = DoCGPlots( fignum, pdeSoln )

% look in default location for base name of input files for discounted
% reward case, otherwise use the argument passed to this function.
%rval = 0;               % assume failure to load unless loading is completed
[ST,~] = dbstack;
routinename = ST.name;  % get function name

if nargin < 1
    basename = 'Matfiles\CG';   % give location
    [rval, pdeSolnStruct] = PDESolnLoad(basename);
elseif isstruct(pdeSoln)
    pdeSolnStruct = pdeSoln;
    PDEparam2=pdeSolnStruct.Header.PDEparam;
    basename = strcat(PDEparam2.matdir,PDEparam2.BaseFileName);
    rval = 1;
else
    basename = pdeSoln;
    [rval, pdeSolnStruct] = PDESolnLoad(basename);
end

%pdeSolnStruct.Header;
%PDEscale2=pdeSolnStruct.Header.PDEscale;
PDEparam2=pdeSolnStruct.Header.PDEparam;

if ~rval
    warning(routinename, ': unable to open ', basename);
	return; 
end

myfontsize=16;
mysmallfontsize=14;
points = 144*3; %spacing between labels on contours - made so that only one label appears per line
fracheight = 0.9;   % take up 90% of screen
square = true;      % make plot format to be 'square' if true

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

%alpha = PDEscale2.alpha;
%beta = PDEscale2.beta;
%gamma = PDEscale2.gamma;
%approxmeth = PDEparam2.approxmethod;


%% Generate plots similar to those in the Electronic companion
% First, figures EC3, EC.4 and EC.5 (adapting to the computations in the local
% install for the PDE solution)
%    for ijk = 1:1
    for ijk = pdeSolnStruct.Header.StartFileVal:pdeSolnStruct.Header.EndFileVal
        svec = pdeSolnStruct.Data(ijk).svec;
        wvec = pdeSolnStruct.Data(ijk).wvec;

        upvec = pdeSolnStruct.Data(ijk).upvec;
        %downvec = pdeSolnStruct.Data(ijk).downvec;
        up1 = pdeSolnStruct.Data(ijk).up1;
        %down1 = pdeSolnStruct.Data(ijk).down1;

%        Vvec = stopfunc(wvec(:),sval,PDEscale2,PDEparam2); % compute value of stopping immediately
        Bwsmatrix = pdeSolnStruct.Data(ijk).Bwsmatrix;
%        ENwsmatrix = pdeSolnStruct.Data(ijk).ENwsmatrix;
%        EPCSwsmatrix = pdeSolnStruct.Data(ijk).EPCSwsmatrix;
        Vvec = zeros(size(Bwsmatrix)); 
        for j=1:length(svec)
            vveccol  = PDEGetVals(pdeSolnStruct,wvec,svec(j)) ;
            Vvec(:,j) = vveccol;
        end
        
        %dw = wvec(2)-wvec(1);

        % try to find some good contour values for the contour plot
        maxbnd = 1.1*max(upvec);
        minV = 10^(floor(log10(2*maxbnd)))/4;
        V = [ minV/5000 minV/1000 minV/200 minV/50 minV/10 minV/4 minV/2 minV:minV:(2*maxbnd) ];
        
        % First, plot a contour plot of the benefit (above 0) of continuing, plus a
        % boundary of the upper and lower continuation set
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(1./svec,wvec,Vvec,V); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        plot(1./svec,upvec,'-','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(max(up1),max(upvec)); tmp(3)=-1.2*tmp(4); axis(tmp);
        xlabel('Scaled time, \tau','FontSize',myfontsize,'FontName','Times'); 
        ylabel('Scaled mean, z_\tau / \tau','FontSize',myfontsize,'FontName','Times')
        title('Value function B_1','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-EC345-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    end 

%% Figure EC.6 A and B (at least for boundary: the B_1(0,s) is not as important here, supplemented with 'improved' analytical approximation for boundary
accumsvec = pdeSolnStruct.Computed.accumsvec;
accumupper = pdeSolnStruct.Computed.accumupper;
%accumlower = pdeSolnStruct.Computed.accumlower;
brezzilaiupper = sqrt(accumsvec) .* bl_git_psi(accumsvec,false);
newcgbound = CGApproxBoundW(accumsvec);

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
hold off
plot(accumsvec,accumupper,'-k',accumsvec,newcgbound,'-.r',accumsvec,brezzilaiupper,'--b');
legend('PDE approx', 'Updated C&G approx', 'Brezzi and Lai approx');
PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
tmp=axis; tmp(4)=1.2*max(accumupper); tmp(3)=0; axis(tmp);
xlabel('Scaled time, \tau','FontSize',myfontsize,'FontName','Times'); 
ylabel('Approximation to boundary, z_\tau / \tau','FontSize',myfontsize,'FontName','Times')
mytitle = strcat(fName,'Fig-LikeEC6a-',int2str(ijk));
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

fignum=fignum+1;figure(fignum);
hold off
loglog(1./accumsvec,accumupper,'-k',1./accumsvec,newcgbound,'-.r',1./accumsvec,brezzilaiupper,'--b');
legend('PDE approx', 'Updated C&G approx', 'Brezzi and Lai approx');
PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
tmp=axis; tmp(4)=1.2*max(accumupper); tmp(3)=min(brezzilaiupper)/2; axis(tmp);
xlabel('Scaled time, \tau','FontSize',myfontsize,'FontName','Times'); 
ylabel('Approximation to boundary, z_\tau / \tau','FontSize',myfontsize,'FontName','Times')
mytitle = strcat(fName,'Fig-LikeEC6b-',int2str(ijk));
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

%%%%%%%% Now Graphs from Main paper

generictermreward=@(wvec,s,p1,p2) PDEsimplereward(wvec);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGApproxValue(wvec,s,p1);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1 };

annualdisc = 0.1; % annual discount rate
eta = 20;   % number of minutes per run

% For figure 1: 
scalevec = {'c', 0, 'sigma', 1e7, 'discrate', eta*annualdisc/365/24/60, 'P', 1 };
CGparamvec = { 't0', 1, 'tEND', 10e5, 'precfactor', 6, 'BaseFileName', 'CG' };
paramvec = [CGparamvec, CGfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param );
disp(scale);
disp(param);

wtest = 10.^(5:7) * scale.beta;
stest = PDEInvBound( upperDisc, wtest );
numrepstest = 1 ./ (scale.gamma * stest);

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
hold off
loglog(1./(scale.gamma*accumsvec),accumupper/scale.beta);
hold on
legend('Implement (c=0), \beta^{-1} b_1(1/\gamma t)');
PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
tmp=axis; tmp(1)=1; tmp(2)=10*max(numrepstest); tmp(4)=2*max(wtest)/scale.beta; tmp(3)=min(accumupper)/10/scale.beta; axis(tmp);
for j=1:length(wtest);
    plot( [1 numrepstest(j)], [1 1]*wtest(j) / scale.beta,'--');
    plot( [numrepstest(j) numrepstest(j)], [1 wtest(j) / scale.beta],'--');    
end
xlabel('Number of replications, t','FontSize',myfontsize,'FontName','Times'); 
ylabel('Output mean, y_t / t','FontSize',myfontsize,'FontName','Times');
        text(1.25* numrepstest(2), 1.25 * wtest(2) / scale.beta, 'Stop to implement system');
        text(1.1*numrepstest(3), 0.5 * wtest(1) / scale.beta, 'Stop to implement system');
mytitle = strcat(fName,'Fig1');
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

% Figure 2 for example 2, similar to Figure EC.2
c=1;
scale.c = c;
wtest = (c/scale.discrate) * scale.beta;
stest = PDEInvBound( upperDisc, wtest );
numrepstest = 1 ./ (scale.gamma * stest);
%[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param );

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
hold off
semilogx(1./(scale.gamma*accumsvec),accumupper/scale.beta,'--');
hold on
semilogx(1./(scale.gamma*accumsvec),accumupper/scale.beta - (c/scale.discrate),'-');
    plot( [1 numrepstest], [0 0],'-.');
    plot( [numrepstest numrepstest], [0 -1e6],'-xk');    
legend('Implement (c=0), \beta^{-1} b_1(1/\gamma t)','Implement (c=1), \beta^{-1} b_1(1/\gamma t)-c/\gamma','Mean 0','Limit of contin set');
PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
tmp=axis; tmp(1)=1; tmp(2)=2e5; tmp(4)=2e6; tmp(3)=-1e6; axis(tmp);
        tmp = axis;
        text(0.95*numrepstest,0.6*c/scale.discrate, 'Stop to implement system');
        text(1.05*numrepstest,tmp(3)/2, 'Stop, take NPV=0');
        text(10, tmp(4)/15, 'Continue sampling');
xlabel('Number of replications, t','FontSize',myfontsize,'FontName','Times'); 
ylabel('Output mean, y_t / t','FontSize',myfontsize,'FontName','Times')
mytitle = strcat(fName,'Fig1');
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

% Figure 3 for example 3
g0 = 250000; % fixed cost of project
u0 = 0.25; % time of project in years
delayfactor = exp(-annualdisc * u0);
%c=1;
%wtest = (c/scale.discrate) * scale.beta;
%stest = PDEInvBound( upperDisc, wtest );
%numrepstest = 1 ./ (scale.gamma * stest);

    for ijk = pdeSolnStruct.Header.StartFileVal:pdeSolnStruct.Header.EndFileVal
        svec = pdeSolnStruct.Data(ijk).svec;
        wvec = pdeSolnStruct.Data(ijk).wvec;

        upvec = pdeSolnStruct.Data(ijk).upvec;
        %downvec = pdeSolnStruct.Data(ijk).downvec;
        %up1 = pdeSolnStruct.Data(ijk).up1;
        %down1 = pdeSolnStruct.Data(ijk).down1;

%        Vvec = stopfunc(wvec(:),sval,PDEscale2,PDEparam2); % compute value of stopping immediately
        Bwsmatrix = pdeSolnStruct.Data(ijk).Bwsmatrix;
        %ENwsmatrix = pdeSolnStruct.Data(ijk).ENwsmatrix;
        %EPCSwsmatrix = pdeSolnStruct.Data(ijk).EPCSwsmatrix;
        Vvec = zeros(size(Bwsmatrix)); 
        stopnowval = wvec' * ones(1,length(svec));  % allocate space and put in default values
        for j=1:length(svec)
            vveccol  = PDEGetVals(pdeSolnStruct,wvec,svec(j)) ;
            Vvec(:,j) = vveccol;
            stopnowval(:,j) = PDEparam2.termrewardfunc(wvec,svec(j),scale, param); % get value of stopping immediately.
         end
        
        eimprovement = abs( -g0 + delayfactor * Vvec / scale.beta - stopnowval / scale.beta);
        bigw = (length(wvec)-1)/2;
        %lowerindx = 0 * svec;  % allocate space to get lower and upper boundaries where one is indifferent between going forward and stopping.
        %upperindx = 0 * svec;
%        for j=1:length(svec)
            [~, lowerindx] = min(eimprovement(1:(bigw+1),:));
            [~, upperindx] = min(eimprovement((bigw+1):(2*bigw+1),:));
%        end       
        upperzero = wvec(bigw+upperindx)/scale.beta;
        lowerzero = wvec(lowerindx)/scale.beta;
        
        dw = wvec(2)-wvec(1);

        % try to find some good contour values for the contour plot
        maxbnd = 1.1*max(upvec)/scale.beta;
        minV = 10^(floor(log10(2*maxbnd)))/4;
       % V = [ minV/5000 minV/1000 minV/200 minV/50 minV/10 minV/4 minV/2 minV:minV:(2*maxbnd) ];
     	V2 = [ 0 minV/200 minV/50 minV/10 minV/4 minV/2 minV:2*minV:(4*maxbnd) ];

        % First, plot a contour plot of the benefit (above 0) of continuing, plus a
        % boundary of the upper and lower continuation set
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        [C, h]=contour(1./svec/scale.gamma,wvec/scale.beta,eimprovement,V2); 
        %[C, h]=contour(1./svec/scale.gamma,wvec/scale.beta,eimprovement); 
        clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        plot(1./svec/scale.gamma,upperzero,'-x','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        plot(1./svec/scale.gamma,lowerzero,'-o','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
%        plot(1./svec/scale.gamma,upvec/scale.beta,'-','LineWidth',2);set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.1*max(upperzero)+dw/scale.beta; tmp(3)=1.1*min(lowerzero)-dw/scale.beta; axis(tmp);
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        tmp = axis;
        if min(upperzero) > max(lowerzero)
            text(tmp(1) + (tmp(2)-tmp(1))/4, 0, 'Build simulation tool');
        end
        text(tmp(1) + (tmp(2)-tmp(1))/3, (tmp(4)+max(upperzero))/2, 'No simulation, implement alternative');
        text(tmp(1) + (tmp(2)-tmp(1))/3, (tmp(3)+min(lowerzero))/2, 'No simulation, keep status quo');
        xlabel('Certainty about unknown size of improvement, t_0 = \sigma^2 / \sigma_0^2','FontSize',myfontsize,'FontName','Times'); 
        ylabel('A priori expected benefit of improvement, \mu_0','FontSize',myfontsize,'FontName','Times')
        title('Improvement in E[NPV], | -g_0 + e^{-\delta u_0} V(\mu_0,t_0) - max(\mu_0,0)|','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-Like3-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    % Figure 4 for example 4
        td = 60; %time in days for analysis
        tddelayfactor = exp(-annualdisc * td / 365);
        r = td * 24 * 60 / eta; % number of replications which are possible
        tp = 1/26 ; % two weeks expressed in units of year, for delay of implementing alternative post selection
        tpdelayfactor = exp(-annualdisc * tp);
 
        t0vec = 1./svec/scale.gamma;
        sigvec = scale.sigma * sqrt( r ./ (t0vec .* (r + t0vec))) ;
        sigmat = ones(length(wvec),1) * sigvec;
        wmat = (wvec' / scale.beta) * ones(1,length(svec))  ;
%        evalfixeddeadline = tddelayfactor * sigmat .* PsiNorm( - wmat ./ sigmat ); 
        evalfixeddeadline = tddelayfactor * (sigmat .* PsiNorm( - wmat ./ sigmat )); 
%        evalfixeddeadline = max(0,evalfixeddeadline); % allow for immediate stopping here - the original paper did not do that
%        evalfixeddeadline = max(0,max(evalfixeddeadline,wmat)); % allow for immediate stopping here - the original paper did not do that
        delsequential = tpdelayfactor * Vvec / scale.beta - evalfixeddeadline;
        
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        MaxV = minV;
%        [C, h]=contour(1./svec/scale.gamma,wvec/scale.beta,eimprovement,V2); 
        V = [(2*MaxV:4*MaxV:10*MaxV)/10000 (2*MaxV:4*MaxV:10*MaxV)/1000 (2*MaxV:4*MaxV:10*MaxV)/100 (2*MaxV:4*MaxV:30*MaxV)/10 ];
        V2 = [-V V];
        [C, h]=contour(1./svec/scale.gamma,wvec/scale.beta,delsequential,V2); 
        clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        axis(tmp);
        PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
        xlabel('Certainty about unknown size of improvement, t_0 = \sigma^2 / \sigma_0^2','FontSize',myfontsize,'FontName','Times'); 
        ylabel('A priori expected benefit of improvement, \mu_0','FontSize',myfontsize,'FontName','Times')
        title('Value in flexible stopping, e^{-\delta t_d} V(\mu_0,t_0) - E[value with fixed]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-Like4-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    end 

    % extras
    fignum = PDEComparePDEboundApproxbound( fignum, pdeSoln );
    
disp('not yet implemented: C&G 2009 Table 1');
disp('not yet implemented: C&G 2009 Figure EC.1');

%%%%%%%%%%%%%%%%%%%%%%%%
    figout = fignum;
end