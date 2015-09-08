function [ rval, figout, pdeSolnStruct ] = DoCGPlots( fignum, basename )

% look in default location for base name of input files for discounted
% reward case, otherwise use the argument passed to this function.
rval = 0;               % assume failure to load unless loading is completed
PDESolnStructure = [];  % default to empty solution structure
[ST,I] = dbstack;
routinename = ST.name;  % get function name

if nargin < 1
    basename = 'Matfiles\CG';   % give location
end

[rval, pdeSolnStruct] = PDESolnLoad(basename);

if ~rval
    warning(routinename, ': unable to open ',basename);
	return; 
end

myfontsize=16;
mysmallfontsize=14;
points = 144*3; %spacing between labels on contours - made so that only one label appears per line
fracheight = 0.9;   % take up 90% of screen
square = true;      % make plot format to be 'square' if true

% do a bit of input parameter checking
figdir = pdeSolnStruct.Header.PDEparam.figdir;
figsave = true;

% create the directory for the figures if it does not exist already and the figures are
% to be saved
if ~isdir(figdir) & figsave
    mkdir(figdir);
end

pdeSolnStruct.Header;
PDEscale=pdeSolnStruct.Header.PDEscale;
PDEparam=pdeSolnStruct.Header.PDEparam;

firsts = pdeSolnStruct.Header.firsts;
lasts = pdeSolnStruct.Header.lasts;
lasthelds = pdeSolnStruct.Header.lasthelds;
fName = pdeSolnStruct.Header.fName;

alpha = pdeSolnStruct.Header.PDEscale.alpha;
beta = pdeSolnStruct.Header.PDEscale.beta;
gamma = pdeSolnStruct.Header.PDEscale.gamma;
approxmeth = pdeSolnStruct.Header.PDEparam.approxmethod;


%% Generate plots similar to those in the Electronic companion
% First, figures EC3, EC.4 and EC.5 (adapting to the computations in the local
% install for the PDE solution)
%    for ijk = pdeSolnStruct.Header.StartFileVal:pdeSolnStruct.Header.EndFileVal
    for ijk = 1:1
        block = ijk;
 
        svec = pdeSolnStruct.Data(ijk).svec;
        wvec = pdeSolnStruct.Data(ijk).wvec;

        upvec = pdeSolnStruct.Data(ijk).upvec;
        downvec = pdeSolnStruct.Data(ijk).downvec;
        up1 = pdeSolnStruct.Data(ijk).up1;
        down1 = pdeSolnStruct.Data(ijk).down1;

%        Vvec = stopfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam); % compute value of stopping immediately
        Bwsmatrix = pdeSolnStruct.Data(ijk).Bwsmatrix;
        ENwsmatrix = pdeSolnStruct.Data(ijk).ENwsmatrix;
        EPCSwsmatrix = pdeSolnStruct.Data(ijk).EPCSwsmatrix;
        Vvec = zeros(size(Bwsmatrix)); 
        for j=1:length(svec)
            vveccol  = PDEGetVals(pdeSolnStruct,wvec,svec(j)) ;
            Vvec(:,j) = vveccol;
        end
        
        dw = wvec(2)-wvec(1);

        maxup1 = max(up1);
        maxupvec = max(upvec); 
        maxwvec = max(wvec);

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
        tmp=axis;tmp(4)=1.2*max(max(up1),max(upvec)); tmp(3)=1.05*min(min(down1),min(downvec)); tmp(3)=-1.2*tmp(4); axis(tmp);
        xlabel('Scaled time, \tau','FontSize',myfontsize,'FontName','Times'); 
        ylabel('Scaled mean, z_\tau / \tau','FontSize',myfontsize,'FontName','Times')
        title('Value function B_1','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(fName,'Fig-EC345-',int2str(ijk));
        PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

    end 

%% Figure EC.6 A and B (at least for boundary: the B_1(0,s) is not as important here, supplemented with 'improved' analytical approximation for boundary
accumsvec = pdeSolnStruct.Computed.accumsvec;
accumupper = pdeSolnStruct.Computed.accumupper;
accumlower = pdeSolnStruct.Computed.accumlower;
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

% For figure 1: 
scalevec = {'c', 0, 'sigma', 1e7, 'discrate', 20*0.1/365/24/60, 'P', 1 };
CGparamvec = { 't0', 1, 'tEND', 10e5, 'precfactor', 6, 'BaseFileName', 'CG' };
paramvec = [CGparamvec, CGfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )

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
ylabel('Output mean, y_t / t','FontSize',myfontsize,'FontName','Times')
mytitle = strcat(fName,'Fig1');
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );

c=1;
wtest = (c/scale.discrate) * scale.beta;
stest = PDEInvBound( upperDisc, wtest );
numrepstest = 1 ./ (scale.gamma * stest);

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
hold off
semilogx(1./(scale.gamma*accumsvec),accumupper/scale.beta,'--');
hold on
semilogx(1./(scale.gamma*accumsvec),accumupper/scale.beta - (c/scale.discrate),'-');
    plot( [1 numrepstest], [0 0],'-.');
    plot( [numrepstest numrepstest], [0 -1e6],'-xk');    
legend('Implement (c=0), \beta^{-1} b_1(1/\gamma t)','Implement (c=0), \beta^{-1} b_1(1/\gamma t)-c/\gamma','mean 0','max reps');
PDEUtilStdizeFigure( fignum, fracheight, mysmallfontsize, square );
tmp=axis; tmp(1)=1; tmp(2)=2e5; tmp(4)=2e6; tmp(3)=-1e6; axis(tmp);


xlabel('Number of replications, t','FontSize',myfontsize,'FontName','Times'); 
ylabel('Output mean, y_t / t','FontSize',myfontsize,'FontName','Times')
mytitle = strcat(fName,'Fig1');
PDEUtilSaveFigEpsPdf( fignum, figdir, mytitle );



%%%%%%%%%%%%%%%%%%%%%%%%
    figout = fignum;
end