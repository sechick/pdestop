function [ rval ] = PDEUtilStdizeFigure( fignum, fracheight, smallfontsize, square )
%UtilStdizeFigure standardizes the size of a given figure in order to help
%make printing from figures more standardized from one figure to the next.
% square is true of the shape of the graphic should be square (the
% default), otherwise the graphic shape is not modified

if nargin < 3
    square = true;
end

scnsize = get(0,'ScreenSize');
h = figure(fignum);
%axhand = axes();

position = get(h,'Position');
outerpos = get(h,'OuterPosition');
borders = outerpos - position;

if square
    axis('square');
    tmp = min(scnsize(3),scnsize(4));
    tmp2 = fracheight*tmp;
    pos1 = [scnsize(3) * (1 - fracheight)/2,...
        scnsize(4) * (1 - fracheight)/2,...
        scnsize(3) * (1 - fracheight)/2 + tmp2,...
        scnsize(4) * (1 - fracheight)/2 + tmp2] ;
else
    pos1 = [scnsize(3) * (1 - fracheight)/2,...
        scnsize(4) * (1 - fracheight)/2,...
        scnsize(3) * (1 + fracheight)/2,...
        scnsize(4) * (1 + fracheight)/2];
end

set(h,'OuterPosition',pos1) 
set(gca,'fontsize',smallfontsize)

rval = true;

end
