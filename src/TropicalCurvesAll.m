function TropicalCurvesAll(F,G,ulim,vlim)
if ~exist('ulim','var') || isempty(ulim) || ~exist('vlim','var')  || isempty(vlim)
    TropicalPhasePlane(F,G)
    ulim = xlim;
    vlim = ylim;
else
    TropicalPhasePlane(F,G,ulim,vlim)
end
TropicalCurvePlot(F,ulim,vlim,'Color','b')
TropicalCurvePlot(G,ulim,vlim,'LineStyle','--','Color','r')
end