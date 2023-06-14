function TropicalCurvesAll(F,G,ulim,vlim)
TropicalPhasePlane(F,G,ulim,vlim)
TropicalCurvePlot(F,ulim,vlim,'Color','b')
TropicalCurvePlot(G,ulim,vlim,'LineStyle','--','Color','r')
end