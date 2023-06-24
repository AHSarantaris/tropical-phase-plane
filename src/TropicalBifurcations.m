function h = TropicalBifurcations(columns,a,symbol,F,G,ulim,vlim,varargin)

na = size(a,1);
handle = tiledlayout(ceil(na/columns),columns,"Padding","tight","TileSpacing","compact");

for i = 1:na
    ax = nexttile;
    ai = a(i,:);
    TropicalPhasePlane(F(ai),G(ai),ulim,vlim)
    str = ai(1);
    if size(a,2) > 1
        str = "(" + str;
        for j = 2:size(a,2)
            str = str + "," + ai(j);
        end
        str = str + ")";
    end
%     title(ax,"{\boldmath$" + symbol +" = " + str + "$}",'Interpreter',"latex")
    title(ax,"$" + symbol +" = " + str + "$",'Interpreter',"latex")
    if nargin > 7
        SimulatedTrajectory(F(ai),G(ai),varargin{:})
    end
end
if nargout
    h=handle;
end
end

