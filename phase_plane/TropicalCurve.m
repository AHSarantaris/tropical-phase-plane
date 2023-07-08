function h = TropicalCurve(TDS, ulim, vlim, varargin)

tol = 1e-6;
nI = TDS.nI;

%%% Initialize axes
ax = gca;
hold(ax,'on')

%%% Create default limits
points = [TDS.points,TDS.nullclinePoints];
[ulim,vlim] = getLimits(points,ulim,vlim);



%%% Polygon variables
polygons = getPolygonVertices(TDS,ulim,vlim);

%%% Plot lines
for i = 1:nI-1
    for j = i+1:nI
        lc = TDS.manifoldType(i,j);
        if ~all(polygons.linePointMap(i,j,:))
            continue
        end
        if TDS.nV == 0
            if TDS.delta(i)*TDS.delta(j) < 0
                lineStyle = '-';
            else
                lineStyle = '--';
            end
        elseif lc == TDS.CROSSING
            if dot(TDS.flowVectors(:,i),TDS.flowVectors(:,j)) ~= 0
                continue
            else
                lineStyle = '--';
            end
        else
            lineStyle = '-';
        end
        Q = polygons.linePoints(:,polygons.linePointMap(i,j,:));
        line(Q(1,:),Q(2,:),'LineStyle',lineStyle,varargin{:})
    end
end

%%% Terminate
xlim(ulim)
ylim(vlim)
if nargout
    h = ax;
end

end


