function h = TropicalCurvePlot(F,ulim, vlim, varargin)

%%% Initialize variables
nF = size(F,1);
TDS = TropicalDynamicalSystem(F,zeros(0,4));
if isempty(TDS)
    return
end
T = F(:,2:4);
nP = size(TDS.points,2);


%%% Initialize axes
ax = gca;
hold(ax,'on')
ax.Units = 'pixels';
axsize = ax.Position(3:4);
ax.Units = 'normalized';
axmean = mean(axsize);

%%% Plot styling
lineColor = 0*ones(1,3);
lineWidth = 0.006*axmean;
markerSize = 3*min(1, lineWidth);
tol = 1e-6;


%%% Create default limits
if ~exist('ulim','var') || isempty(ulim) || ~exist('vlim','var')  || isempty(vlim)
    Q = [TDS.points,TDS.nullclinePoints];
    nQ = size(Q,2);
    if nQ > 1
        qlim = [min(Q,[],2), max(Q,[],2)];
        midU = sum(qlim(1,:))/2;
        midV = sum(qlim(2,:))/2;
        rangeU = diff(qlim(1,:));
        rangeV = diff(qlim(2,:));
        rangeU = max(rangeU,rangeV);
        rangeV = rangeU;
        if rangeU < tol && rangeV < tol
            rangeU = 1;
            rangeV = 1;
%         elseif rangeU < tol
%             rangeU = rangeV;
%         elseif rangeV < tol
%             rangeV = rangeU;
        end
    elseif nQ == 1
        midU = Q(1,1);
        midV = Q(2,1);
        rangeU = 1;
        rangeV = 1;
    end
    ulim = midU + rangeU*[-1,1];
    vlim = midV + rangeV*[-1,1];
end



%%% Polygon variables
corners = [ulim, ulim([2,1]); vlim, vlim];
vertices = corners;
vertexMap = zeros(nF,1);
nVertices = 4;
numVertices = zeros(1,nF);
linePoints = zeros(2,2*nP);
linePointMap = zeros(nF,nF,2);
nLP = 0;

%%% Find vertices
for i = 1:nF-1
    for j = i+1:nF
        if TDS.numEndpoints(i,j) == 2
            P = TDS.points(:, TDS.pointMap(i,j,1:2));
            if all(abs(diff(P,1,2)) < tol)
                continue
            end
            Q = pointsInside(P,ulim,vlim);
            if isempty(Q)
                continue
            end
            nLP = nLP + 2;
            linePoints(:,[nLP-1,nLP]) = Q;
            linePointMap(i,j,:) = [nLP-1,nLP];
            nVertices = nVertices + 2;
            vertices(:,[nVertices-1,nVertices]) = Q;
            numVertices(i) = numVertices(i) + 2;
            numVertices(j) = numVertices(j) + 2;
            vertexMap(i, [numVertices(i)-1,numVertices(i)]) = [nVertices-1,nVertices];
            vertexMap(j, [numVertices(j)-1,numVertices(j)]) = [nVertices-1,nVertices];
        elseif TDS.numEndpoints(i,j) == 1
            p1 = TDS.points(:, TDS.pointMap(i,j,1));
            p2 = boundaryPoint(p1, TDS.manifoldDirection(:,TDS.manifoldMap(i,j)), ulim, vlim);
            Q = pointsInside([p1 p2],ulim,vlim);
            if isempty(Q)
                continue
            end
            T2 = T(j,:)*[1;Q(:,2)];
            Tmax = max(T*[1;Q(:,2)]);
            if T2 + tol >= Tmax
                nLP = nLP + 2;
                linePoints(:,[nLP-1,nLP]) = Q;
                linePointMap(i,j,:) = [nLP-1,nLP];
                nVertices = nVertices + 2;
                vertices(:,[nVertices-1 nVertices]) = Q;
                numVertices(i) = numVertices(i) + 2;
                numVertices(j) = numVertices(j) + 2;
                vertexMap(i, [numVertices(i)-1 numVertices(i)]) = [nVertices-1 nVertices];
                vertexMap(j, [numVertices(j)-1 numVertices(j)]) = [nVertices-1 nVertices];
            end
        end
    end
end

%%% Add corners as vertices
for ic = 1:4
    Tc = T*[1;corners(:,ic)];
    Tmax = max(Tc);
    [jmax,~] = find(abs(Tc-Tmax) < tol);
    for j = jmax'
        numVertices(j) = numVertices(j) + 1;
        vertexMap(j, numVertices(j)) = ic;
    end
end

%%% Plot dominating areas
for i = 1:nF
    V = vertices(:, vertexMap(i, 1:numVertices(i)));
    uniqueVertices = uniquePoints(V,tol);
    nUV = size(uniqueVertices,2);
    if nUV < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);
    end
end

%%% Plot lines
for i = 1:nF-1
    for j = i+1:nF
        lc = TDS.manifoldType(i,j);
        if ~all(linePointMap(i,j,:))
            continue
        end
        if lc == TDS.CROSSING
            lineStyle = '--';
        else
            lineStyle = '-';
        end
        Q = linePoints(:,linePointMap(i,j,:));
        line(Q(1,:),Q(2,:),'Color',lineColor,'LineWidth',lineWidth,'LineStyle',lineStyle,varargin{:})
    end
end
for i = 1:nP
    plot(TDS.points(1,i),TDS.points(2,i),'o','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize);
end



%%% Terminate
xlim(ulim)
ylim(vlim)
% hold(ax,'off')
if nargout
    h = ax;
end

end


