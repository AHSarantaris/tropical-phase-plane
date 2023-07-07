function h = PlotTropicalCurve(TDS, ulim, vlim, varargin)

nI = TDS.nI;

%%% Initialize axes
ax = gca;
hold(ax,'on')

%%% Plot styling
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
vertexMap = zeros(nI,1);
nVertices = 4;
numVertices = zeros(1,nI);
linePoints = zeros(2,2*size(TDS.points,2));
linePointMap = zeros(nI,nI,2);
nLP = 0;

%%% Determine polygon vertices for each region
for i = 1:nI-1
    for j = i+1:nI
        if TDS.numEndpoints(i,j) == 2
            P = TDS.points(:, TDS.pointMap(i,j,1:2));
            if all(abs(diff(P,1,2)) < tol)
                continue
            end
            Q = pointsInside(P,ulim,vlim);
            if isempty(Q)
                continue
            end
        elseif TDS.numEndpoints(i,j) == 1
            p1 = TDS.points(:, TDS.pointMap(i,j,1));
            p2 = boundaryPoint(p1, TDS.manifoldDirection(:,TDS.manifoldMap(i,j)), ulim, vlim);
            Q = pointsInside([p1 p2],ulim,vlim);
            if isempty(Q)
                continue
            end
            z = TDS.alpha + TDS.deg*Q(:,2);
            if z(j) + tol < max(z)
                continue
            end
        elseif TDS.manifoldMap(i,j) ~= 0
            manifoldDirection = TDS.manifoldDirection(:,TDS.manifoldMap(i,j));
            if manifoldDirection(1) == 0
                vmid = mean(vlim);
                umid = (TDS.alpha(j)-TDS.alpha(i))/(TDS.deg(i,1)-TDS.deg(j,1));
            else
                umid = mean(ulim);
                vmid = ((TDS.deg(j,1)-TDS.deg(i,1))*umid + TDS.alpha(j)-TDS.alpha(i))/(TDS.deg(i,2)-TDS.deg(j,2));
            end
            pmid = [umid; vmid];
            p1 = boundaryPoint(pmid, manifoldDirection, ulim, vlim);
            p2 = boundaryPoint(pmid, -manifoldDirection, ulim, vlim);
            Q = pointsInside([p1 p2],ulim,vlim);
            if isempty(Q)
                continue
            end
        else
            continue
        end
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

%%% Add corners as vertices
for ic = 1:4
    zc = TDS.alpha + TDS.deg*corners(:,ic);
    [jmax,~] = find(abs(zc-max(zc)) < tol);
    for j = jmax'
        numVertices(j) = numVertices(j) + 1;
        vertexMap(j, numVertices(j)) = ic;
    end
end

%%% Plot dominating areas
for i = 1:nI
    V = vertices(:, vertexMap(i, 1:numVertices(i)));
    uniqueVertices = uniquePoints(V,tol);
    nUV = size(uniqueVertices,2);
    if nUV < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);
    end
end

%%% Plot lines
for i = 1:nI-1
    for j = i+1:nI
        lc = TDS.manifoldType(i,j);
        if ~all(linePointMap(i,j,:))
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
        Q = linePoints(:,linePointMap(i,j,:));
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


