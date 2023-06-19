function h = TropicalCurvePlot2(F, G,ulim, vlim, varargin)

%%% Initialize variables
nF = size(F,1);
if isempty(G)
    G = zeros(0,4);
end
nG = size(G,1);
TL = TropicalCurveData(F,G);
if isempty(TL)
    return
end
T = [F(:,2:4); G(:,2:4)];
dT = [F(:,1)', zeros(1,nG); zeros(1,nF), G(:,1)'];
nT = nF + nG;
nP = size(TL.points,2);


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

%%% Line classification constants
BENT_CROSSING = 1;
STRAIGHT_CROSSING = 2;
SLIDING = 3;
NULLCLINE = 4;
CROSS_N_SLIDE = 5;

%%% Create default limits
if ~exist('ulim','var') || isempty(ulim) || ~exist('vlim','var')  || isempty(vlim)
    Q = [TL.points,TL.tropicalPoints];
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
vertexMap = zeros(nT,1);
nVertices = 4;
numVertices = zeros(1,nT);
linePoints = zeros(2,2*nP);
linePointMap = zeros(nT,nT,2);
nLP = 0;

%%% Find vertices
for i = 1:nT-1
    for j = i+1:nT
        if TL.numPoints(i,j) == 2
            P = TL.points(:, TL.pointMap(i,j,1:2));
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
        elseif TL.numPoints(i,j) == 1
            p1 = TL.points(:, TL.pointMap(i,j,1));
            p2 = OuterPoint(p1, TL.vectors(:,TL.vectorMap(i,j)), ulim, vlim);
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
for i = 1:nT
    V = vertices(:, vertexMap(i, 1:numVertices(i)));
    uniqueVertices = uniquePoints(V,tol);
    nUV = size(uniqueVertices,2);
    if nUV < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);
    end
end

%%% Plot lines
for i = 1:nT-1
    for j = i+1:nT
        lc = TL.lineClassification(i,j);
        if ~all(linePointMap(i,j,:))
            continue
        end
        if nG == 0
            if F(i,1)*F(j,1) < 0
                lineStyle = '-';
            else
                lineStyle = '--';
            end
        else
            n_ij = T(j,2:3) - T(i,2:3);
            di = TL.flowVectors(:,i);
            dj = TL.flowVectors(:,j);
            if (n_ij*di)*(n_ij*dj) > 0 && isempty(G)
                lineStyle = '--';
            else
                lineStyle = '-';
            end
        end
        Q = linePoints(:,linePointMap(i,j,:));
        line(Q(1,:),Q(2,:),'Color',lineColor,'LineWidth',lineWidth,'LineStyle',lineStyle,varargin{:})
        
    end
end
% for i = 1:nP
%     plot(TL.points(1,i),TL.points(2,i),'o','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize);
% end



%%% Terminate
xlim(ulim)
ylim(vlim)
% hold(ax,'off')
if nargout
    h = ax;
end

end


