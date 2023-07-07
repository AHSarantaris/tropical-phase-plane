function h = TropicalPhasePlane(F, G, ulim, vlim, showLabels)
%%% Options
datacursormode off
if ~exist('showLabels','var')
    showLabels = 0;
end
if ~exist('ulim','var')
    ulim = [];
end
if ~exist('vlim','var')
    vlim = [];
end

%%% Get systems
TDS = TropicalDynamicalSystem(F,G);
TDS_F = TropicalDynamicalSystem(F,[]);
TDS_G = TropicalDynamicalSystem(G,[]);

%%% Initialize axes
ax = gca;
axis square
hold(ax,'on')
ax.Units = 'pixels';
axsize = ax.Position(3:4);
ax.Units = 'normalized';
ax.TickDir = 'both';
fontsize(gcf,9,'points');

%%% Plot styling
colors = [0.6 0.7 1;
        0.8 0.9 1;
        1 0.5 0.5;
        1 0.75 0.6];
lineColor = 0*ones(1,3);
markerSize = 4;
tol = 1e-6;
fontSize = 8;

%%% Plot handles
nA = 0;
nE = 0;
arrowPlots = cell(1,TDS.nI);
equilibriumPlots = cell(1,TDS.nI);
labelHandles = cell(2,TDS.nI);

%%% Create default limits
Q = [TDS.points,TDS.nullclinePoints,TDS_F.points,TDS_G.points];
qlim = [min(Q,[],2), max(Q,[],2)];
nQ = size(Q,2);
if nQ > 1
    midU = sum(qlim(1,:))/2;
    midV = sum(qlim(2,:))/2;
    rangeU = diff(qlim(1,:));
    rangeV = diff(qlim(2,:));
else
    if nQ == 1
        midU = Q(1,1);
        midV = Q(2,1);
    else
        midU = 0;
        midV = 0;
    end
    rangeU = 1;
    rangeV = 1;
end
if rangeU == 0
    rangeU = rangeV;
end
if rangeV == 0
    rangeV = rangeU;
end
rangeMax = max(rangeU,rangeV);
if isempty(ulim)
    rangeU = rangeMax;
    ulim = midU + rangeU*[-1,1];
    ulim = [floor(ulim(1)*10)/10 ceil(ulim(2)*10)/10];
end
if isempty(vlim)
    rangeV = rangeMax;
    vlim = midV + rangeV*[-1,1];
    vlim = [floor(vlim(1)*10)/10 ceil(vlim(2)*10)/10];
end


%%% Polygon variables
corners = [ulim, ulim([2,1]); vlim, vlim];
vertices = corners;
vertexMap = zeros(TDS.nI,1);
nVertices = 4;
numVertices = zeros(1,TDS.nI);
polygons = cell(TDS.nI,1);
linePoints = zeros(2,2*size(TDS.points,2));
linePointMap = zeros(TDS.nI,TDS.nI,2);
nLP = 0;

%%% Determine polygon vertices for each region
for i = 1:TDS.nI-1
    for j = i+1:TDS.nI
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

%%% Plot regions
for i = 1:TDS.nI
    j = TDS.equalDegreeIndex(i);
    if j ~= 0 && j < i
        continue
    end
    V = vertices(:, vertexMap(i, 1:numVertices(i)));
    uniqueVertices = uniquePoints(V,tol);
    nUV = size(uniqueVertices,2);
    if nUV < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);        
    else
        convhullVertices = uniqueVertices(:,convhull(uniqueVertices'));
        polygons{i} = polyshape(convhullVertices(1,:), convhullVertices(2,:));
        d = TDS.flowVectors(:,i);
        if d(1) > 0
            c1 = colors(1,:);
        else
            c1 = colors(2,:);
        end
        if d(2) > 0
            c2 = colors(3,:);
        else
            c2 = colors(4,:);
        end
        if d(2) == 0
            c = c1;
        elseif d(1) == 0
            c = c2;
        else
            c = 0.9*[1 1 1];
        end
        hp = plot(polygons{i},'FaceColor',c,'EdgeColor','none','FaceAlpha',1);
        [Cu,Cv] = centroid(polygons{i});
        nA = nA + 1;
        arrowPlots{nA} = LineArrow([Cu;Cv],d,ulim,vlim,'k',1.3,true,true);
        s = sign(d);
        ds1 = 15;
        ds2 = 10;
        du = (-abs(s(2))*ds1 + s(1)*ds2) / axsize(1)*diff(ulim);
        dv = (abs(s(1))*ds1 + s(2)*ds2) / axsize(2)*diff(vlim);
        hc = plot(Cu-du,Cv-dv,'o','MarkerEdgeColor',0.2*ones(1,3),'MarkerFaceColor',0.98*ones(1,3),...
            'Visible',showLabels,'MarkerSize',fontSize+2,'LineWidth',1);
        ht = text(Cu-du,Cv-dv, "" +i, 'Visible',showLabels,'FontSize',fontSize,'FontWeight','bold',...
            'HorizontalAlignment','center');
        labelHandles(:,i) = {hc,ht};
%         hp.ButtonDownFcn = @(~,~) toggleVisibility(hc,ht) ;
    end
end


%%% Plot nullclines
PlotTropicalCurve(TDS,ulim,vlim,'Color','k','LineWidth',2)
PlotTropicalCurve(TDS_G,ulim,vlim,'Color','r','LineWidth',1)
PlotTropicalCurve(TDS_F,ulim,vlim,'Color',[0.1 0.4 1],'LineWidth',1)
for j = 1:2
    if j == 1
        TDS2 = TDS_F;
    else
        TDS2 = TDS_G;
    end
    for i = 1:size(TDS2.points,2)
        plot(TDS2.points(1,i),TDS2.points(2,i),'ko','MarkerFaceColor','k','MarkerSize',markerSize);
    end
end

%%% Plot lines
for i = 1:TDS.nI-1
    for j = i+1:TDS.nI
        if ~all(linePointMap(i,j,:))
            continue
        end
        Q = linePoints(:,linePointMap(i,j,:));
        I = TDS.equalDegreeIndex(i);
        J = TDS.equalDegreeIndex(j);
        if I ~= 0 || J ~= 0 && (I==0 || J==0)
            continue
        end
        type = TDS.manifoldType(i,j);
        manifoldFlow = TDS.manifoldDirection(:,TDS.manifoldMap(i,j));
        if type == TDS.TRANSVERSAL_FILIPPOV || type == TDS.TANGENTIAL_FILIPPOV
            qmid = [mean(Q(1,:)); mean(Q(2,:))];
            si = dot(manifoldFlow, TDS.flowVectors(:,i));
            sj = dot(manifoldFlow, TDS.flowVectors(:,j));
            if si < 0 || sj < 0
                manifoldFlow = -manifoldFlow;
            end
            nA = nA + 1;
            arrowPlots{nA} = LineArrow(qmid,manifoldFlow,ulim,vlim,'k',1.5,false,true);
        elseif type == TDS.TRANSVERSAL_NULLCLINE || type == TDS.TANGENTIAL_NULLCLINE
            direction = (i <= TDS.nU) + 1;
            if i <= TDS.nU
                ks = TDS.nU+1;
                ke = TDS.nI;
            else
                ks = 1;
                ke = TDS.nU;
            end
            K = ks:ke;
            Q = [Q(:,1) TDS.nullclinePoints(:,TDS.nullclinePointMap{i,j}) Q(:,2)];
            for iq = 2:size(Q,2)
                q1 = Q(:,iq-1);
                q2 = Q(:,iq);
                if all(abs(q1-q2) < tol)
                    continue
                end
                P = [q1 q2];
                pmid = [mean(P(1,:)); mean(P(2,:))];
                zp = TDS.alpha(K) + TDS.deg(K,:)*pmid;
                zmaxK = max(zp);
                [kmax,~] = find(abs(zp-zmaxK) < tol);
                kmax = K(kmax);
                sign_d = sum(sign(TDS.flowVectors(:,kmax)),2);
                if isempty(kmax) || (length(kmax) > 1 && sign_d(direction) == 0)
                    c = [1,1,1];
                else
                    if sign_d(direction) > 0
                        c = colors(1+2*(direction-1),:);
                    else
                        c = colors(2+2*(direction-1),:);
                    end
                    manifoldFlow = dot(sign_d,manifoldFlow)*manifoldFlow;
                    if ~all(manifoldFlow)
                        manifoldFlow = sign_d;
                    end
                    nA = nA + 1;
                    arrowPlots{nA} = FilledArrow(pmid,manifoldFlow,ulim,vlim,'k',0);
                    nA = nA + 1;
                    arrowPlots{nA} = FilledArrow(pmid,manifoldFlow,ulim,vlim,c,0);
                end
                line(P(1,:),P(2,:),'Color','k','LineWidth',4)
                line(P(1,:),P(2,:),'Color',c,'LineWidth',2)
            end
        end
    end
end
for i = 1:size(TDS.points,2)
    plot(TDS.points(1,i),TDS.points(2,i),'o',...
        'Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize)
end
for i = 1:nA
    uistack(arrowPlots{i},'top')
end

%%% Plot nullcline points
for i = 1:size(TDS.nullclinePoints,2)
    pt = TDS.nullclinePoints(:,i);
    if TDS.isEquilibrium(i)
        nE = nE + 1;
        equilibriumPlots{nE} = plot(pt(1),pt(2),'ko',...
            'MarkerFaceColor','w','MarkerSize',2*markerSize,'LineWidth',1.5);
    else
        plot(pt(1),pt(2),'o','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize);
    end
end
for i = 1:nE
    uistack(equilibriumPlots{i},'top')
end

%%% Draw labels on top
for i = 1:TDS.nI
    if ~isempty(labelHandles{1,i})
        uistack(labelHandles{1,i},'top')
        uistack(labelHandles{2,i},'top')
    end
end

%%% Set ticks at tropical points and nullcline points
ticks = zeros(2,2+size(TDS.points,2)+size(TDS.nullclinePoints,2)+size(TDS_F.points,2)+size(TDS_G.points,2));
iticks = 0;
for i = 1:size(TDS.points,2)
    iticks = iticks + 1;
    ticks(:,iticks) = TDS.points(:,i);
end
for i = 1:size(TDS.nullclinePoints,2)
    iticks = iticks + 1;
    ticks(:,iticks) = TDS.nullclinePoints(:,i);
end
for i = 1:size(TDS_F.points,2)
    iticks = iticks + 1;
    ticks(:,iticks) = TDS_F.points(:,i);
end
for i = 1:size(TDS_G.points,2)
    iticks = iticks + 1;
    ticks(:,iticks) = TDS_G.points(:,i);
end
ticks(:,[end-1 end]) = [ulim;vlim];
ticks = round(ticks*100)/100;
xticks(unique(ticks(1,:)))
yticks(unique(ticks(2,:)))


%%% Terminate
xlim(ulim)
ylim(vlim)
hold(ax,'off')
if nargout
    h = ax;
end

end





