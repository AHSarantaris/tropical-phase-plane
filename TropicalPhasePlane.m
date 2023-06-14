function h = TropicalPhasePlane(F, G, ulim, vlim, showLabels)
%%% Options
datacursormode off
if ~exist('showLabels','var')
    showLabels = 0;
end

%%% Initialize variables
nF = size(F,1);
nG = size(G,1);
nT = nF + nG;
TL = TropicalCurves(F,G);
T = [F(:,2:4); G(:,2:4)];
dT = [F(:,1)', zeros(1,nG); zeros(1,nF), G(:,1)'];
nV = size(TL.vectors,2);
nP = size(TL.points,2);

%%% Line classification
BENT_CROSSING = 1;
STRAIGHT_CROSSING = 2;
SLIDING = 3;
NULLCLINE = 4;
CROSS_N_SLIDE = 5;

%%% Initialize axes
ax = gca;
hold(ax,'on')
ax.Units = 'pixels';
axsize = ax.Position(3:4);
ax.Units = 'normalized';
axmean = mean(axsize);
ax.TickDir = 'both';

%%% Plot styling
colors = [0.6 0.7 1; 0.8 0.9 1; 1 0.5 0.5; 1 0.75 0.6];
lineColor = 0*ones(1,3);
lineWidth = 0.006*axmean;
tol = 1e-6;
fontsize = 0.02*axmean;

%%% Plot handles
labelHandles = cell(2,nT);

%%% Create default limits
if ~exist('vlim','var') && exist('ulim','var')
    vlim = ulim;
elseif ~exist('ulim','var') || ~exist('vlim','var')
    Q = [TL.points,TL.tropicalPoints];
    nQ = size(Q,2);
    if nQ > 1
        qlim = [min(Q,[],2), max(Q,[],2)];
        midU = sum(qlim(1,:))/2;
        midV = sum(qlim(2,:))/2;
        rangeU = diff(qlim(1,:));
        rangeV = diff(qlim(2,:));
        if rangeU < tol && rangeV < tol
            rangeU = 1;
            rangeV = 1;
        elseif rangeU < tol
            rangeU = rangeV;
        elseif rangeV < tol
            rangeV = rangeU;
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
polygons = cell(nT,1);
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
    uniqueVertices = UniquePoints(V,tol);
    nUV = size(uniqueVertices,2);
    if nUV < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);        
    else
        convhullVertices = uniqueVertices(:,convhull(uniqueVertices'));
        polygons{i} = polyshape(convhullVertices(1,:), convhullVertices(2,:));
        [~,imax] = max(abs(dT(:,i)));
        if imax == 1
            if dT(imax,i) > 0
                c = colors(1,:);
            else
                c = colors(2,:);
            end
        else
            if dT(imax,i) > 0
                c = colors(3,:);
            else
                c = colors(4,:);
            end
        end
        hp = plot(polygons{i},'FaceColor',c,'EdgeColor','none','FaceAlpha',1);
        [Cu,Cv] = centroid(polygons{i});
        ha = FilledArrow([Cu;Cv],dT(:,i),ulim,vlim,'k',1);
        s = sign(dT(:,i));
        ds1 = 15;
        ds2 = 10;
        du = (-abs(s(2))*ds1 + s(1)*ds2) / axsize(1)*diff(ulim);
        dv = (abs(s(1))*ds1 + s(2)*ds2) / axsize(2)*diff(vlim);
        hc = plot(Cu-du,Cv-dv,'o','MarkerEdgeColor',0.2*ones(1,3),'MarkerFaceColor',0.98*ones(1,3),...
            'Visible',showLabels,'MarkerSize',fontsize+2);
        ht = text(Cu-du,Cv-dv+0.002*diff(vlim), "" +i, 'Visible',showLabels,'FontSize',fontsize,'FontWeight','bold',...
            'HorizontalAlignment','center');
        hp.ButtonDownFcn = @(~,~) toggleVisibility(hc,ht) ;
        labelHandles(:,i) = {hc,ht};
        %{
%         if i <= nF
%             hc.DataTipTemplate.DataTipRows(1) = dataTipTextRow("F["+i+"] = ",T(i,:));
%         else
%             hc.DataTipTemplate.DataTipRows(1) = dataTipTextRow("G["+(i+nF)+"] = ",T(i,:));
%         end
%         hc.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Term: ',T(i,:));
        %}
    end
end

%%% Plot lines
arrowPlots = cell(1,nT);
plotTropicalPoints = ones(1,size(TL.tropicalPoints,2));
nA = 0;
for i = 1:nT-1
    for j = i+1:nT
        lc = TL.lineClassification(i,j);
        if ~all(linePointMap(i,j,:))
            continue
        end
        Q = linePoints(:,linePointMap(i,j,:));
        if lc ~= NULLCLINE
            if lc == SLIDING
                lineStyle = '-';
            elseif lc == CROSS_N_SLIDE
                lineStyle = '-';
            else
                lineStyle = '--';
            end
%             lineStyle = '-';
            line(Q(1,:),Q(2,:),'LineStyle',lineStyle,'Color',lineColor,'LineWidth',lineWidth);
        else % Nullcline
            qmid = [mean(Q(1,:)); mean(Q(2,:))];
            direction = (i <= nF) + 1;
            if i <= nF
                ks = nF+1;
                ke = nT;
            else
                ks = 1;
                ke = nF;
            end
            K = ks:ke;
            nPE = TL.numTropicalPoints(i,j);
            Q = [Q(:,1) TL.tropicalPoints(:,TL.tropicalPointMap(i,j,1:nPE)) Q(:,2)];
            for iq = 2:size(Q,2)
                q1 = Q(:,iq-1);
                q2 = Q(:,iq);
                if all(abs(q1-q2) < tol)
                    continue
                end
                P = [q1 q2];
                pmid = [mean(P(1,:)); mean(P(2,:))];
                Tp = T(K,:)*[1;pmid];
                TmaxK = max(Tp);
                [kmax,~] = find(abs(Tp-TmaxK) < tol);
                kmax = K(kmax);
                dTsign = sum(sign(dT(:,kmax)),2);
                if isempty(kmax) || (length(kmax) > 1 && dTsign(direction) == 0)
                    c = [1,1,1];
                    plotTropicalNullcline(P,c)
                    if nPE
                        plotTropicalPoints(TL.tropicalPointMap(i,j,1:nPE)) = zeros(1,nPE);
                    end
                else
                    if dTsign(direction) > 0
                        c = colors(1+2*(direction-1),:);
                    else
                        c = colors(2+2*(direction-1),:);
                    end
                    dV = [T(i,3)-T(j,3); T(j,2)-T(i,2)];
                    dV = dot(dTsign,dV)*dV;
                    if ~all(dV)
                        dV = dTsign;
                    end
                    plotTropicalNullcline(P,c)
                    nA = nA + 1;
                    arrowPlots{nA} = FilledArrow(pmid,dV,ulim,vlim,c,0);
                end
            end
        end
    end
end
markerSize = min(3, 3*lineWidth);
for i = 1:nP
    plot(TL.points(1,i),TL.points(2,i),'o','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize)
end
for i = 1:nA
    uistack(arrowPlots{i},'top')
end
for i = 1:nT
    if ~isempty(labelHandles{1,i})
        uistack(labelHandles{1,i},'top')
        uistack(labelHandles{2,i},'top')
    end
end
%%% Plot tropical points
for i = 1:size(TL.tropicalPoints,2)
    if plotTropicalPoints(i)
        pt = TL.tropicalPoints(:,i);
        if TL.isEquilibrium(i)
            plot(pt(1),pt(2),'ko','MarkerFaceColor','w','MarkerSize',2.5*markerSize,'LineWidth',1.5);
        else
            plot(pt(1),pt(2),'o','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',markerSize);
        end
%         plot(pe(1),pe(2),'k.','MarkerSize',1.2*fontsize);
    end
end

%%% Terminate
xlim(ulim)
ylim(vlim)
hold(ax,'off')
if nargout
    h = ax;
end
end


function plotTropicalNullcline(Q,c)
%     ax = gca;
%     ax.Units = 'pixels';
%     line(Q(1,:),Q(2,:),'Color',0*ones(1,3),'LineWidth',0.006*mean(ax.Position(3:4))) % uncomment
    
    line(Q(1,:),Q(2,:),'Color',0*ones(1,3),'LineWidth',4)
    line(Q(1,:),Q(2,:),'Color',c,'LineWidth',2.5)
end

