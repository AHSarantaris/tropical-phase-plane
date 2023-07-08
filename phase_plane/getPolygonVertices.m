function polygons = getPolygonVertices(TDS,ulim,vlim)

tol = 1e-6;
nI = TDS.nI;

%%% Polygon variables
corners = [ulim, ulim([2,1]); vlim, vlim];
vertices = corners;
vertexMap = zeros(nI,1);
numVertices = zeros(1,nI);
linePoints = zeros(2,2*size(TDS.points,2));
linePointMap = zeros(nI,nI,2);
uniqueVertices = cell(1,nI);
nUniqueVertices = zeros(1,nI);

%%% Counters
nVertices = 4;
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

for i = 1:nI
    V = vertices(:, vertexMap(i, 1:numVertices(i)));
    uniqueVertices{i} = uniquePoints(V,tol);
    nUniqueVertices(i) = size(uniqueVertices{i},2);
    if nUniqueVertices < 3
        linePointMap(i,:,:) = 0*linePointMap(i,:,:);
        linePointMap(:,i,:) = 0*linePointMap(:,i,:);
    end
end

polygons.linePoints = linePoints;
polygons.linePointMap = linePointMap;
polygons.vertices = uniqueVertices;

end


%% Auxiliary function 1

function U = uniquePoints(P,tol)
    U = zeros(size(P));
    nP = size(P,2);
    nU = 0;
    for i = 1:nP
        dP = abs(P(:,i) - P(:,i+1:end));
        if sum(all(dP < tol)) == 0
            nU = nU + 1;
            U(:,nU) = P(:,i);
        end
    end
    U = U(:,1:nU);
end

%% Auxiliary function 2

function p2 = boundaryPoint(p1, r, ulim, vlim)
if (p1(1) < ulim(1) && r(1) < 0) ...
        || (p1(1) > ulim(2) && r(1) > 0) ...
        || (p1(2) < vlim(1) && r(2) < 0) ...
        || (p1(2) > vlim(2) && r(2) > 0)
    p2 = p1;
    return
end
if r(2) == 0
    if r(1) < 0
        p2 = [ulim(1); p1(2)];
    else
        p2 = [ulim(2); p1(2)];
    end
    return
elseif r(1) == 0
    if r(2) < 0
        p2 = [p1(1); vlim(1)];
    else
        p2 = [p1(1); vlim(2)];
    end
    return
end
if r(1) < 0
    t1 = (ulim(1) - p1(1)) / r(1);
else
    t1 = (ulim(2) - p1(1)) / r(1);
end
if r(2) < 0
    t2 = (vlim(1) - p1(2)) / r(2);
else
    t2 = (vlim(2) - p1(2)) / r(2);
end
[tmin,~] = min([t1,t2]);
p2 = p1 + tmin * r;
end

%% Auxiliary function 3

function Q = pointsInside(P, ulim, vlim)
    Q = [];
    p1 = P(:,1);
    p2 = P(:,2);
    pos1 = position(p1);
    pos2 = position(p2);
    if pos1 == 5 && pos2 == 5
        Q = [p1 p2];
    elseif pos1 == 5
        Q(:,1) = p1;
        Q(:,2) = intersectionFromInside(p1,p2,pos2);
    elseif pos2 == 5
        Q(:,1) = intersectionFromInside(p2,p1,pos1);
        Q(:,2) = p2;
    else
        qu1 = intersection(1,p1,p2,ulim(1),vlim);
        qu2 = intersection(1,p1,p2,ulim(2),vlim);
        qv1 = intersection(2,p1,p2,vlim(1),ulim);
        qv2 = intersection(2,p1,p2,vlim(2),ulim);
        Q4 = {qu1 qu2 qv1 qv2};
        for i = 1:3
            for j = i+1:4
                if ~isempty(Q4{i}) && ~isempty(Q4{j})
                    Q = cell2mat(Q4([i j]));
                end                
            end
        end
    end
    
    function q = intersection(ip,p1,p2,line,lim)
        r = p2 - p1;
        if r(ip) == 0
            q = [];
            return
        end
        tp = norm(r);
        tq = (line - p1(ip)) / r(ip);
        q = p1 + tq * r;
        iq = mod(ip,2) + 1;
    %         if tq < 0 || tp < tq || q(iq) < lim(1) || lim(2) < q(iq)
    %             q = [];
    %         end
    end
    
    function q = intersectionFromInside(pi,po,pos)
        if pos == 4
            q = intersection(1,pi,po,ulim(1),vlim);
        elseif pos == 6
            q = intersection(1,pi,po,ulim(2),vlim);
        elseif pos == 2
            q = intersection(2,pi,po,vlim(1),ulim);
        elseif pos == 8
            q = intersection(2,pi,po,vlim(2),ulim);
        elseif pos == 1
            q1 = intersection(1,pi,po,ulim(1),vlim);
            q2 = intersection(2,pi,po,vlim(1),ulim);
            if q1(2) > vlim(1)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 3
            q1 = intersection(1,pi,po,ulim(2),vlim);
            q2 = intersection(2,pi,po,vlim(1),ulim);
            if q1(2) > vlim(1)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 7
            q1 = intersection(1,pi,po,ulim(1),vlim);
            q2 = intersection(2,pi,po,vlim(2),ulim);
            if q1(2) < vlim(2)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 9
            q1 = intersection(1,pi,po,ulim(2),vlim);
            q2 = intersection(2,pi,po,vlim(2),ulim);
            if q1(2) < vlim(2)
                q = q1;
            else
                q = q2;
            end
        else
            error('Third argument must be between 1 and 9 but not 5')
        end
    end
    
    function pos = position(p)
        if p(2) < vlim(1)
            if p(1) < ulim(1)
                pos = 1;
            elseif p(1) > ulim(2)
                pos = 3;
            else
                pos = 2;
            end
        elseif p(2) > vlim(2)
            if p(1) < ulim(1)
                pos = 7;
            elseif p(1) > ulim(2)
                pos = 9;
            else
                pos = 8;
            end
        else
            if p(1) < ulim(1)
                pos = 4;
            elseif p(1) > ulim(2)
                pos = 6;
            else
                pos = 5;
            end  
        end
    end

end
