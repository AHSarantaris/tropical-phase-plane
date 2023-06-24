function TDS = TropicalDynamicalSystem(F,G)


%%% Initialize
if isempty(F)
    F = zeros(0,4);
end
if isempty(G)
    G = zeros(0,4);
end
tol = 1e-12;
nU = size(F,1);
nV = size(G,1);
nI = nU + nV;


%%% Throw error if the flow vector of a monomial is the zero vector 
for i = 1:nU
    if abs(F(i,1)) < tol
        error('The sign of monomial F[%d] is zero.',i)
    end
end
for i = 1:nV
    if abs(G(i,1)) < tol
        error('The sign of monomial G[%d] is zero.',i)
    end
end
%%% Throw error if the monomials are not unique for each polynomial
for i = 1:nU-1
    for j = i+1:nU
        if isequal(F(i,2:4), F(j,2:4))
            error('The monomials in F are not unique: F[%d] and F[%d] are equal.\n',i,j)
        end
    end
end
for i = 1:nV-1
    for j = i+1:nV
        if isequal(G(i,2:4), G(j,2:4))
            error('The monomials in G are not unique: G[%d] and G[%d] are equal.\n',i,j)
        end
    end
end


%%% Initialize Counters
nM = 0;     % Number of manifold direction vectors
nP = 0;     % Number of tropical points
nNP = 0;    % Number of nullcline points


%%% Output structure
% Monomial variables
TDS.delta = [F(:,1); G(:,1)];
TDS.alpha = [F(:,2); G(:,2)];
TDS.deg = [F(:,3:4); G(:,3:4)];   
% Number of monomials
TDS.nU = size(F,1);
TDS.nV = size(G,1);
TDS.nI = nU + nV;
% Enumeration of manifold type
TDS.CROSSING = 1;
TDS.TRANSVERSAL_FILIPPOV = 2;
TDS.TANGENTIAL_FILIPPOV = 3;
TDS.TRANSVERSAL_NULLCLINE = 4;
TDS.TANGENTIAL_NULLCLINE = 5;
% Flow vectors
TDS.flowVectors = zeros(2,nI);
TDS.equalDegreeIndex = zeros(1,nI);
% Tropical manifolds
TDS.manifoldDirection = zeros(2,nI);
TDS.manifoldMap = zeros(nI,nI);
TDS.numEndpoints = zeros(nI,nI);
TDS.manifoldType = zeros(nI,nI);
% Tropical points
TDS.points = zeros(2,nI);
TDS.pointIndices = zeros(3,nI);
TDS.pointMap = zeros(nI,nI,2);
% Tropical nullcline points
TDS.nullclinePoints = zeros(2,nI);
TDS.nullclinePointIndices = zeros(4,nI);
TDS.nullclinePointMap = cell(nI,nI);
TDS.distanceToNullclinePoint = zeros(nI,nI,nI);
TDS.isEquilibrium = zeros(1,nI);


% Auxiliary matrix
A = [-TDS.deg, ones(nI,1)];


%%% Create flow vectors
for i = 1:nI
    if i <= nU
        TDS.flowVectors(:,i) = [TDS.delta(i); 0];
    else
        TDS.flowVectors(:,i) = [0; TDS.delta(i)];
    end
end
for i = 1:nU
    for j = 1:nV
        J = j + nU;
        if isequal(F(i,2:4), G(j,2:4))
            TDS.equalDegreeIndex(i) = J;
            TDS.equalDegreeIndex(J) = i;
            TDS.flowVectors(:,i) = sum(TDS.flowVectors(:,[i j]),2);
            TDS.flowVectors(:,J) = TDS.flowVectors(:,i);
        end
    end
end


%%% Find tropical points
for i = 1:nI-2
    for j = i+1:nI-1
        for k = j+1:nI
            I = [i,j,k];
            if abs(det(A(I,:))) < tol
                continue
            end
            p = A(I,:) \ TDS.alpha(I);
            zmax = max(TDS.deg * p(1:2) + TDS.alpha);
            if nI == 3 || p(3)+tol > zmax
                nP = nP + 1;
                TDS.pointIndices(:,nP) = I;
                TDS.points(:,nP) = p(1:2);
                for m = 1:3
                    iP = I(mod(m,3)+1);
                    jP = I(mod(m+1,3)+1);                    
                    if TDS.numEndpoints(iP,jP) == 1 && all(abs(p(1:2) - TDS.points(:,TDS.pointMap(iP,jP,1))) < tol)
                        % The endpoint has already been added for this manifold.
                        continue
                    elseif TDS.numEndpoints(iP,jP) == 2
                        % A manifold cannot have more than 2 endpoints
                        continue
                    end
                    TDS.numEndpoints(iP,jP) = TDS.numEndpoints(iP,jP) + 1;
                    TDS.numEndpoints(jP,iP) = TDS.numEndpoints(iP,jP);
                    TDS.pointMap(iP,jP,TDS.numEndpoints(iP,jP)) = nP;
                    TDS.pointMap(jP,iP,TDS.numEndpoints(iP,jP)) = nP;
                end
            end
        end
    end
end
TDS.points = TDS.points(:,1:nP);



%%% Tropical manifolds
for i = 1:nI-1
    for j = i+1:nI
        if nP == 0
            pi = [TDS.deg(i,:), TDS.alpha(i)];
            pj = [TDS.deg(j,:), TDS.alpha(j)];
            w_ij = pj-pi;
            for k = 1:nI
                if k == i || k == j
                    continue
                end
                pk = [TDS.deg(k,:), TDS.alpha(k)];
                w_ik = pk-pi;
                ratio = norm(w_ij)/norm(w_ik);
                p = pi + ratio*w_ij;
                if dot(w_ij(1:2), w_ik(1:2)) > 0 && ratio > 1 && TDS.alpha(k) > p(3)
                    k = k-1;
                    break
                end
            end
            if k == nI
                nM = nM + 1;
                TDS.manifoldDirection(:,nM) = [-w_ij(2); w_ij(1)];
                TDS.manifoldMap(i,j) = nM;
                TDS.manifoldMap(j,i) = nM;
            end
        else
            iP = TDS.pointMap(i,j,1);
            if ~iP
                continue
            end
            I = TDS.pointIndices(:,iP);
            for m = 1:3
                nM = nM + 1;
                iV = I(mod(m,3)+1);
                jV = I(mod(m+1,3)+1);
                r = det(A(I,:)) * cross(A(iV,:), A(jV,:));
                TDS.manifoldDirection(:,nM) = 1/norm(r(1:2)) * r(1:2);
                TDS.manifoldMap(iV,jV) = nM;
                TDS.manifoldMap(jV,iV) = nM;
            end
        end
        
        %%% Determine manifold type
        d_i = TDS.flowVectors(:,i);
        d_j = TDS.flowVectors(:,j);
        n_ij = [0 -1; 1 0]*TDS.manifoldDirection(:,TDS.manifoldMap(i,j));
        classifier = dot(d_i,n_ij)*dot(d_j,n_ij);
        if classifier > tol
            TDS.manifoldType(i,j) = TDS.CROSSING;
        elseif abs(det([d_i'; d_j'])) > tol
            if abs(classifier) < tol
                TDS.manifoldType(i,j) = TDS.TANGENTIAL_FILIPPOV;
            else
                TDS.manifoldType(i,j) = TDS.TRANSVERSAL_FILIPPOV;
            end
        else
            if abs(classifier) < tol
                TDS.manifoldType(i,j) = TDS.TANGENTIAL_NULLCLINE;
            else
                TDS.manifoldType(i,j) = TDS.TRANSVERSAL_NULLCLINE;
            end
            if i <= nU
                k_start = nU+1;
                k_end = nI;
            else
                k_start = 1;
                k_end = nU;
            end
            % Find nullcline points
            distances = [];
            for k = k_start:k_end-1
                for l = k_start+1:k_end
                    B = [diff(TDS.deg([i,j],:)); diff(TDS.deg([k,l],:))];
                    b = [diff(TDS.alpha([j,i])); diff(TDS.alpha([l,k]))];
                    if abs(det(B)) < tol
                        continue
                    end
                    p = B \ b;
                    if TDS.numEndpoints(i,j) == 2
                        Q = TDS.points(:,TDS.pointMap(i,j,:));
                        qmin = min(Q,[],2);
                        qmax = max(Q,[],2);
                        pqmin = p - qmin;
                        pqmax = p - qmax;
                        if dot(pqmin,pqmax) > tol
                            % Point is outside the endpoints
                            continue
                        end
                    elseif TDS.numEndpoints(i,j) == 1
                        q = TDS.points(:,TDS.pointMap(i,j,1));
                        rp = p-q;
                        rq = TDS.manifoldDirection(:,TDS.manifoldMap(i,j));
                        if dot(rp,rq) < -tol
                            % Point is outside the half line
                            continue
                        end
                    else
                        continue
                    end
                    TmaxK = max(TDS.deg(k_start:k_end,:) * p + TDS.alpha(k_start:k_end));
                    TmaxI = max(TDS.deg * p + TDS.alpha);
                    pk3 = TDS.deg(k,:) * p + TDS.alpha(k);
                    pi3 = TDS.deg(i,:) * p + TDS.alpha(i);
                    if length(k_start:k_end) == 2 || ...
                            (pk3 + tol >= TmaxK && pi3 + tol >= TmaxI)
                        nNP = nNP + 1;
                        TDS.nullclinePoints(:,nNP) = p;
                        TDS.nullclinePointIndices(:,nNP) = [i j k l];
                        TDS.nullclinePointMap{i,j} = [TDS.nullclinePointMap{i,j} nNP];
                        TDS.isEquilibrium(nNP) = TDS.delta(k) * TDS.delta(l) < 0;
                        distances = [distances, norm(p - TDS.points(:,TDS.pointMap(i,j,1)))];
                    end
                end
            end
            [~,id] = sort(distances,'ComparisonMethod','abs');
            TDS.nullclinePointMap{i,j} = TDS.nullclinePointMap{i,j}(id);
        end
    end
end
TDS.nullclinePoints = TDS.nullclinePoints(:,1:nNP);
TDS.nullclinePointIndices = TDS.nullclinePointIndices(:,1:nNP);
TDS.manifoldDirection = TDS.manifoldDirection(:,1:nM);
end

