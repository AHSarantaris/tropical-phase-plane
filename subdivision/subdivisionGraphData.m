function subdivision = subdivisionGraphData(P)

tol = 1e-6;
nP = size(P,1);

%%% Triangulation of convex hull
warning('off')
isplanar = false;
E2 = zeros(0,2);
try
    k = convhull(P,"Simplify",false);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
    k = k(N(:,3)>0,:);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
    %%% Find edges
    E1 = edges(T);
    nE = 0;
    for i = 1:size(E1,1)
        iT = edgeAttachments(T,E1(i,1),E1(i,2));
        iT = iT{:};
        if length(iT) < 2 || norm(N(iT(1),:) - N(iT(2),:)) > tol
            nE = nE+1;
            E2(nE,:) = E1(i,:);

        end
    end
catch
    isplanar = true;
    k = convhull(P(:,1:2),"Simplify",false);
    E2 = [k k([2:length(k), 1])];
end
warning('on')

ku = unique(k(:));

subdivision.isplanar = isplanar;
subdivision.faces = k;
subdivision.vertices = ku;
subdivision.edges = E2;

end