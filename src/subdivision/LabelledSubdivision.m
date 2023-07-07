function LabelledSubdivision(F,G,varargin)

[P,D] = getPointsAndLabels(F,G);

tol = 1e-6;
isplanar = rank(P,tol) == 2;
nP = size(P,1);

%%% Initialize axes
ax = gca;
hold(ax,'on')
grid on


%%% Triangulation of convex hull
warning('off')
if isplanar
    k = convhull(P(:,1:2),"Simplify",false);
    ku = k;
    N = [0,0,1];
else
    k = convhull(P,"Simplify",false);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
    k = k(N(:,3)>0,:);
    ku = unique(k);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
    warning('on')

    %%% Find edges
    E1 = edges(T);
    E2 = zeros(0,2);
    nE = 0;
    for i = 1:size(E1,1)
        iT = edgeAttachments(T,E1(i,1),E1(i,2));
        iT = iT{:};
        if length(iT) < 2 || norm(N(iT(1),:) - N(iT(2),:)) > tol
            nE = nE+1;
            E2(nE,:) = E1(i,:);
        end
    end
end


%%% Plot settings
hold on
%     D = max(P) - min(P);
%     pbaspect([D(1) D(2) min(D(1:2))])
%     pbaspect([D(1) D(2) max(D(3),1)])
nlim = [min(P(:,1)) max(P(:,1))];
mlim = [min(P(:,2)) max(P(:,2))];
xlim(nlim + [-1 1])
ylim(mlim + [-1 1])

%%% Plot faces
if isplanar
    patch(P(k,1),P(k,2),'white','LineStyle','none',varargin{:})
else
    patch('Faces',k,'Vertices',P(:,1:2),'FaceColor','white','LineStyle','none',varargin{:})
end

%%% Plot edges
if isplanar
    plot(P(k,1),P(k,2),'.-k','LineWidth',1)
else
    for i = 1:size(E2,1)
        pE = P(E2(i,:),:);
        plot(pE(:,1),pE(:,2),'-k','LineWidth',1)
    end
end

%%% Plot flow vectors
nD = 0;
dD = zeros(0,2);
dS = zeros(0,2);
for i = 1:nP
    hasSubdominant = false;
    imax = i;
    for j = 1:nP
        if i == j
            continue
        elseif P(i,1) == P(j,1) && P(i,2) == P(j,2)
            nD = nD + 1;
            if abs(P(i,3)-P(j,3)) < tol
                dD(nD,:) = D(i,:) + D(j,:);
                dS(nD,:) = [0 0];
            elseif P(i,3) > P(j,3)
                dD(nD,:) = D(i,:);
                dS(nD,:) = D(j,:);
            else % P(i,3) < P(j,3)
                dD(nD,:) = D(j,:);
                dS(nD,:) = D(i,:);
                imax = j;
            end
            hasSubdominant = true;
            break
        end
    end
    if ~hasSubdominant
        nD = nD + 1;
        dD(nD,:) = D(i,:);
        dS(nD,:) = [0 0];
    end
    markerEdgeColor = 'black';
    if ~isplanar
        attachments = vertexAttachments(T,i);
        attachments = attachments{:};
        if isempty(attachments)
            markerEdgeColor = 0.4*[1 1 1];
        end
    end
    plot(P(i,1),P(i,2),'s','MarkerSize',20,'MarkerFaceColor','white','MarkerEdgeColor',markerEdgeColor)
    arrow = FilledArrow(P(i,1:2)',dD(nD,:)',xlim,ylim,markerEdgeColor,true);
    arrow.EdgeColor = markerEdgeColor;
end
hold off
xticks(nlim(1):nlim(2))
yticks(mlim(1):mlim(2))

end