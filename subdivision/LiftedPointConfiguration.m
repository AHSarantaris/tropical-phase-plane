function LiftedPointConfiguration(P,h0,varargin)

tol = 1e-6;
nP = size(P,1);
if ~exist('h0','var') || isempty(h0)
    h0 = min(P(:,3))-1;
end

%%% Plot settings
liftedPointSize = 25;
defaultColor = [0.5 0.7 1];
grey = 0.4*[1 1 1];
nlim = [min(P(:,1)) max(P(:,1))];
mlim = [min(P(:,2)) max(P(:,2))];
xlim(nlim + [-1 1])
ylim(mlim + [-1 1])


%%% Initialize axes
ax = gca;
hold(ax,'on')

%%% Triangulation of convex hull
warning('off')
isplanar = false;
try
    k = convhull(P,"Simplify",false);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
    k = k(N(:,3)>0,:);
    T = triangulation(k,P(:,1),P(:,2),P(:,3));
    N = faceNormal(T);
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
catch
    k = convhull(P(:,1:2),"Simplify",false);
    isplanar = true;
end
warning('on')

ku = unique(k(:));
notk = setdiff(1:nP,ku);


%%% Plot lifted faces
if isplanar
    patch(P(k,1),P(k,2),P(k,3),defaultColor,'LineStyle','none','AmbientStrength',0.9,'SpecularStrength',0,'DiffuseStrength',0.5,varargin{:})
else
    patch('Faces',k,'Vertices',P,'LineStyle','none','FaceColor',defaultColor,'AmbientStrength',0.9,'SpecularStrength',0,'DiffuseStrength',0.5,varargin{:})
end

%%% Plot lifted edges
if isplanar
    plot3(P(k,1),P(k,2),P(k,3),'.-k','LineWidth',1)
else
    for i = 1:size(E2,1)
        pE = P(E2(i,:),:);
        plot3(pE(:,1),pE(:,2),pE(:,3),'-k','LineWidth',1)
    end
end

%%% Plot lifted points
plot3(P(:,1),P(:,2),P(:,3),'.k','MarkerSize',liftedPointSize)

%%% Plot line between lifted and projected points
for i = 1:nP
    line([P(i,1) P(i,1)],[P(i,2) P(i,2)],[h0 P(i,3)],'Color',grey)
end


%%% Plot projected faces
if isplanar
    patch(P(k,1),P(k,2),h0*ones(size(k)),defaultColor,'LineStyle','none','FaceLighting','none',varargin{:})
else
    patch('Faces',k,'Vertices',[P(:,1:2) h0*ones(nP,1)],'FaceColor',defaultColor,'LineStyle','none','FaceLighting','none',varargin{:})
end

%%% Plot projected edges
if isplanar
    plot3(P(k,1),P(k,2),h0*ones(size(k)),'.-k','LineWidth',1)
else
    for i = 1:size(E2,1)
        pE = P(E2(i,:),:);
        plot3(pE(:,1),pE(:,2),h0*[1 1],'-k','LineWidth',1)
    end
end

%%% Plot projected points
for i = 1:nP
    squareColor = 'black';
    if ~isplanar
        attachments = vertexAttachments(T,i);
        attachments = attachments{:};
        if isempty(attachments)
            squareColor = 'white';
        end
    end
    Square3(P(i,1:2)',h0,squareColor);
end

light('Position',[0 0 1])
view([-80 25])
grid on
hold off
xticks(nlim(1):nlim(2))
yticks(mlim(1):mlim(2))

end


%% Plot flat square in 3D
function Square3(p,z,color)
r = 0.15;
q = zeros(2,4);
q(:,1) = p + r*[-1;1];
q(:,2) = p + r*[1;1];
q(:,3) = p + r*[1;-1];
q(:,4) = p + r*[-1;-1];
h = patch(q(1,:),q(2,:), z*ones(1,4)+1e-3,'k','FaceColor',color,'LineWidth',0.5,'FaceLighting','none');
end