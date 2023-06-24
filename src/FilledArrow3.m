function varargout = FilledArrow3(p,dD,dS,ulim,vlim,cD)
colors = [0.6 0.7 1; 0.8 0.9 1; 1 0.5 0.5; 1 0.75 0.6];
ax = gca;
ax.Units = "pixels";
pos = ax.Position;
ax.Units = "normalized";
wu = diff(ulim);
hu = diff(vlim);
wp = pos(3);
hp = pos(4);
pMu = [wp/wu;hp/hu];
uMp = 1./pMu;
dD = dD/norm(dD);
dS = dS/norm(dS);
dD = pMu.*dD;
dS = pMu.*dS;
dD = 8/norm(dD)*dD;
dS = 8/norm(dS)*dS;
dD = uMp.*dD;
dS = uMp.*dS;
nD = pMu.*dD;
nS = pMu.*dS;
nD = [-nD(2);nD(1)];
nS = [-nS(2);nS(1)];
nD = 0.5*uMp.*nD;
nS = 0.5*uMp.*nS;
% dD1 = 1/3*dD;
% dS1 = 1/3*dS;
% dD2 = 2/3*dD;
% dS2 = 2/3*dS;
dD1 = 0;
dD2 = 0.8*dD;
p2 = p(1:2);

nt = 0.2*nD;
dt = 0.7*dD;
qD = zeros(2,7);
qD(:,1) = p2-dD1+nD;
qD(:,2) = p2+dD2;
qD(:,3) = p2-dD1-nD;
qD(:,4) = p2-dD1-nt;
qD(:,5) = p2-dt-nt;
qD(:,6) = p2-dt+nt;
qD(:,7) = p2-dD1+nt;

% s = 1.4;
% qS = zeros(2,4);
% qS(:,1) = p2-dS1+s*nS;
% qS(:,2) = p2+s*dS2+0.1*nS;
% qS(:,3) = p2+s*dS2-0.1*nS;
% qS(:,4) = p2-dS1-s*nS;

% q = zeros(2,3);
% q(:,1) = p+dh;
% q(:,2) = p-dh-n;
% q(:,3) = p-dh+n;
qs = zeros(2,4);
qs(:,1) = min(qD,[],2)-1.2*norm(nD);
qs(:,3) = max(qD,[],2)+1.2*norm(nD);
qs(:,2) = [qs(1,1); qs(2,3)];
qs(:,4) = [qs(1,3); qs(2,1)];

qS = zeros(2,3);
s = 0.8;
if dS(2) == 0
    if dS(1) > 0
        qS(:,1) = qs(:,1);
        qS(:,2) = qs(:,2);
        qS(:,3) = qs(:,1) + 0.5*(qs(:,2) - qs(:,1)) + s*(qs(:,4)-qs(:,1)).*dS;
    elseif dS(1) < 0
        qS(:,1) = qs(:,3);
        qS(:,2) = qs(:,4);
        qS(:,3) = qs(:,4) + 0.5*(qs(:,3) - qs(:,4)) + s*(qs(:,4)-qs(:,1)).*dS;
    end
else
    if dS(2) > 0
        qS(:,1) = qs(:,1);
        qS(:,2) = qs(:,4);
        qS(:,3) = qs(:,1) + 0.5*(qs(:,4) - qs(:,1)) + s*(qs(:,2)-qs(:,1)).*dS;
    elseif dS(2) < 0
        qS(:,1) = qs(:,2);
        qS(:,2) = qs(:,3);
        qS(:,3) = qs(:,2) + 0.5*(qs(:,3) - qs(:,2)) + s*(qs(:,2)-qs(:,1)).*dS;
    end
end

% hs = patch(qs(1,:),qs(2,:), p(3)*ones(1,4)+1e-3,'k','FaceColor','white','LineWidth',0.5);
% uistack(hs,"top")
% if any(dS)
%     hS = patch(qS(1,:),qS(2,:), p(3)*ones(1,size(qS,2))+2e-3,'k','FaceColor','k','LineWidth',0.5);
%     uistack(hS,"top")
% end
plot3(p(1),p(2),p(3), 's','MarkerFaceColor','w','MarkerSize',20,'LineWidth',1,'Color',0.2*[1 1 1])
h = patch(qD(1,:),qD(2,:), p(3)*ones(1,size(qD,2))+3e-3,'k','FaceColor',cD,'LineWidth',1);
uistack(h,"top")
if nargout
    varargout{1} = h;
end
% plot(p(1),p(2),'k.')
end