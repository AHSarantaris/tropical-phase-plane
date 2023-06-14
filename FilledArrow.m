function varargout = FilledArrow(p,d,ulim,vlim,c,hasTail)
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
d = pMu.*d;
d = 8/norm(d)*d;
d = uMp.*d;
n = pMu.*d;
n = [-n(2);n(1)];
n = 0.5*uMp.*n;
d1 = 1/3*d;
d2 = 2/3*d;
if hasTail
    nt = 0.2*n;
    dt = d;
    q = zeros(2,7);
    q(:,1) = p-d1+n;
    q(:,2) = p+d2;
    q(:,3) = p-d1-n;
    q(:,4) = p-d1-nt;
    q(:,5) = p-dt-nt;
    q(:,6) = p-dt+nt;
    q(:,7) = p-d1+nt;
else
    s = 1.4;
    q = zeros(2,4);
    q(:,1) = p-d1+s*n;
    q(:,2) = p+s*d2+0.1*n;
    q(:,3) = p+s*d2-0.1*n;
    q(:,4) = p-d1-s*n;
end

% q = zeros(2,3);
% q(:,1) = p+dh;
% q(:,2) = p-dh-n;
% q(:,3) = p-dh+n;
h = patch(q(1,:),q(2,:),'k','FaceColor',c,'LineWidth',0.5);
if nargout
    varargout{1} = h;
end
% plot(p(1),p(2),'k.')
end