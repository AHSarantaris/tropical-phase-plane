function varargout = FilledArrow2(p,d,ulim,vlim,c,lineWidth,hasTail,isDoubleHeaded)
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
d = 0.05*(lineWidth+15)*d;
d = uMp.*d;


n = pMu.*d;
n = [-n(2);n(1)];
n = 0.5*uMp.*n;


d1 = 1/3*d;
d2 = 2/3*d;


n = 1.8*n;
v1 = -d1+n;
v2 = -d1-n;
q(:,1) = p+v1;
q(:,2) = p+d2;
q(:,3) = p+v2;

h = plotArrow(q);

if isDoubleHeaded
    p = p-1.1*d2;
    q(:,1) = p+v1;
    q(:,2) = p+d2;
    q(:,3) = p+v2;
    h = plotArrow(q);
end

function h = plotArrow(q)
    h=plot(q(1,:),q(2,:),'Color',c,'LineWidth',lineWidth,'LineJoin','miter');
    hold on
    if hasTail
        q2 = [q(:,2) q(:,2)-5*d1];
        h = plot(q2(1,:),q2(2,:),'Color',c,'LineWidth',lineWidth,'LineJoin','miter');
    end

end

% 
% % q = zeros(2,3);
% % q(:,1) = p+dh;
% % q(:,2) = p-dh-n;
% % q(:,3) = p-dh+n;
% h = patch(q(1,:),q(2,:),'k','FaceColor',c,'LineWidth',0.5);
if nargout
    varargout{1} = h;
end
% plot(p(1),p(2),'k.')
end

