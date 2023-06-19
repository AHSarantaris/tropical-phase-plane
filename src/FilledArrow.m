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
t = 0.9;
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

%     n = 1.2*n;
%     nt = 0.05*n;
%     dt = 2*d;
%     d1 = 0.6*d;
%     s = 2;
%     d2 = 1/6*d;
%     q = zeros(2,0);
%     q(:,1) = p-d1+n;
%     q(:,2) = p+d2;
%     q(:,3) = p-d1-n;
%     q(:,4) = p + t*(q(:,3)-q(:,2));
%     q(:,5) = p-nt;
%     p2 = p - d1;
%     q(:,6) = p2+d2-nt;
%     q(:,7) = p2-d1-n;
%     q(:,8) = p2 + t*(q(:,7)-q(:,6));
%     q(:,9) = p2-nt;
%     q(:,10) = p-dt-nt;
%     q(:,11) = p-dt+nt;
%     q(:,12) = p2+nt;
%     q(:,13) = p2 + t*(q(:,1)-q(:,2));
%     q(:,14) = p2-d1+n;
%     q(:,15) = p2+d2+nt;
%     q(:,16) = p+nt;
%     q(:,17) = p + t*(q(:,1)-q(:,2));

else
    s = 1.4;
    q = zeros(2,4);
    q(:,1) = p-d1+s*n;
    q(:,2) = p+s*d2+0.1*n;
    q(:,3) = p+s*d2-0.1*n;
    q(:,4) = p-d1-s*n;
    
%     d1 = 0.5*d;
%     s = 2;
%     d2 = 1/3*d;
%     q = zeros(2,7);
%     q(:,1) = p-d1+s*n;
%     q(:,2) = p+s*d2+0.1*n;
%     q(:,3) = p+s*d2-0.1*n;
%     q(:,4) = p-d1-s*n;
%     q(:,5) = p + t*(q(:,4)-q(:,3));
%     q(:,6) = p;
%     q(:,7) = p + t*(q(:,1)-q(:,2));
%     
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