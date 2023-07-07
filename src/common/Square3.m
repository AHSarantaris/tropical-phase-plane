function varargout = Square3(p,color)

r = 0.15;
qs = zeros(2,4);
qs(:,1) = p(1:2)+r*[-1;1];
qs(:,2) = p(1:2)+r*[1;1];
qs(:,3) = p(1:2)+r*[1;-1];
qs(:,4) = p(1:2)+r*[-1;-1];

hs = patch(qs(1,:),qs(2,:), p(3)*ones(1,4)+1e-3,'k','FaceColor',color,'LineWidth',0.5,'FaceLighting','none');
uistack(hs,"top")
if nargout
    varargout{1} = hs;
end
end