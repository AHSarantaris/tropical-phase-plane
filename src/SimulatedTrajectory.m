function [h,w1] = SimulatedTrajectory(F,G,w0,epsilon,tspan,maxseconds,color)

% Initialize variables
lineWidth = 2;
ax = gca;
hold(ax,'on')
nW = size(w0,1);

if ~exist('epsilon','var')
    epsilon = 0.2;
end
if ~exist('tspan','var')
    tspan = [0 1e30];
end
if ~exist('maxseconds','var')
    maxseconds = 4;
end
if ~exist('color','var')
    color = [0 0.6 0; 0.6 0 0.6; 0.7 0.5 0; 0 0.5 0.7];
end

%%% ODE system
% udot = @(w) F(:,1)' * exp(1/epsilon * F(:,2:end) * [1;w]);
% vdot = @(w) G(:,1)' * exp(1/epsilon * G(:,2:end) * [1;w]);
udot = @(w) epsilon * F(:,1)' * exp(1/epsilon * F(:,2:end) * [1;w]);
vdot = @(w) epsilon * G(:,1)' * exp(1/epsilon * G(:,2:end) * [1;w]);
wdot = @(t,w) [udot(w); vdot(w)];

% wdot = @(t,w) TropicalSystem(F,G,w,eps);

%%% Jacobian
J11 = @(w) F(:,3)' .* F(:,1)' * exp(1/epsilon * F(:,2:end) * [1;w]);
J12 = @(w) F(:,4)' .* F(:,1)' * exp(1/epsilon * F(:,2:end) * [1;w]);
J21 = @(w) G(:,3)' .* G(:,1)' * exp(1/epsilon * G(:,2:end) * [1;w]);
J22 = @(w) G(:,4)' .* G(:,1)' * exp(1/epsilon * G(:,2:end) * [1;w]);
J = @(t,w) [J11(w), J12(w); J21(w), J22(w)];

for i = 1:nW
    % Solve and plot
    initialTime = cputime;
%     options = odeset('Jacobian',J,...
%         'Events',@(t,w)event(w,ax.XLim,ax.YLim,initialTime,maxseconds));
    options = odeset('Events',@(t,w)event(w,ax.XLim,ax.YLim,initialTime,maxseconds),'RelTol',1e-9,'AbsTol',1e-9);
    [~,wID1] = lastwarn;
    warning('off','all')
%     [~,w] = ode23s(wdot,tspan,w0(i,:),options);
%     [~,w] = ode23tb(wdot,tspan,w0(i,:),options);
%     [~,w] = ode15s(wdot,tspan,w0(i,:),options);
    [~,w] = ode45(wdot,tspan,w0(i,:),options);
    warning('on','all')
    [~,wID2] = lastwarn;
    if ~isequal(wID1, wID2)
        warning(lastwarn)
    end
    handle = plot(w(:,1),w(:,2),'Color',color(i,:),'LineWidth',(nW-i)+lineWidth,'LineStyle','-');
    plot(w0(i,1),w0(i,2),'ok','MarkerSize',6,'MarkerFaceColor',color(i,:),'LineWidth',1)

end
if nargout >= 1
    h = handle;
end
if nargout >= 2
    w1 = w(end,:); 
end
hold(ax,'off')
end

function [value, isterminal, direction] = event(w,ulim,vlim,initialTime,maximumTime)
dt = cputime - initialTime;
timeleft = maximumTime - dt;
du = (ulim - w(1)) .* [-1 1];
dv = (vlim - w(2)) .* [-1 1];
wmin = zeros(1,2);
wmin(1) = min(du);
wmin(2) = min(dv);
value = min(timeleft,min(wmin));
isterminal = 1;
direction = 0;
end

