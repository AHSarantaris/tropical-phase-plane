function wdot = tropicalODE(~,w,F,G,eps)
Fw = F(:,2:4)*[1;w];
Gw = G(:,2:4)*[1;w];
FGmax = max([Fw;Gw]);
epsFw = exp((Fw-FGmax)./eps);
epsGw = exp((Gw-FGmax)./eps);
epsFG = sum([epsFw;epsGw]);
udot = 1/epsFG*F(:,1)'*epsFw;
vdot = 1/epsFG*G(:,1)'*epsGw;
wdot = [udot;vdot];
end