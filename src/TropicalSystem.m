function wdot = TropicalSystem(F,G,w,eps)
Fw = F(:,2:4)*[1;w];
Gw = G(:,2:4)*[1;w];
Fmax = max(Fw);
Gmax = max(Gw);
Hmax = max(Fmax,Gmax);
udot = eps*F(:,1)'*exp(Fw/eps);
vdot = eps*G(:,1)'*exp(Gw/eps);
wdot = [udot;vdot];
end