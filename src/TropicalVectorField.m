function [udot,vdot] = TropicalVectorField(F,G,epsilon)
udot = @(w) epsilon * F(:,1)' * exp(1/epsilon * F(:,2:end) * [1;w]);
vdot = @(w) epsilon * G(:,1)' * exp(1/epsilon * G(:,2:end) * [1;w]);
% wdot = @(t,w) [udot(w); vdot(w)];
end