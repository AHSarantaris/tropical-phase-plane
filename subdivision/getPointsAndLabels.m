function [P,L] = getPointsAndLabels(F,G)
    Pu = F(:,[3 4 2]);
    Pv = G(:,[3 4 2]);
    P = [Pu;Pv];
    dF = [F(:,1), zeros(size(F,1), 1)];
    dG = [zeros(size(G,1), 1), G(:,1)];
    L = [dF;dG]; 
end