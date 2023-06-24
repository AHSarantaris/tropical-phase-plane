function [F,G] = tropicalPolynomialsAllQuadrants2(F1,G1)
F = {F1 F1 F1 F1};
G = {G1 G1 G1 G1};

nF = F1(:,3)';
nG = G1(:,3)';
mF = F1(:,4)';
mG = G1(:,4)';

F{2}(:,1) = (-1).^nF .* F1(:,1)';
G{2}(:,1) = (-1).^nG .* G1(:,1)';

F{4}(:,1) = (-1).^mF .* F1(:,1)';
G{4}(:,1) = (-1).^mG .* G1(:,1)';

F{3}(:,1) = (-1).^(nF + mF) .* F1(:,1)';
G{3}(:,1) = (-1).^(nG + mG) .* G1(:,1)';
end