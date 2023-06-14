function [F1] = TranslateTropicalPolynomial(du,dv,F1)
F1(:,2) = F1(:,2) - du*F1(:,3) - dv*F1(:,4) - max(- du*F1(:,3) - dv*F1(:,4)) - 1;

end