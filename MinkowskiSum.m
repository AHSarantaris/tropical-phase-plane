function [M,dD,dS] = MinkowskiSum(P,Q,dP,dQ)
nP = size(P,1);
nQ = size(Q,1);
M = zeros(nP*nQ,size(P,2));
dD = zeros(size(M,1), 2);
dS = zeros(size(dD));
nM = 0;
for i = 1:nP
    for j = 1:nQ
        nM = nM + 1;
        M(nM,:) = P(i,:) + Q(j,:);
        dD(nM,:) = dP(i,:) + dQ(j,:);
%         if P(i,3) > Q(j,3)
%             dD(nM,:) = dP(i,:);
%             dS(nM,:) = dQ(j,:);
%         elseif P(i,3) < Q(j,3)
%             dD(nM,:) = dQ(j,:);
%             dS(nM,:) = dP(i,:);
%         else
% %             dD(nM,:) = dQ(j,:);
% %             dS(nM,:) = dP(i,:);
%             dD(nM,:) = dP(i,:) + dQ(j,:);
%         end
    end
end
end