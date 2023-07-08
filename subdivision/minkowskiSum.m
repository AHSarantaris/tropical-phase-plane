function M = minkowskiSum(P,Q)
nP = size(P,1);
nQ = size(Q,1);
M = zeros(nP*nQ,size(P,2));
nM = 0;
for i = 1:nP
    for j = 1:nQ
        nM = nM + 1;
        M(nM,:) = P(i,:) + Q(j,:);
    end
end
end