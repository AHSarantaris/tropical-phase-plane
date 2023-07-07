function U = uniquePoints(P,tol)
    U = zeros(size(P));
    nP = size(P,2);
    nU = 0;
    for i = 1:nP
        dP = abs(P(:,i) - P(:,i+1:end));
        if sum(all(dP < tol)) == 0
            nU = nU + 1;
            U(:,nU) = P(:,i);
        end
    end
    U = U(:,1:nU);
end