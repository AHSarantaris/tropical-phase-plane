function [ulim,vlim] = getLimits(points,ulim,vlim)

if ~isempty(ulim) && ~isempty(vlim)
    return
end
minmax = [min(points,[],2), max(points,[],2)];
nP = size(points,2);
if nP > 1
    midU = sum(minmax(1,:))/2;
    midV = sum(minmax(2,:))/2;
    rangeU = diff(minmax(1,:));
    rangeV = diff(minmax(2,:));
else
    if nP == 1
        midU = points(1,1);
        midV = points(2,1);
    else
        midU = 0;
        midV = 0;
    end
    rangeU = 1;
    rangeV = 1;
end
if rangeU == 0
    rangeU = rangeV;
end
if rangeV == 0
    rangeV = rangeU;
end
rangeMax = max(rangeU,rangeV);
if isempty(ulim)
    rangeU = rangeMax;
    ulim = midU + rangeU*[-1,1];
    ulim = [floor(ulim(1)*10)/10 ceil(ulim(2)*10)/10];
end
if isempty(vlim)
    rangeV = rangeMax;
    vlim = midV + rangeV*[-1,1];
    vlim = [floor(vlim(1)*10)/10 ceil(vlim(2)*10)/10];
end

end