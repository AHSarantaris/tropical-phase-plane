function [hT,axs] = TropicalPhaseQuadrants(F1,G1,ulim,vlim)
f = gcf;
% f.Units = 'normalized';
w = 0.5;
h = 0.8;
% f.OuterPosition = [(1-w)/2 (1-h)/2 w h];´
f.Units = 'pixels';

[F,G] = TropicalQuadrants(F1,G1);

ulims = zeros(4,2);
vlims = zeros(4,2);
if ~exist('ulim','var') || ~exist('vlim','var')
    for i = 1:4
        TL = TropicalCurves(F{i},G{i});
        Q = [TL.points,TL.tropicalPoints];
        nQ = size(Q,2);
        if nQ > 1
            qlim = [min(Q,[],2), max(Q,[],2)];
            midU = sum(qlim(1,:))/2;
            midV = sum(qlim(2,:))/2;
            rangeU = diff(qlim(1,:));
            rangeV = diff(qlim(2,:));
            if rangeU == 0 && rangeV == 0
                rangeU = 1;
                rangeV = 1;
            elseif rangeU == 0
                rangeU = rangeV;
            elseif rangeV == 0
                rangeV = rangeU;
            end
        elseif nQ == 1
            midU = Q(1,1);
            midV = Q(2,1);
            rangeU = 1;
            rangeV = 1;
        end
        r = 1.5;
        ulims(i,:) = midU + rangeU*[-r r];
        vlims(i,:) = midV + rangeV*[-r r];
    end
    
    for i = [[1;4] [2;3]]
        I = i';
        umax = [min(ulims(I,1)) max(ulims(I,2))];
        for j = I
            ulims(j,:) = umax;
        end
    end
    for i = [[1;2] [3;4]]
        I = i';
        vmax = [min(vlims(I,1)) max(vlims(I,2))];
        for j = I
            vlims(j,:) = vmax;
        end
    end

    udiff = diff(ulims,1,2);
    vdiff = diff(vlims,1,2);

    udiffmax = max(diff(ulims,1,2));
    vdiffmax = max(diff(vlims,1,2));

    for i = 1:4
        ulims(i,:) = udiffmax / udiff(i) * ulims(i,:);
        vlims(i,:) = vdiffmax / vdiff(i) * vlims(i,:);
    end
else
    for i = [1 4]
        ulims(i,:) = ulim;
    end
    for i = [2 3]
        ulims(i,:) = sort(-ulim);
    end
    for i = [1 2]
        vlims(i,:) = vlim;
    end
    for i = [3 4]
        vlims(i,:) = sort(-vlim);
    end
end



hf = tiledlayout(f,2,2,"Padding","compact","TileSpacing","compact");
ax = cell(1,4);
for i=[2 1 3 4]
    ax{i} = nexttile;
    TropicalPhasePlane(F{i},G{i},ulims(i,:),vlims(i,:))
%     axis(ax{i}, "equal")
end
ax{1}.XAxisLocation = 'top';
ax{2}.XAxisLocation = 'top';
ax{1}.YAxisLocation = 'right';
ax{4}.YAxisLocation = 'right';

if nargout >= 1
    hT = hf;
end
if nargout >= 2
    axs = ax;
end

end

