function h = TropicalPhasePlaneAllQuadrants(F1,G1,ulim,vlim)

[F,G] = tropicalPolynomialsAllQuadrants(F1,G1);

ulims = zeros(4,2);
vlims = zeros(4,2);
if ~exist('ulim','var') || ~exist('vlim','var')
    if ~exist('ulim','var')
        ulim = [];
    end
    if ~exist('vlim','var')
        vlim = [];
    end
    for i = 1:4
        TDS = TropicalDynamicalSystem(F{i},G{i});
        points = [TDS.points,TDS.nullclinePoints];
        [ulim1,vlim1] = getLimits(points,ulim,vlim);
        ulims(i,:) = ulim1;
        vlims(i,:) = vlim1;
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



hf = tiledlayout(2,2,"Padding","compact","TileSpacing","compact");
ax = cell(1,4);
for i=[2 1 3 4]
    ax{i} = nexttile;
    TropicalPhasePlane(F{i},G{i},ulims(i,:),vlims(i,:))
end
ax{1}.XAxisLocation = 'top';
ax{2}.XAxisLocation = 'top';
ax{1}.YAxisLocation = 'right';
ax{4}.YAxisLocation = 'right';

if nargout
    h = hf;
end

end

%% Get the tropical polynomials from each quadrant

function [F,G] = tropicalPolynomialsAllQuadrants(F1,G1)
F = {F1 F1 F1 F1};
G = {G1 G1 G1 G1};

nF = F1(:,3)';
nG = G1(:,3)';
mF = F1(:,4)';
mG = G1(:,4)';

F{2}(:,3) = -nF;
G{2}(:,3) = -nG;
F{2}(:,1) = (-1).^(nF + 1) .* F1(:,1)';
G{2}(:,1) = (-1).^nG .* G1(:,1)';

F{4}(:,4) = -mF;
G{4}(:,4) = -mG;
F{4}(:,1) = (-1).^mF .* F1(:,1)';
G{4}(:,1) = (-1).^(mG + 1) .* G1(:,1)';

F{3}(:,3) = -nF;
F{3}(:,4) = -mF;
G{3}(:,3) = -nG;
G{3}(:,4) = -mG;
F{3}(:,1) = (-1).^(nF + mF + 1) .* F1(:,1)';
G{3}(:,1) = (-1).^(nG + mG + 1) .* G1(:,1)';
end

