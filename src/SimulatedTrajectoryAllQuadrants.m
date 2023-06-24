function SimulatedTrajectoryAllQuadrants(n,iq0,axs,F1,G1,w0,epsilon,tspan,maxseconds)

if ~exist('epsilon','var')
    epsilon = 0.2;
end
if ~exist('tspan','var')
    tspan = [0 1e4];
end
if ~exist('maxseconds','var')
    maxseconds = 1;
end


[F,G] = tropicalPolynomialsAllQuadrants(F1,G1);


iq = iq0;
w1 = w0;
W = zeros(n,2);
iw = 0;
for i = 1:n
    ax = axs{iq0};
    axes(ax)
    [~,w1] = SimulatedTrajectory(F{iq0},G{iq0},w1,epsilon,tspan,maxseconds);
    iq1 = nextQuadrant(iq0,w1,ax.XLim,ax.YLim);
    if iq1 == 0 || norm(w1 - w0) < 1e-6
        break
    end
    if (iq0 == 1 && iq1 == 2) || (iq0 == 4 && iq1 == 3)
        w1(1) = axs{iq1}.XLim(2);
    elseif (iq0 == 2 && iq1 == 1) || (iq0 == 3 && iq1 == 4)
        w1(1) = axs{iq1}.XLim(1);
    elseif (iq0 == 1 && iq1 == 4) || (iq0 == 2 && iq1 == 3)
        w1(2) = axs{iq1}.YLim(2);
    elseif (iq0 == 4 && iq1 == 1) || (iq0 == 3 && iq1 == 2)
        w1(2) = axs{iq1}.YLim(1);
    end
    if iq1 == iq
        iw = iw + 1;
        W(iw,:) = w1;
    end
    iq0 = iq1;
end
W = W(1:iw,:)


end

function iq1 = nextQuadrant(iq0,w1,ulim,vlim)
wmin = zeros(1,2);
[wmin(1),iu] = min(abs(ulim - w1(1)));
[wmin(2),iv] = min(abs(vlim - w1(2)));
[~,iw] = min(wmin);
iq1 = 0;
if iq0 == 1
    if iw == 1
        if iu == 1
            iq1 = 2;
        end
    else
        if iv == 1
            iq1 = 4;
        end
    end
elseif iq0 == 2
    if iw == 1
        if iu == 2
            iq1 = 1;
        end
    else
        if iv == 1
            iq1 = 3;
        end
    end
elseif iq0 == 3
    if iw == 1
        if iu == 2
            iq1 = 4;
        end
    else
        if iv == 2
            iq1 = 2;
        end
    end
elseif iq0 == 4
    if iw == 1
        if iu == 1
            iq1 = 3;
        end
    else
        if iv == 2
            iq1 = 1;
        end
    end            
end
end

