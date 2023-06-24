function p2 = boundaryPoint(p1, r, ulim, vlim)
if (p1(1) < ulim(1) && r(1) < 0) ...
        || (p1(1) > ulim(2) && r(1) > 0) ...
        || (p1(2) < vlim(1) && r(2) < 0) ...
        || (p1(2) > vlim(2) && r(2) > 0)
    p2 = p1;
    return
end
if r(2) == 0
    if r(1) < 0
        p2 = [ulim(1); p1(2)];
    else
        p2 = [ulim(2); p1(2)];
    end
    return
elseif r(1) == 0
    if r(2) < 0
        p2 = [p1(1); vlim(1)];
    else
        p2 = [p1(1); vlim(2)];
    end
    return
end
if r(1) < 0
    t1 = (ulim(1) - p1(1)) / r(1);
else
    t1 = (ulim(2) - p1(1)) / r(1);
end
if r(2) < 0
    t2 = (vlim(1) - p1(2)) / r(2);
else
    t2 = (vlim(2) - p1(2)) / r(2);
end
[tmin,imin] = min([t1,t2]);
p2 = p1 + tmin * r;
end