function Q = pointsInside(P, ulim, vlim)
Q = [];
p1 = P(:,1);
p2 = P(:,2);
pos1 = position(p1);
pos2 = position(p2);
if pos1 == 5 && pos2 == 5
    Q = [p1 p2];
elseif pos1 == 5
    Q(:,1) = p1;
    Q(:,2) = intersectionFromInside(p1,p2,pos2);
elseif pos2 == 5
    Q(:,1) = intersectionFromInside(p2,p1,pos1);
    Q(:,2) = p2;
else
    qu1 = intersection(1,p1,p2,ulim(1),vlim);
    qu2 = intersection(1,p1,p2,ulim(2),vlim);
    qv1 = intersection(2,p1,p2,vlim(1),ulim);
    qv2 = intersection(2,p1,p2,vlim(2),ulim);
    Q4 = {qu1 qu2 qv1 qv2};
    for i = 1:3
        for j = i+1:4
            if ~isempty(Q4{i}) && ~isempty(Q4{j})
                Q = cell2mat(Q4([i j]));
            end                
        end
    end
end

    function q = intersection(ip,p1,p2,line,lim)
        r = p2 - p1;
        if r(ip) == 0
            q = [];
            return
        end
        tp = norm(r);
        tq = (line - p1(ip)) / r(ip);
        q = p1 + tq * r;
        iq = mod(ip,2) + 1;
%         if tq < 0 || tp < tq || q(iq) < lim(1) || lim(2) < q(iq)
%             q = [];
%         end
    end

    function q = intersectionFromInside(pi,po,pos)
        if pos == 4
            q = intersection(1,pi,po,ulim(1),vlim);
        elseif pos == 6
            q = intersection(1,pi,po,ulim(2),vlim);
        elseif pos == 2
            q = intersection(2,pi,po,vlim(1),ulim);
        elseif pos == 8
            q = intersection(2,pi,po,vlim(2),ulim);
        elseif pos == 1
            q1 = intersection(1,pi,po,ulim(1),vlim);
            q2 = intersection(2,pi,po,vlim(1),ulim);
            if q1(2) > vlim(1)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 3
            q1 = intersection(1,pi,po,ulim(2),vlim);
            q2 = intersection(2,pi,po,vlim(1),ulim);
            if q1(2) > vlim(1)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 7
            q1 = intersection(1,pi,po,ulim(1),vlim);
            q2 = intersection(2,pi,po,vlim(2),ulim);
            if q1(2) < vlim(2)
                q = q1;
            else
                q = q2;
            end
        elseif pos == 9
            q1 = intersection(1,pi,po,ulim(2),vlim);
            q2 = intersection(2,pi,po,vlim(2),ulim);
            if q1(2) < vlim(2)
                q = q1;
            else
                q = q2;
            end
        else
            error('Third argument must be between 1 and 9 but not 5')
        end
    end

    function pos = position(p)
        if p(2) < vlim(1)
            if p(1) < ulim(1)
                pos = 1;
            elseif p(1) > ulim(2)
                pos = 3;
            else
                pos = 2;
            end
        elseif p(2) > vlim(2)
            if p(1) < ulim(1)
                pos = 7;
            elseif p(1) > ulim(2)
                pos = 9;
            else
                pos = 8;
            end
        else
            if p(1) < ulim(1)
                pos = 4;
            elseif p(1) > ulim(2)
                pos = 6;
            else
                pos = 5;
            end  
        end
    end

end
