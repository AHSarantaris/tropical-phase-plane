function SubdivisionHeightPlot(P,dP,showLabels)
    
    tol = 1e-6;
    isplanar = ~any(P(:,3)-P(1,3));    
    nP = size(P,1);

    

    %%% Initialize axes
    ax = gca;
    hold(ax,'on')
    ax.Units = 'pixels';
    axsize = ax.Position(3:4);
    ax.Units = 'normalized';
    axmean = mean(axsize);
    
    lineWidth = 0.003*axmean;
    markerSize = min(10,10*lineWidth);

    %%% Triangulation of convex hull
    warning('off')
    if isplanar
        k = convhull(P(:,1:2),"Simplify",true);
        ku = k;
        N = [0,0,1];
    else
        k = convhull(P,"Simplify",true);
        T = triangulation(k,P(:,1),P(:,2),P(:,3));
        N = faceNormal(T);
        k = k(N(:,3)>0,:);
        ku = unique(k);
        T = triangulation(k,P(:,1),P(:,2),P(:,3));
        N = faceNormal(T);
        warning('on')
    
        %%% Find edges
        E1 = edges(T);
        E2 = zeros(0,2);
        nE = 0;
        for i = 1:size(E1,1)
            iT = edgeAttachments(T,E1(i,1),E1(i,2));
            iT = iT{:};
            if length(iT) < 2 || norm(N(iT(1),:) - N(iT(2),:)) > tol
                nE = nE+1;
                E2(nE,:) = E1(i,:);
            end
        end
    end
    

    %%% Face color
    azimuth = (atan2(N(:,2),N(:,1)) + pi)/(2*pi);
    S = azimuth;
    H = 0.6*ones(size(azimuth));
%     minI = find(N(:,3) > 0);
%     minZ = min(N(minI,3));
%     z0 = minZ/2;
    z0 = 0.2;
    elevation = 2*atan2((N(:,3)+1)/2, sqrt(N(:,1).^2+N(:,2).^2))/pi;
%     S = (1-z0)*azimuth+z0;
    L = (1-z0)*elevation+z0;
%     L = elevation;
%     L = 1/(1-z0)*elevation-z0/(1-z0);
%     hsl = [H, S, L];
%     hsl = [ones(size(H)), (L+H)/2, L];
    hsl = [H, S, L];

    %%% Plot settings
    hold on
    view(2)
%     D = max(P) - min(P);
%     pbaspect([D(1) D(2) min(D(1:2))])
%     pbaspect([D(1) D(2) max(D(3),1)])
    nlim = [min(P(:,1)) max(P(:,1))];
    mlim = [min(P(:,2)) max(P(:,2))];
    xlim(nlim + [-1 1])
    ylim(mlim + [-1 1])

    %%% Plot faces
    patch('Faces',k,'Vertices',P,'FaceVertexCData',hsl2rgb(hsl),'FaceColor','flat','LineStyle','none')

    %%% Plot edges
    if isplanar
        plot3(P(k,1),P(k,2),P(k,3),'.-k','LineWidth',lineWidth)
    else
        for i = 1:size(E2,1)
            pE = P(E2(i,:),:);
            plot3(pE(:,1),pE(:,2),pE(:,3),'.-k','LineWidth',lineWidth,'MarkerSize',10)
        end
    end

    %%% Plot flow vectors
    nD = 0;
    dD = zeros(0,2);
    dS = zeros(0,2);
    for i = 1:nP
        hasSubdominant = false;
        imax = i;
        for j = i+1:nP
            if P(i,1) == P(j,1) && P(i,2) == P(j,2)
                nD = nD + 1;
                if P(i,3) > P(j,3)
                    dD(nD,:) = dP(i,:);
                    dS(nD,:) = dP(j,:);
                elseif P(i,3) < P(j,3)
                    dD(nD,:) = dP(j,:);
                    dS(nD,:) = dP(i,:);
                    imax = j;
                else
                    dD(nD,:) = dP(i,:) + dP(j,:);
                    dS(nD,:) = [0 0];
                end
                hasSubdominant = true;
                break
            end
        end
        if ~hasSubdominant
            nD = nD + 1;
            dD(nD,:) = dP(i,:);
            dS(nD,:) = [0 0];
        end
%         if showLabels
%             FilledArrow3([P(i,1); P(i,2); max(P(:,3)) + 0.01],dD(nD,:)',dS(nD,:)',xlim,ylim,'k')
%         end
    end
    if showLabels
        for i = 1:length(ku)
            FilledArrow3([P(ku(i),1); P(ku(i),2); max(P(ku,3)) + tol],dD(ku(i),:)',dS(ku(i),:)',xlim,ylim,'k')
        end
    end
    hold off
    xticks(nlim(1):nlim(2))
    yticks(mlim(1):mlim(2))

end