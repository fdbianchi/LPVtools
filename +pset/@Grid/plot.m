function plot(obj)

% PLOT illustrates a parameter set 
%
% Use: 
%   plot(set)
%
% For 1, 2 or 3 parameters.
%
% See also PSET.GRID

% fbianchi - 2021-03-30

% set dimensions
[np,nv] = size(obj);

% figure setting
pColor = [0.8500    0.3250    0.0980];
aColor = pColor;
delete(gca);
Ha = axes();
Ha.NextPlot = 'add';

if (np == 1)
    % 1D case

    % set limits
    xd = obj.range(1,2) - obj.range(1,1);   
    if (xd == 0)
        xo = 0.1;
    else
        xo = 0.1*xd;
    end
    yd = 1;                                 %yo = 0.1*yd; 
   
    Ha.XLim = obj.range(1,:) + xo*[-1 1];
    Ha.YLim = 0.5*[-1 1];
                
    % set area
    xa = obj.range(1,[1 2]);
    ya = [0 0];
    patch(xa,ya,aColor,...
        'EdgeColor',aColor,...
        'FaceAlpha',0.3);
    
    % set points
    x = obj.points(1,:);
    y = zeros(1,nv);
    plot(Ha,x,y,'Marker','s',...
                'MarkerEdgeColor',pColor,...
                'MarkerFaceColor',pColor,...
                'LineStyle','none')
    for ii = 1:nv
        text(obj.points(1,ii) + 0.02*xd, 0.04*yd,...
            sprintf('v%1.0f',ii))
    end
    xlabel(obj.ParameterNames{1})
    
elseif (np == 2)
    % 2D case

    % set limits
    xd = obj.range(1,2) - obj.range(1,1);
    if (xd == 0)
        xo = 0.1;
    else
        xo = 0.1*xd;
    end
    yd = obj.range(2,2) - obj.range(2,1); 
    if (yd == 0)
        yo = 0.1;
    else
        yo = 0.1*yd;
    end
    
    Ha.XLim = obj.range(1,:) + xo*[-1 1];
    Ha.YLim = obj.range(2,:) + yo*[-1 1];
                
    % set area
    xa = obj.range(1,[1 2 2 1 1]);
    ya = obj.range(2,[1 1 2 2 1]);
    patch(xa,ya,aColor,...
        'EdgeColor',aColor,...
        'FaceAlpha',0.3);
    
    % set points
    x = obj.points(1,:);
    y = obj.points(2,:);
    plot(Ha,x,y,'Marker','s',...
                'MarkerEdgeColor',pColor,...
                'MarkerFaceColor',pColor,...
                'LineStyle','none')
    for ii = 1:nv
        text(obj.points(1,ii) + 0.02*xd,...
             obj.points(2,ii) + 0.04*yd,...
            sprintf('v%1.0f',ii))
    end
    
    % point grid
    xt = unique(obj.points(1,:));
    for ii = 1:length(xt)
        plot(xt(ii)*[1 1],obj.range(2,:),'k:')
    end
    yt = unique(obj.points(2,:));
    for ii = 1:length(yt)
        plot(obj.range(1,:),yt(ii)*[1 1],'k:')
    end

    xlabel(obj.ParameterNames{1})
    ylabel(obj.ParameterNames{2})
    
%     Ha.NextPlot = 'replace';

            
elseif (np == 3)
    % 3D case

    % set limits
    xd = obj.range(1,2) - obj.range(1,1);
    if (xd == 0)
        xo = 0.1;
    else
        xo = 0.1*xd;
    end
    yd = obj.range(2,2) - obj.range(2,1);
    if (yd == 0)
        yo = 0.1;
    else
        yo = 0.1*yd;
    end
    zd = obj.range(3,2) - obj.range(3,1);
    if (zd == 0)
        zo = 0.1;
    else
        zo = 0.1*zd;
    end
   
    Ha.XLim = obj.range(1,:) + xo*[-1 1];
    Ha.YLim = obj.range(2,:) + yo*[-1 1];
    Ha.ZLim = obj.range(3,:) + zo*[-1 1];
                
    % set area
    vert = pgrid(obj.range)';
    if (sum([xd yd  zd] == 0) == 3)
        face = 1;
    elseif (sum([xd yd  zd] == 0) == 2)
        face = [1 2];
    elseif (sum([xd yd  zd] == 0) == 1)
        face = [1 2 4 3 1];
    else
        face = [1 2 4 3; 5 6 8 7; 
                1 2 6 5; 2 4 8 6;
                4 3 7 8; 3 1 5 7];
    end
    patch('Vertices',vert,'Faces',face,...
        'FaceColor',aColor,...
        'EdgeColor',aColor,...
        'FaceAlpha',0.3);
    view(3);
    
    % set points
    x = obj.points(1,:);
    y = obj.points(2,:);
    z = obj.points(3,:);
    plot3(Ha,x,y,z,'Marker','s',...
                'MarkerEdgeColor',pColor,...
                'MarkerFaceColor',pColor,...
                'LineStyle','none')
    for ii = 1:nv
        text(obj.points(1,ii) + 0.02*xd,...
             obj.points(2,ii) + 0.04*yd,...
             obj.points(3,ii) + 0.04*zd,...
            sprintf('v%1.0f',ii))
    end
    
    % point grid
    xt = unique(obj.points(1,:));
    yt = unique(obj.points(2,:));
    zt = unique(obj.points(3,:));
    for ii = 1:length(xt)
        for jj = 1:length(zt)
            plot3(xt(ii)*[1 1],obj.range(2,:),zt(jj)*[1 1],'k:')
        end
    end
    for ii = 1:length(yt)
        for jj = 1:length(zt)
            plot3(obj.range(1,:),yt(ii)*[1 1],zt(jj)*[1 1],'k:')
        end
    end
    for ii = 1:length(xt)
        for jj = 1:length(yt)
            plot3(xt(ii)*[1 1],yt(jj)*[1 1],obj.range(3,:),'k:')
        end
    end

    xlabel(obj.ParameterNames{1})
    ylabel(obj.ParameterNames{2})
    zlabel(obj.ParameterNames{3})
    
else
    error('PLOT:inputError','PLOT is no available for PSET of more than 3 parameters')
end






