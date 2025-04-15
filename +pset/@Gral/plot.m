function plot(obj)

% PLOT illustrates a parameter set 
%
% Use: 
%   plot(set)
%
% For 1, 2 or 3 parameters.
%
% See also pset.Gral

% fbianchi - 2021-03-30


% set dimensions
[np,nv] = size(obj);

% figure setting
pColor = [0.4940    0.1840    0.5560];
aColor = pColor;
delete(gca);
Ha = axes();
Ha.NextPlot = 'add';

if (np == 1)
    % 1D case

    % plot limits
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

    % plot limits
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
    k = obj.hullIndex;
    xa = obj.points(1,k);
    ya = obj.points(2,k);
    patch(xa,ya,aColor,...
        'EdgeColor',aColor,...
        'FaceAlpha',0.3);
    
    % set points
    x = obj.points(1,:);
    y = obj.points(2,:);
    if (nv > np)
        triplot(obj.simplices,x,y,'Color',[pColor 0.6]);
    end
    plot(Ha,x,y,'Marker','s',...
                'MarkerEdgeColor',pColor,...
                'MarkerFaceColor',pColor,...
                'LineStyle','none')
    for ii = 1:nv
        text(obj.points(1,ii) + 0.02*xd,...
             obj.points(2,ii) + 0.04*yd,...
            sprintf('v%1.0f',ii))
    end
    
    % set range
    plot(obj.range(1,[1 2 2 1 1]),obj.range(2,[1 1 2 2 1]),'k:')

    % labels
    xlabel(obj.ParameterNames{1})
    ylabel(obj.ParameterNames{2})

            
elseif (np == 3)
    % 3D case

    % plot limits
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
    face = obj.hullIndex;
    vert = obj.points';
    patch('Vertices',vert,'Faces',face,...
        'FaceColor',aColor,...
        'EdgeColor',aColor,...
        'FaceAlpha',0.3);
    view(3);
    
    % set points
    x = obj.points(1,:);
    y = obj.points(2,:);
    z = obj.points(3,:);
    if (nv > np)
        trimesh(obj.simplices,x,y,z,'EdgeColor',pColor,'FaceAlpha',0.6);
    end
    plot3(Ha,x,y,z,'Marker','s',...
                'MarkerEdgeColor',pColor,...
                'MarkerFaceColor',pColor,...
                'LineStyle','none')
    for ii = 1:nv
        text(obj.points(1,ii) + 0.02*xd,...
             obj.points(2,ii) + 0.04*yd,...
             obj.points(3,ii) + 0.04*zd,...
             sprintf('v%2.0f',ii))
    end
    
    % set range
    plot3(obj.range(1,[1 2 2 1 1]),obj.range(2,[1 1 2 2 1]),obj.range(3,[1 1 1 1 1]),'k:')
    plot3(obj.range(1,[1 2 2 1 1]),obj.range(2,[1 1 2 2 1]),obj.range(3,2*[1 1 1 1 1]),'k:')
    plot3(obj.range(1,[1 1]),obj.range(2,[1 1]),obj.range(3,[1 2]),'k:')
    plot3(obj.range(1,[1 1]),obj.range(2,2*[1 1]),obj.range(3,[1 2]),'k:')
    plot3(obj.range(1,2*[1 1]),obj.range(2,[1 1]),obj.range(3,[1 2]),'k:')
    plot3(obj.range(1,2*[1 1]),obj.range(2,2*[1 1]),obj.range(3,[1 2]),'k:')

    % labels
    xlabel(obj.ParameterNames{1})
    ylabel(obj.ParameterNames{2})
    zlabel(obj.ParameterNames{3})
    
else
    error('PSET:GRAL:PLOT:inputError','PLOT is no available for PSET of more than 3 parameters')
end






