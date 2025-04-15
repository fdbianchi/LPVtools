function disp(obj)

% DISP method for class pset.Hull

% fbianchi - 2021-03-30

if isempty(obj)
    str = '  Empty pset.Hull object';
    
else
    % dimensions of the set
    [np,nv] = size(obj);
    
    str = sprintf('Convex hull parameter set of %d parameters and %d points:\n',...
                    np,nv);
    if isempty(obj.rate)
        for ii = 1:np
            str = [str, sprintf('\t%s:\t range: [%-g, %-g]\n',...
                        obj.ParameterNames{ii},obj.range(ii,:))];
        end
    else
        for ii = 1:np
            str = [str, sprintf('\t%s:\t range: [%-g, %-g],\t\trate: [%-g, %-g]\n',...
                        obj.ParameterNames{ii},obj.range(ii,:),obj.rate(ii,:))];
        end
    end        
end

disp(str)

end

 