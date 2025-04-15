function disp(obj)

% DISP method for class pset.Box

% fbianchi - 2021-03-29

if isempty(obj)
    str = '  Empty pset.Box class';
    
else
    % dimensions of the set
    np = size(obj);
    
    str = sprintf('Box parameter set of %d parameter(s):\n',np);
    
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

 