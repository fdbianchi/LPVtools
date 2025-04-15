function newSet = subsref(oldSet,S)

% SUBSREF function for pset.Grid
%
%   newSet = oldSet(idx)
%
% returns a subset corresponding to the subscripts in idx
%
% See also pset.Grid

% fbianchi - 2021-03-30

if (nargin == 1)
    newSet = oldSet;
    return
end

switch S(1).type
    case '()'
        
        % only the first subscript is used    
        if (length(S(1).subs) == 1)
            idx = S(1).subs{1};
        else
            error('PSET:GRID:SUBSREF:inputError',...
                'Too many subscripts')
        end
        
        % checking if the subscript is in range
        np = size(oldSet);
        if ~(isnumeric(idx) || strcmp(idx,':'))
            error('PSET:GRID:SUBSREF:inputError',...
                'The subscript must be numeric')
        elseif (~strcmp(idx,':') && (max(idx) > np || min(idx) <= 0))
            error('PSET:GRID:SUBSREF:inputError',...
                'Index exceeds set dimensions.')
        end

        % new set
        newSet = pset.Grid;
        % removing repeated points
        newSet.points = unique(oldSet.points(idx,:)','rows','stable')';
        % parameter rate
        if ~isempty(oldSet.rate)
            newSet.rate = oldSet.rate(idx,:);
        else
            newSet.rate = oldSet.rate;
        end    
        % names
        newSet.ParameterNames = oldSet.ParameterNames(idx);
        
        
    case '.'
        % leave other subref without definition
        newSet = builtin('subsref',oldSet,S);
        
    otherwise
        error('PSET:GRID:SUBSREF:inputError','Indexing type not supported')

end