function newSet = subsref(oldSet,S)

% SUBSREF function for pset.Box 
%
%   newSet = oldSet(idx)
%
% returns a subset corresponding to the subscripts in idx
%
% See also pset.Box

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
            error('PSET:BOX:SUBSREF:inputError',...
                'Too many subscripts')
        end

        % checking if the subscript is in range
        np = size(oldSet);
        if (~isnumeric(idx) && ~strcmp(idx,':'))
            error('PSET:BOX:SUBSREF:inputError',...
                'The subscript must be numeric')
            
        elseif (~strcmp(idx,':') && (max(idx) > np || min(idx) <= 0))
            error('PSET:BOX:SUBSREF:inputError',...
                'Index exceeds set dimensions.')
            
        end
        
        % new range
        range = oldSet.range(idx,:);
        % new rate
        if ~isempty(oldSet.rate)
            rate   = oldSet.rate(idx,:);
        else
            rate   = oldSet.rate;
        end    
        % new names
        names   = oldSet.ParameterNames(idx);
        
        % new set
        newSet = pset.Box(range, rate, names);
        

    case '.'
        % leave other subref without definition
        newSet = builtin('subsref',oldSet,S);
        
    otherwise
        error('PSET:BOX:SUBSREF:inputError','Indexing type not supported')

end