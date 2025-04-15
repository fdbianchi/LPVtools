function newFcn = subfunc(oldFcn,index,ne)

% newFcn = SUBFUNC(oldFcn,index) returns an anonymous functions with only the
% elements indicated in index
%
% Example:
%   fcnOld =@(x) [x(1);x(2);x(1)^2;x(1)*x(2)];
%
%   newFcn = SUBFUNC(oldFcn,[3 2])
% 
% produce:
%
%   newFcn =@(x) [x(1)^2;x(2)];
%

% fbianchi - 2021-03-23

if (nargin < 2)
    error('SUBFUNC:inputerror','insufficient input arguments')
end

if ~isa(oldFcn,'function_handle')
    error('SUBFUNC:inputerror','oldFcn must be an anonymous function')
end

if isempty(index) 
    index = 0;
elseif ~isnumeric(index) || ~isvector(index) || (min(index) < 0)
    error('SUBFUNC:inputerror','index must be a numeric vector with elements >= 0')
end
    
try 
    % string with old function
    strFcn = func2str(oldFcn);

    % getting function workspace (not recommended by Matlab documentation
    % see help functions)
    fInfo  = functions(oldFcn);
    workspacevars = fInfo.workspace{1};
    varnames = fieldnames(workspacevars);
    for ii = 1:length(varnames)
        evalstr = sprintf('%s = %s;', varnames{ii}, func2str(workspacevars.(varnames{ii})));
        eval(evalstr);
    end
    
    % get the argument
    idx1 = find(strFcn=='@',1,'first');
    idx2 = find(strFcn==')',1,'first');
    strCtrlFcn = strFcn(idx1:idx2);

    % elements in oldFcn
    idx1 = find(strFcn=='[',1,'first') + 1;
    idx2 = find(strFcn==']',1,'first') - 1;
    strElems = strsplit(strFcn(idx1:idx2),{';',',',' '});

    ne = length(strElems);
    if (max(index) > ne)
        error('SUBFUNC:inputerror','index elements must be <= %2.0f',ne)
    end

    % selecting elements of oldFcn
    if (length(index) == 1) && (index == 0)
        strCtrlFcn = [strCtrlFcn '[]'];
    else
        index(index==0) = [];
        strCtrlFcn = [strCtrlFcn ' ['];
        for ii = index
            strCtrlFcn = [strCtrlFcn strElems{ii} ';'];
        end
        strCtrlFcn(end) = ']';
    end
    newFcn = eval(strCtrlFcn);
    
    % check if it works
    aux = newFcn(ones(ne,1));
    
catch
    % in case the other method did not work
    % the problem with this option is that the results is not easy to
    % understand
    if (max(index) > ne)
        error('SUBFUNC:inputerror','index elements must be <= %2.0f',ne)
    end
    I = eye(ne);
    f = oldFcn;
    
    newFcn =@(p) I(index,:)*f(p);
    
end

