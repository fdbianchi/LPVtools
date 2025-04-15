function fout = multfunc(fcn,parset)

% computing the new parfcn when there are cross products

if ~isa(fcn,'function_handle')
    error('SUBFUNC:inputerror','oldFcn must be an anonymous function')
end

try
    nv = length(fcn(parset.points(:,1)));
catch err
    error('PGSS:PGSS:inputError',...
        'Function definition does not correspond with the parameter set')
end
    
% string with old function
strFcn = func2str(fcn);

% get the argument
idx1 = find(strFcn=='@',1,'first');
idx2 = find(strFcn==')',1,'first');
strCtrlFcn = strFcn(idx1:idx2);

% elements in oldFcn
idx1 = find(strFcn=='[',1,'first') + 1;
idx2 = find(strFcn==']',1,'first') - 1;
strElems = strsplit(strFcn(idx1:idx2),{';',',',' '});
ne = length(strElems);

if (ne ~= nv)
    
    % "linear" terms
    strFcn = ['@(f)[' sprintf('f(%d);',1:nv)];
    % "multi-affine" terms
    idx = nchoosek(1:nv,2);
    for ii = 1:size(idx,1)
        strFcn = [strFcn, sprintf('f(%d).*f(%d);',idx(ii,:))];
    end
    % "quadratic" terms
    for jj = 1:nv
        strFcn = [strFcn, sprintf('f(%d).^2;',jj)];
    end
    strFcn(end) = ']';
    f = str2func(strFcn);
    fout =@(p) f(fcn(p));
    
else
    
    % "linear" terms
    fi = strElems;
    
    % "multi-affine" terms
    idx = nchoosek(1:ne,2);
    for ii = 1:size(idx,1)
        fi{ii+ne} = [strElems{idx(ii,1)} '*' strElems{idx(ii,2)}];
    end
    
    % "quadratic" terms
    for jj = 1:ne
        ii = ne + ii + 1;
        fi{ii} = [strElems{jj} '^2'];
    end
    
    strNewFcn = [strCtrlFcn ' ['];
    for ii = 1:length(fi)
        strNewFcn = [strNewFcn fi{ii} ';'];
    end
    strNewFcn(end) = ']';
    
    fout = str2func(strNewFcn);
end