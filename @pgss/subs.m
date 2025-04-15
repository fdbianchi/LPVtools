function [sys,msg] = subs(obj,par)

% SUBS(pdG,par) evaluates an LPV model pdG at a frozen parameter value PAR
%
% Use:
%   sys = SUBS(pdG,par)
%
% where
%   - pdG:  a pass object
%   - par:  parameter value
%   - sys:  ss-object
%
% See also pgss, pset

% fbianchi - 2021-03-31

if (nargin < 2)
    error('PGSS:SUBS:notEnoughInputs','use sys = subs(pdG,par)')
end

% check if PAR is consistent with the parameter set
[bool,msg] = checkval(obj.parset,par);

if (bool > 0)

    n   = size(par,2);
    
    % model dimensions
    [no,ni,ns] = size(obj);

    % pre-allocate memory
    a = zeros(ns,ns,n);
    b = zeros(ns,ni,n);
    c = zeros(no,ns,n);
    d = zeros(no,ni,n);
    
    for ii = 1:n

        % parameter evaluated at parfcn
        fpar = obj.parfcn(par(:,ii));
        
        % independent terms
        a(:,:,ii) = obj.A(:,:,1);   
        b(:,:,ii) = obj.B(:,:,1);
        c(:,:,ii) = obj.C(:,:,1);   
        d(:,:,ii) = obj.D(:,:,1);

        for jj = 1:nsys(obj)-1
            % parameter dependant terms
            a(:,:,ii) = a(:,:,ii) + obj.A(:,:,jj+1)*fpar(jj);
            b(:,:,ii) = b(:,:,ii) + obj.B(:,:,jj+1)*fpar(jj);
            c(:,:,ii) = c(:,:,ii) + obj.C(:,:,jj+1)*fpar(jj);
            d(:,:,ii) = d(:,:,ii) + obj.D(:,:,jj+1)*fpar(jj);
        end
    end
    
    % new SS object with the same names
    warning('off','Control:ltiobject:RepeatedChannelNames')
    sys = ss(a,b,c,d,...
        'StateName',obj.StateName,...
        'InputName',obj.InputName,...
        'OutputName',obj.OutputName,...
        'InputGroup',obj.InputGroup,...
        'OutputGroup',obj.OutputGroup);
    warning('on','Control:ltiobject:RepeatedChannelNames')

else
    if (nargout < 2)
        error('PGSS:SUBS:outputError',msg)
        
    else
        sys = [];
        
    end        
    
end
