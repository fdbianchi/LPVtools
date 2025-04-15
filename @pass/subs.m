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
% See also pass, pset

% fbianchi - 2021-03-31

if (nargin < 2)
    error('PASS:SUBS:notEnoughInputs','use sys = subs(pdG,par)')
end

% # of points
n = size(par,2);
% model dimensions
[no,ni,ns,~,nv] = size(obj);

if (nv == 1)
    % trivial case (one model)
    sys = ss(obj.A,obj.B,obj.C,obj.D,...
        'StateName',obj.StateName,...
        'InputName',obj.InputName,...
        'OutputName',obj.OutputName,...
        'InputGroup',obj.InputGroup,...
        'OutputGroup',obj.OutputGroup);
    msg = [];

else

    % check if PAR is consistent with the parameter set
    [bool,msg] = checkval(obj.parset,par);
    
    if (bool > 0)
        
        % using reshape to implement
        %   S = S0 + sum_{i=1}^{n} par_i*S_i
        auxPar = [ones(1,n);par];
        a = reshape(reshape(obj.A,ns*ns,nv,1)*auxPar,ns,ns,n);
        b = reshape(reshape(obj.B,ns*ni,nv,1)*auxPar,ns,ni,n);
        c = reshape(reshape(obj.C,no*ns,nv,1)*auxPar,no,ns,n);
        d = reshape(reshape(obj.D,ni*no,nv,1)*auxPar,no,ni,n);
        
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
            error('PASS:SUBS:outputError',msg)
        else
            sys = [];
        end
        
    end
end


