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
% See also ppss, pset

% fbianchi - 2020-02-20

if (nargin < 2)
    error('PPSS:SUBS:notEnoughInputs','use sys = subs(pdG,par)')
end

% dimensions
n = size(par,2);
[no,ni,ns,~,nv] = size(obj);

if (nv == 1)
    % trivial case (one model)
    warning('off','Control:ltiobject:RepeatedChannelNames')
    sys = ss(obj.A,obj.B,obj.C,obj.D,...
        'StateName',obj.StateName,...
        'InputName',obj.InputName,...
        'OutputName',obj.OutputName,...
        'InputGroup',obj.InputGroup,...
        'OutputGroup',obj.OutputGroup);
    warning('on','Control:ltiobject:RepeatedChannelNames')
    msg = [];
    
else
    
    % check if PAR is consistent with the parameter set
    [bool,msg] = checkval(obj.parset,par);
    
    if (bool > 0)
        
        % memory allocation
        a = zeros(ns,ns,n);
        b = zeros(ns,ni,n);
        c = zeros(no,ns,n);
        d = zeros(no,ni,n);
        
        for ii = 1:n
            
            % in case of polytopic first it is necessary to obtain a convex
            % decomposition
            try
                [alpha, idx] = cvxdec(obj.parset,par(:,ii));
            catch ME
                error('PPSS:SUBS:inputError',ME.message);
            end
            
            % using reshape to implement
            %   S = sum_{i=1}^{n} alpha_i*S_i
            nv = length(alpha);
            a(:,:,ii) = reshape(reshape(obj.A(:,:,idx),ns*ns,nv,1)*alpha,ns,ns,1);
            b(:,:,ii) = reshape(reshape(obj.B(:,:,idx),ns*ni,nv,1)*alpha,ns,ni,1);
            c(:,:,ii) = reshape(reshape(obj.C(:,:,idx),no*ns,nv,1)*alpha,no,ns,1);
            d(:,:,ii) = reshape(reshape(obj.D(:,:,idx),ni*no,nv,1)*alpha,no,ni,1);
            
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
            error('PPSS:SUBS:outputError',msg)
        else
            sys = [];
        end
        
    end
end


