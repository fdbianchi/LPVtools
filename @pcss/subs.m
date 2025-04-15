function [K,msg] = subs(obj,par)

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
% See also pcss

% fbianchi - 2021-07-02

if (nargin < 2)
    error('PCSS:SUBS:notEnoughInputs','use sys = subs(pdG,par)')
end

% # of points
n = size(par,2);

% check if PAR is consistent with the parameter set
[bool,msg] = checkval(obj.ctrller.parset,par);

if (bool > 0)
    
    if isempty(obj.xFcn)    % State feedback
        
        % Lyapunov function
        if isnumeric(obj.yFcn)
            Y = obj.yFcn;
        else
            Y = subs(obj.yFcn,par);
            Y = Y.D;
        end
        
        % Controller variables at p
        Khat = subs(obj.ctrller,par);
        [no,ni] = size(obj.ctrller);
        
        Dk = Khat.D;
        
        % undoing variable change
        dk = zeros(no,ni,n);
        for ii = 1:n
            dk(:,:,ii) = Dk(:,:,ii)/Y(:,:,ii);
        end
        
        K = ss(dk);
        
    else    % Output feedback
        
        % model dimensions
        [no,ni,ns] = size(obj.ctrller);
        
        % pre-allocate memory
        ak = zeros(ns,ns,n);
        bk = zeros(ns,ni,n);
        ck = zeros(no,ns,n);
        dk = zeros(no,ni,n);
        
        for ii = 1:n
            
            % Lyapunov function factorization
            if isnumeric(obj.xFcn) && (size(obj.xFcn,3) == 1)
                X = obj.xFcn;
            else
                [~,~,~,X] = ssdata(subs(obj.xFcn,par(:,ii)));
            end
            
            if isnumeric(obj.yFcn) && (size(obj.yFcn,3) == 1)
                Y = obj.yFcn;
            else
                [~,~,~,Y] = ssdata(subs(obj.yFcn,par(:,ii)));
            end
            
            ns = size(X,1);
            
            % Controller variables at p
            Khat = subs(obj.ctrller,par(:,ii));
            [Ak,Bk,Ck,Dk] = ssdata(Khat);
            
            % Plant matrices at p
            G22 = subs(obj.plant,par(:,ii));
            [A,B2,C2] = ssdata(G22);

            % Compute N & M
            if isnumeric(obj.xFcn)
                Mt = (eye(ns) - X*Y);
                % N = eye(ns);
                
                % undoing variable change
                ak(:,:,ii) = (Ak - X*(A-B2*Dk*C2)*Y - Bk*C2*Y - X*B2*Ck)/Mt;
                bk(:,:,ii) = (Bk - X*B2*Dk);
                ck(:,:,ii) = (Ck - Dk*C2*Y)/Mt;

            elseif isnumeric(obj.yFcn)
                % Mt  = eye(ns);
                N = (eye(ns) - X*Y);
                
                % undoing variable change
                ak(:,:,ii) = N\(Ak - X*(A-B2*Dk*C2)*Y - Bk*C2*Y - X*B2*Ck);
                bk(:,:,ii) = N\(Bk - X*B2*Dk);
                ck(:,:,ii) = (Ck - Dk*C2*Y);

            else
                [U,S,V] = svd(eye(ns) - X*Y);
                Si = diag(1./diag(S))^(1/2);
                Ni = Si*U';
                Mti = V*Si;
                % Ni = eye(ns)/(eye(ns) - X*Y);
                % Mti = eye(ns);
                
                % undoing variable change
                ak(:,:,ii) = Ni*(Ak - X*(A-B2*Dk*C2)*Y - Bk*C2*Y - X*B2*Ck)*Mti;
                bk(:,:,ii) = Ni*(Bk - X*B2*Dk);
                ck(:,:,ii) = (Ck - Dk*C2*Y)*Mti;
                
            end
            
            dk(:,:,ii) = Dk;
            
        end
        K = ss(ak,bk,ck,dk);
    end
else
    if (nargout < 2)
        error('PCSS:SUBS:outputError',msg)
    else
        K = [];
    end
end

