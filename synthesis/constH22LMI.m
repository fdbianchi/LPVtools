function [lmis,obj,lmisChecked,Q] = constH22LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% *** for internal use ***
%
% CONSTH22LMI returns a set of LMIs for H2 performance constraint
%
% Use:
%   [lmis,obj] = CONSTH22LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% fbianchi - 2021-06-29

% objective initial value
obj = 0;

% optimization settings
eigtol = options.eigtol;

% number of systems in the description
nv = size(synSet.LPVatPoints,3);

% system order
ns = order(pdG);

% disturbance and performance input-output map
ioChannel.perf = constraint.outputs;
nz = length(ioChannel.perf);
ioChannel.dist = constraint.inputs;
nw = length(ioChannel.dist);

% additional variable
Q = sdpvar(nz);

% mu is an optimization variable or not depending on soft or hard
% constraints
if isa(constraint.bound,'sdpvar')
    mu2 = constraint.bound;
else
    mu2 = constraint.bound^2;
end

% LMIs
lmis = [];
lmisChecked = [];
switch ctrlSet.fback             % type of feedback
    
    case 'output'
        
        switch synSet.type       % type of projection
            
            case 'prFull'
                % ToDo: Implement the projected versions
                
            case 'prA'
                % ToDo: Implement the projected versions
                
                
            case 'basic'
                
                track('H2: Basic (change of variables)',options.debug)
                
                
                for ii = 1:nv
                    
                    % system matrices at par
                    [A,B1,b2,C1,c2,D11,d12,d21] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
                    
                    if (all(all(d12)) == 0 || all(all(d21)) == 0) && any(any(D11))
                        error('LPVSYN:constH22LMI:LMI',...
                            'The constraint Dcl=0 is infeasible')
                    end
                    
                    % controller vars at par
                    [Aki,Bki,Cki,Dki,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
                    % ToDo: Implement constraint to force ctrller p
                    % dependence
                    
                    % Lyapunov matrices at par
                    [Y,X,dY,dX] = evalLyapFcn(lyapSet,ii);
                    
                    for jj = 1:max(length(dY),length(dX))
                        
                        if isscalar(dY)
                            dy = 0;
                        else
                            dy = dY{jj};
                        end
                        if isscalar(dX)
                            dx = 0;
                        else
                            dx = dX{jj};
                        end
                        
                        % H2 constraints
                        matXY = blkvar;
                        matXY(1,1) = (-dy/2 + A*Y + b2*Cki) + (-dy/2 + A*Y + b2*Cki)';
                        matXY(2,1) = (Aki + (A + b2*Dki*c2)');
                        matXY(2,2) = (dx/2 + X*A + Bki*c2) + (dx/2 + X*A + Bki*c2)';
                        matXY(3,1) = (B1 + b2*Dki*d21)';
                        matXY(3,2) = (X*B1 + Bki*d21)';
                        matXY(3,3) = -eye(nw);
                        matXY = sdpvar(matXY);
                        lmiName = sprintf('H2: LMIXY @(%d) - Basic',ii);
                        if isnumeric(matXY)
                            aux.lmiValue = all(eig(matXY) <= -eigtol);
                            aux.lmiName  = lmiName;
                            lmisChecked = [lmisChecked, aux];
                        else
                            lmis = [lmis, (matXY <= -eigtol*eye(size(matXY))):lmiName];
                        end
                        
                    end
                    
                    lmiQ = blkvar;
                    lmiQ(1,1) = Y;
                    lmiQ(1,2) = eye(ns);
                    lmiQ(1,3) = (C1*Y + d12*Cki)';
                    lmiQ(2,2) = X;
                    lmiQ(2,3) = (C1 + d12*Dki*c2)';
                    lmiQ(3,3) = Q;
                    lmiQ = sdpvar(lmiQ);
                    lmiName = sprintf('H2: LMIQ @(%d) - Basic',ii);
                    if isnumeric(lmiQ)
                        aux.lmiValue = all(eig(lmiQ) <= -eigtol);
                        aux.lmiName = lmiName;
                        lmisChecked = [lmisChecked, aux];
                    else
                        lmis = [lmis, (lmiQ >= eigtol*eye(size(lmiQ))):lmiName];
                    end
                    
                    % Dcl=0
                    Dcl = D11 + d12*Dki*d21;
                    lmiName = sprintf('H2: Dcl==0 @(%d) - Basic',ii);
                    if isnumeric(Dcl)
                        aux.lmiValue = all(all(Dcl==0));
                        aux.lmiName = lmiName;
                        lmisChecked = [lmisChecked, aux];
                    else
                        lmis = [lmis, (Dcl == 0):lmiName];
                    end
                    
                end
        end 
        
    case 'state'
        
        track('H2: State feedback',options.debug)
        
        lmis = [];
        for ii = 1:nv
            
            % system matrices at par
            [A,B1,b2,C1,~,D11,d12] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
            
            if any(any(D11))
                error('LPVSYN:constH22LMI:LMI',...
                    'The constraint Dcl=0 is infeasible')
            end
            
            % controller vars at par
            [~,~,~,Wi] = evalCtrlVars(ctrlSet,ii);%par,ii);
            
            % Lyapunov matrices at par
            [Y,~,dY] = evalLyapFcn(lyapSet,ii);%par);
            
            for jj = 1:length(dY)
                
                if isscalar(dY)
                    dy = 0;
                else
                    dy = dY{jj};%(:,:,jj);
                end
                
                % H2 constraints
                lmiY = blkvar;
                lmiY(1,1) = dy + ((A*Y + b2*Wi) + (A*Y + b2*Wi)');
                lmiY(1,2) = B1;
                lmiY(2,2) = -eye(nw);
                lmiY = sdpvar(lmiY);
                lmiName = sprintf('H2: LMIY @(%d) - SF',ii);
                if isnumeric(lmiY)
                    aux.lmiValue = all(eig(lmiY) <= -eigtol);
                    aux.lmiName = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmis = [lmis, (lmiY <= -eigtol*eye(size(lmiY))):lmiName];
                end
                
                lmiQ = blkvar;
                lmiQ(1,1) = Q;
                lmiQ(1,2) = C1*Y + d12*Wi;
                lmiQ(2,2) = Y;
                lmiQ = sdpvar(lmiQ);
                lmiName = sprintf('H2: LMIQ @(%d) - SF',ii);
                if isnumeric(lmiQ)
                    aux.lmiValue = all(eig(lmiQ) >= eigtol);
                    aux.lmiName = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmis = [lmis, (lmiQ >= eigtol*eye(size(lmiQ))):lmiName];
                end
                
            end
            
        end
        
    otherwise
        error('LPVSYN:constH22LMI:LMI',...
            'Type of feedback invalid')
end


lmiName = 'H2: Trace(Q) < mu';
if isnumeric(lmiQ)
    aux.lmiValue = (trace(Q) <= mu2);
    aux.lmiName = lmiName;
    lmisChecked = [lmisChecked, aux];
else
    lmis = [lmis, (trace(Q) <= mu2):lmiName];
end

% objective function
if isa(constraint.bound,'sdpvar')
    % soft constraint
    obj = obj + mean(constraint.bound)*constraint.factor;
  
end


% =======================================================
% Local functions

function track(msg,dbg)

if dbg
    disp(msg)
end
