function [lmis,obj,lmisChecked] = constH2g2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% *** for internal use ***
%
% CONSTH22LMI returns a set of LMIs for generalized H2 performance constraint
%
% Use:
%   [lmis,obj] = CONSTH2g2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% fbianchi - 2021-06-29

% objective initial value
obj = 0;

% optimization settings
eigtol  = options.eigtol;

% number of systems in the description
nv = size(synSet.LPVatPoints,3);

% system order
ns = order(pdG);

% disturbance and performance input-output map
ioChannel.perf = constraint.outputs;
nz = length(ioChannel.perf);
ioChannel.dist = constraint.inputs;
nw = length(ioChannel.dist);

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
                % ToDo: implement prjected versions
                
            case 'prA'
                % ToDo: implement prjected versions
                
            case 'basic'
                
                track('H2g: Basic (change of variables)',options.debug)
                
                for ii = 1:nv
                    
                    % system matrices at par
                    [A,B1,b2,C1,c2,D11,d12,d21] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
                    
                    if (all(all(d12)) == 0 || all(all(d21)) == 0) && any(any(D11))
                        error('LPVSYN:constH2g2LMI:LMI',...
                            'The constraint Dcl=0 is infeasible')
                    end
                    
                    % controller vars at par
                    [Aki,Bki,Cki,Dki,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
                    % ToDo: Implement constraints for forcing ctrller
                    % parameter dependence
                    
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
                        
                        % H2g constraints
                        matXY = blkvar;
                        matXY(1,1) = (-dy/2 + A*Y + b2*Cki) + (-dy/2 + A*Y + b2*Cki)';
                        matXY(2,1) = (Aki + (A + b2*Dki*c2)');
                        matXY(2,2) = (dx/2 + X*A + Bki*c2) + (dx/2 + X*A + Bki*c2)';
                        matXY(3,1) = (B1 + b2*Dki*d21)';
                        matXY(3,2) = (X*B1 + Bki*d21)';
                        matXY(3,3) = -eye(nw);
                        matXY = sdpvar(matXY);
                        lmiName = sprintf('H2g: LMIXY @(%d) - Basic',ii);
                        if isnumeric(matXY)
                            aux.lmiValue = all(eig(matXY) <= -eigtol);
                            aux.lmiName = lmiName;
                            lmisChecked = [lmisChecked, aux];
                        else
                            lmis = [lmis, (matXY <= -eigtol*eye(size(matXY))):lmiName];
                        end
                        
                    end
                    
                    lmiMu = blkvar;
                    lmiMu(1,1) = Y;
                    lmiMu(1,2) = eye(ns);
                    lmiMu(1,3) = (C1*Y + d12*Cki)';
                    lmiMu(2,2) = X;
                    lmiMu(2,3) = (C1 + d12*Dki*c2)';
                    lmiMu(3,3) = mu2*eye(nz);
                    lmiMu = sdpvar(lmiMu);
                    
                    lmiName = sprintf('H2g: LMI_mu @(%d) - Basic',ii);
                    if isnumeric(lmiMu)
                        aux.lmiValue = all(eig(lmiMu) <= -eigtol);
                        aux.lmiName = lmiName;
                        lmisChecked = [lmisChecked, aux];
                    else
                        lmis = [lmis, (lmiMu >= eigtol*eye(size(lmiMu))):lmiName];
                    end
                    
                    % Dcl=0
                    Dcl = D11 + d12*Dki*d21;
                    lmiName = sprintf('H2g: Dcl==0 @(%d) - Basic',ii);
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
        
        track('H2g: State feedback',options.debug)
        
        lmis = [];
        for ii = 1:nv
            
            % system matrices at par
            [A,B1,b2,C1,~,D11,d12] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
            
            if any(any(D11))
                error('LPVSYN:constH2g2LMI:LMI',...
                    'The constraint Dcl=0 is infeasible')
            end
            
            % controller vars at par
            [~,~,~,Wi] = evalCtrlVars(ctrlSet,ii);
            
            % Lyapunov matrices at par
            [Y,~,dY] = evalLyapFcn(lyapSet,ii);
            
            for jj = 1:length(dY)
                
                if isscalar(dY)
                    dy = 0;
                else
                    dy = dY{jj};
                end
                
                % H2 constraints
                lmiY = blkvar;
                lmiY(1,1) = dy + ((A*Y + b2*Wi) +(A*Y + b2*Wi)');
                lmiY(2,3) = B1;
                lmiY(3,3) = -eye(nw);
                lmiY = sdpvar(lmiY);
                lmiName = sprintf('H2g: LMIY @(%d) - SF',ii);
                if isnumeric(matY)
                    aux.lmiValue = all(eig(lmiY) <= -eigtol);
                    aux.lmiName = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmis = [lmis, (lmiY <= -eigtol*eye(size(lmiY))):lmiName];
                end
                
                lmiMu = blkvar;
                lmiMu(1,1) = mu2*eye(nw);
                lmiMu(1,2) = C1 + d12*Wi;
                lmiMu(2,2) = Y;
                lmiMu = sdpvar(lmiMu);
                lmiName = sprintf('H2g: LMI_mu @(%d) - SF',ii);
                if isnumeric(lmiMu)
                    aux.lmiValue = all(eig(lmiMu) <= -eigtol);
                    aux.lmiName = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmis = [lmis, (lmiMu >= eigtol*eye(size(lmiMu))):lmiName];
                end
                
            end
            
        end
        
    otherwise
        error('LPVSYN:constH2g2LMI:LMI',...
            'Type of feedback invalid')
        
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

