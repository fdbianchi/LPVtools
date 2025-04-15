function [lmis,obj,lmisChecked] = constHinf2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% *** for internal use ***
%
% CONSTHINF2LMI returns a set of LMIs for Hinfinity performance constraint
%
% Use: 
%   [lmis,obj] = CONSTHINF2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% fbianchi - 2020-04-23


% objective initial value
obj = 0;

% optimization settings
eigtol = options.eigtol;

% number of systems in the description
nv = size(synSet.LPVatPoints,3);

% number of terms
nt = nsys(pdG);

% disturbance and performance input-output map
ioChannel.perf = constraint.outputs;
nz = length(ioChannel.perf);
ioChannel.dist = constraint.inputs;
nw = length(ioChannel.dist);

% gamma is an optimization variable or not depending on soft or hard
% constraints
if isa(constraint.bound,'sdpvar')
    g2 = constraint.bound;
else
    g2 = constraint.bound^2;
end

% LMIs
lmis = [];
lmisChecked = [];
switch ctrlSet.fback            % type of feedback
    
    case 'output'
        
        switch synSet.type      % type of projection
            
            case 'prFull'
                % full elimination of controller matrices
                track('Hinf: Full projection',options.debug)
                
                % null spaces
                [~,~,b2,~,c2,~,d12,d21] = parsysdata(synSet.LPVatPoints(:,:,1),ioChannel);
                NY = null([b2;d12]');  NY = blkdiag(NY,eye(nw));
                NX = null([c2,d21]);   NX = blkdiag(NX,eye(nz));
                
                for ii = 1:nv
                    
                    % system matrices at par
                    [A,B1,C1,D11] = ssdata(synSet.LPVatPoints(ioChannel.perf,ioChannel.dist,ii));
                    
                    % Lyapunov matrices at par
                    [Y,X,dY,dX] = evalLyapFcn(lyapSet,ii);

                    for jj = 1:length(dX)
                        
                        matX = blkvar;
                        matX(1,1) = (0.5*dX{jj} + X*A) + (0.5*dX{jj} + X*A)';
                        matX(1,2) = X*B1;
                        matX(1,3) = C1';
                        matX(2,2) = -eye(nw);
                        matX(2,3) = D11';
                        matX(3,3) = -g2*eye(nz);
                        matX = sdpvar(matX);
                        auxMat = NX'*matX*NX;
                        auxMat = 0.5*(auxMat + auxMat'); % to avoid non-symmetric matrices
                        lmiName = sprintf('Hinf: BRL-X @(%d) - prFull',ii);
                        if isnumeric(auxMat)
                            aux.lmiValue = all(eig(auxMat) <= -eigtol);
                            aux.lmiName = lmiName;
                            lmisChecked = [lmisChecked, aux];
                        else
                            lmis = [lmis, (auxMat <= -eigtol*eye(size(auxMat))):lmiName];
                        end
                    end
                        
                    for jj = 1:length(dY)

                        matY = blkvar;
                        matY(1,1) = (-0.5*dY{jj} + A*Y) + (-0.5*dY{jj} + A*Y)';
                        matY(1,2) = Y*C1';
                        matY(1,3) = B1;
                        matY(2,2) = -g2*eye(nz);
                        matY(2,3) = D11;
                        matY(3,3) = -eye(nw);
                        matY = sdpvar(matY);
                        auxMat = NY'*matY*NY;
                        auxMat = 0.5*(auxMat + auxMat'); % to avoid non-symmetric matrices
                        lmiName = sprintf('Hinf: BRL-Y @(%d) - prFull',ii);
                        if isnumeric(auxMat)
                            aux.lmiValue = all(eig(auxMat) <= -eigtol);
                            aux.lmiName = lmiName;
                            lmisChecked = [lmisChecked, aux];
                        else
                            lmis = [lmis, (auxMat <= -eigtol*eye(size(auxMat))):lmiName];
                        end
                    end
                end
                % ToDo: check full proj when Lyapunov is affine
                
                
            case 'prA'
                % only Ak elimination
                % ToDo: Hinf: implement Ak projection
                
            case 'basic'
                % only change of variables

                track('Hinf: Basic (change of variables)',options.debug)
                
                lyapType = lyapSet.type;
                if strcmp(lyapSet.type(1:3), 'pwa') && (nv > 1)
                    % in case of PWA Lyapunov fcn, the affine constraints
                    % are repeated for each simplex
                    Simplices = pdG.parset.simplices;
                    if isempty(Simplices)
                        Simplices = 1;
                    end
                    nx = size(Simplices, 1);
                    lyapType(1:3) = 'aff';
                else
                    Simplices = 1:nv;
                    nx = 1;
                end
                
                for ll = 1:nx
                    
                    simplex = Simplices(ll, :);
                    
                    if ~strcmp(lyapType,'cte')
                        if strcmp(lyapSet.type(1:3), 'pwa') && (nv > 1)
                            % PWA: the affine constraints are based on local
                            % model for each simplex
                            pdGlocal = ppss2affloc(pdG, ll);
                            [~,~,~,~,Yaff,Xaff] = evalLyapFcn(lyapSet, 1, simplex);
                        else
                            pdGlocal = pdG;
                            Yaff = lyapSet.Y;
                            Xaff = lyapSet.X;
                        end
                        nt = nsys(pdGlocal);
                    end
                    
                    % extra constraints for affine cases
                    switch lyapType
                        
                        case 'aff'
                            
                            mu = sdpvar(1,nt-1,'full');
                            
                            for jj = 2:nt
                                [A,B1,~,C1] = parsysdata(pdGlocal(jj),ioChannel);
                                
                                Y = Yaff(:,:,jj);
                                matAffY = [A*Y + (A*Y)' Y*C1'; C1*Y zeros(nz)];
                                lmiName = sprintf('Hinf: aff-Y (term: %d:%d)',ll,jj);
                                lmis = [lmis, (matAffY + mu(jj-1)*eye(size(matAffY)) >= 0):lmiName];
                                
                                X = Xaff(:,:,jj);
                                matAffX = [X*A + (X*A)' X*B1; B1'*X zeros(nw)];
                                lmiName = sprintf('Hinf: aff-X (term: %d:%d)',ll,jj);
                                lmis = [lmis, (matAffX + mu(jj-1)*eye(size(matAffX)) >= 0):lmiName];
                                
                            end
                            lmiName = 'Hinf : aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                            
                        case {'affX','affdX'}
                            mu = sdpvar(1,nt-1,'full');
                            for jj = 2:nt
                                [A,B1] = parsysdata(pdGlocal(jj),ioChannel);
                                
                                X = Xaff(:,:,jj);
                                matAffX = [X*A + (X*A)' X*B1; B1'*X zeros(nw)];
                                lmiName = sprintf('Hinf: aff-X (term: %d:%d)',ll,jj);
                                lmis = [lmis, (matAffX + mu(jj-1)*eye(size(matAffX)) >= 0):lmiName];
                                
                            end
                            lmiName = 'Hinf : aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                            
                        case {'affY','affdY'}
                            mu = sdpvar(1,nt-1,'full');
                            for jj = 2:nt
                                [A,~,~,C1] = parsysdata(pdGlocal(jj),ioChannel);
                                
                                Y = Yaff(:,:,jj);
                                matAffY = [A*Y + (A*Y)' Y*C1'; C1*Y zeros(nz)];
                                lmiName = sprintf('Hinf: aff-Y (term: %d:%d)',ll,jj);
                                lmis = [lmis, (matAffY + mu(jj-1)*eye(size(matAffY)) >= 0):lmiName];
                            end
                            lmiName = 'Hinf: aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                    end
                    
                    for ii = simplex
                        
                        % one point of PARSET.POINTS
                        if strcmp(lyapType(1:3),'aff')
                            par = pdGlocal.parset.points(:,ii==simplex);
                            Mu = mu*(par.^2);
                        else
                            Mu = 0;
                        end
                        
                        % system matrices at par
                        [A,B1,b2,C1,c2,D11,d12,d21] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
                        
                        % controller vars at par
                        [Aki,Bki,Cki,Dki,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
                        
                        % Lyapunov matrices at par
                        [Y,X,dY,dX] = evalLyapFcn(lyapSet,ii,simplex);
                        
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
                            
                            % Hinf constraints
                            matXY = blkvar;
                            matXY(1,1) = (-dy/2 + A*Y + b2*Cki) + (-dy/2 + A*Y + b2*Cki)';
                            matXY(2,1) = (Aki + (A + b2*Dki*c2)');
                            matXY(2,2) = (dx/2 + X*A + Bki*c2) + (dx/2 + X*A + Bki*c2)';
                            matXY(3,1) = (B1 + b2*Dki*d21)';
                            matXY(3,2) = (X*B1 + Bki*d21)';
                            matXY(4,1) = C1*Y + d12*Cki;
                            matXY(4,2) = C1 + d12*Dki*c2;
                            matXY(4,3) = D11 + d12*Dki*d21;
                            if (ctrlSet.extFcn == 1)
                                % this option produce smaller norm for X
                                % and Y
                                matXY(3,3) = -g2*eye(nw);
                                matXY(4,4) = -eye(nz);
                            else
                                matXY(3,3) = -eye(nw);
                                matXY(4,4) = -g2*eye(nz);
                            end
                            % extra elements for imposing controller dependence
                            % using relaxation
                            na = size(E,2);
                            if (ctrlSet.extFcn == 1) && (na > 0) && (ctrlSet.XArYbnd == 0)
                                matXY(1,5) = Y*F;
                                matXY(2,5) = X*E;
                                matXY(5,5) = -eye(na);
                            end
                            matXY = sdpvar(matXY);
                            lmiName = sprintf('Hinf: Basic:BRL @(%d)',ii);
                            if isnumeric(matXY)
                                aux.lmiValue = all(eig(auxMat) <= 0);
                                aux.lmiName = lmiName;
                                lmisChecked = [lmisChecked, aux];
                            else
                                lmis = [lmis, (matXY <= -(eigtol + Mu)*eye(size(matXY))):lmiName];
                            end
                            
                        end
                        
                        
                    end
                    
                end
                
            otherwise
                error('LPVSYN:constHinf2LMI:LMI',...
                    'Projection type invalid')
        end
        
    case 'state'
        
        track('Hinf: State feedback',options.debug)
        
        for ii = 1:nv
            
            % system matrices at par
            [A,B1,b2,C1,~,D11,d12] = parsysdata(synSet.LPVatPoints(:,:,ii),ioChannel);
            
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
                
                % Hinf constraints
                matY = blkvar;
                matY(1,1) = dy + (A*Y + b2*Wi) + (A*Y + b2*Wi)';
                matY(1,2) = B1;
                matY(1,3) = (C1*Y + d12*Wi)';
                matY(2,2) = -eye(nw);
                matY(2,3) = D11';
                matY(3,3) = -g2*eye(nz);
                matY = sdpvar(matY);
                lmiName = sprintf('Hinf: BRL @(%d) - SF',ii);
                if isnumeric(matY)
                    aux.lmiValue = all(eig(matY) <= -eigtol);
                    aux.lmiName = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmis = [lmis, (matY <= -eigtol*eye(size(matY))):lmiName];
                end
            end
        end
        
    otherwise
        error('LPVSYN:constHinf2LMI:LMI',...
            'Type of feedback invalid')

end

% objective function
if isa(constraint.bound,'sdpvar')
    % soft constraint
    obj = obj + constraint.bound*constraint.factor;
end
    

% =======================================================
% Local functions

function track(msg,dbg)

if dbg
    disp(msg)
end

