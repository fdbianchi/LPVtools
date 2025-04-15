function [lmis,obj] = constHinf2LMI_a(pdG,lyapSet,constraint,options)

% *** for internal use ***
%
% CONSTHINF2LMI_A returns a set of LMIs for Hinfinity performance constraint
%
% Use: 
%   [lmis,obj] = CONSTHINF2LMI_A(LPVatPoints,lyapSet,constraint,optionss)

% fbianchi - 2021-07-05
% fbianchi - 2024-01-26 - rev


% objective initial value
obj = 0;

% optimization settings
eigtol = options.eigtol;

% number of systems in the description
nv = size(lyapSet.LPVatPoints,3);

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
    % iteration by simplex
    
    simplex = Simplices(ll, :);
    
    if ~strcmp(lyapType,'cte')
        if strcmp(lyapSet.type(1:3), 'pwa') && (nv > 1)
            % PWA: the affine constraints are based on local
            % model for each simplex
            pdGlocal = ppss2affloc(pdG, ll);
            [~,~,Xaff] = evalLyapFcn_a(lyapSet, 1, simplex);
        else
            pdGlocal = pdG;
            Xaff = lyapSet.X;
        end
        nt = nsys(pdGlocal);
    end
    
    % extra constraints for affine/pwa cases
    if strcmp(lyapType(1:3), 'aff')
        
        mu = sdpvar(1,nt-1,'full');
        
        for jj = 2:nt
            % iteration by term in pdG
            
            if lyapSet.inv
                % using inverse Lyapunov function
                [A,~,C] = ssdata(pdGlocal(ioChannel.perf,ioChannel.dist,jj));
                
                X = Xaff(:,:,jj);
                matAff = [A*X + (A*X)' X*C'; C*X zeros(nz)];
                lmiName = sprintf('Hinf: affX-i (term: %d:%d)',ll,jj);
                lmis = [lmis, (matAff + mu(jj-1)*eye(size(matAff)) >= 0):lmiName];
                
            else
                [A,B] = ssdata(pdGlocal(ioChannel.perf,ioChannel.dist,jj));
                
                X = Xaff(:,:,jj);
                matAff = [X*A + (X*A)' X*B; B'*X zeros(nw)];
                lmiName = sprintf('Hinf: affX (term: %d:%d)',ll,jj);
                lmis = [lmis, (matAff + mu(jj-1)*eye(size(matAff)) >= 0):lmiName];
            end
        end
        lmiName = 'Hinf : aff mu > 0';
        lmis = [lmis, (mu >= 0):lmiName];
        
    end
    
    for ii = simplex
        % iteration by points in the simplex
        
        % one point of PARSET.POINTS
        if strcmp(lyapType(1:3),'aff')
            % par = pdG.parset.points(:,ii);
            par = pdGlocal.parset.points(:,ii==simplex);
            Mu  = mu*(par.^2);
        else
            Mu = 0;
        end
        
        % system matrices at par
        [A,B,C,D] = ssdata(lyapSet.LPVatPoints(ioChannel.perf,ioChannel.dist,ii));
        
        % Lyapunov matrices at par
        [X,dX] = evalLyapFcn_a(lyapSet,ii,simplex);
        
        % Bounded Real Lemma
        for jj = 1:length(dX)
            
            if lyapSet.inv
                % using inverse Lyapunov function
                BRL = blkvar;
                BRL(1,1) = (-0.5*dX{jj} + A*X) + (-0.5*dX{jj} + A*X)';
                BRL(1,2) = X*C';
                BRL(2,2) = -g2*eye(nz);
                BRL(1,3) = B;
                BRL(2,3) = D;
                BRL(3,3) = -eye(nw);
                BRL = sdpvar(BRL);
                
                lmiName = sprintf('Hinf: BRL-i @(%d)',ii);
                
            else
                BRL = blkvar;
                BRL(1,1) = (0.5*dX{jj} + X*A) + (0.5*dX{jj} + X*A)';
                BRL(1,2) = X*B;
                BRL(2,2) = -eye(nw);
                BRL(1,3) = C';
                BRL(2,3) = D';
                BRL(3,3) = -g2*eye(nz);
                BRL = sdpvar(BRL);
                
                lmiName = sprintf('Hinf: BRL @(%d)',ii);
                
            end
            lmis = [lmis, (BRL <= -eigtol*eye(size(BRL))):lmiName];
        end
        
    end
    
end

% objective function
if isa(constraint.bound,'sdpvar')
    % soft constraint
    obj = obj + mean(constraint.bound)*constraint.factor;
end
    
