 function lmis = constPoles2LMI_a(pdG,lyapSet,constraint,options)

% *** for internal use ***
%
% CONSTPOLES2LMI_A returns a set of LMIs for poles location constraint
%
% Use: 
%   [lmis,obj] = CONSTH22LMI_A(LPVatPoints,lyapSet,constraint,optionss)

% fbianchi - 2021-07-05
% fbianchi - 2024-01-26 - rev

      
% optimization settings
eigtol  = options.eigtol;

% number of systems in the description
nv = size(lyapSet.LPVatPoints,3);

% LMIs
lmis = [];
    
% Terms (one for each grid point)
% extra constraints for affine cases
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
    
    mu = [];
    for kk = 1:length(constraint.region.L)
        % iteration by constraint
        
        M = constraint.region.M{kk};
        
        if any(any(M)) && strcmp(lyapSet.type(1:3),'aff')
            mu_k = sdpvar(1,nt-1,'full');
            mu = cat(3,mu,mu_k);
            for jj = 2:nt
                % iteration by term
                A = pdGlocal(jj).A;
                X = Xaff(:,:,jj);
                
                if lyapSet.inv
                    % using inverse Lyapunov function
                    MAclXcl = kron(M,A*X) + kron(M',(A*X)');
                else
                    MAclXcl = kron(M,X*A) + kron(M',(X*A)');
                end
                
                lmiName = sprintf('Poles: aff-XY (Region:%d, term: %d)',kk,jj);
                lmis = [lmis, (MAclXcl + mu_k(jj-1)*eye(size(MAclXcl)) >= 0):lmiName];
            end
            lmiName = 'Poles: aff mu > 0';
            lmis = [lmis, (mu_k >= 0):lmiName];
            
        end
    end
    
    % Terms (one for each grid point)
    for ii = simplex
        
        % system matrices at par
        A = ssdata(lyapSet.LPVatPoints(:,:,ii));
        
        % Lyapunov matrices at par
        X = evalLyapFcn_a(lyapSet,ii,simplex);
        
        % Pole placement constraints
        for kk = 1:length(constraint.region.L)
            
            L = constraint.region.L{kk};
            M = constraint.region.M{kk};
            
            if any(any(L))
                LX = kron(L,X);
            else
                LX = 0;
            end
            if any(any(M))
                if lyapSet.inv
                    % using inverse Lyapunov function
                    MAclXcl = kron(M,A*X) + kron(M',(A*X)');
                else
                    MAclXcl = kron(M,X*A) + kron(M',(X*A)');
                end
            else
                MAclXcl = 0;
            end
            PoleMat = LX + MAclXcl;
            
            lmiName = sprintf('Poles: Region:%d, @(%d)',kk,ii);
            if strcmp(lyapSet.type(1:3),'aff')
                par = pdGlocal.parset.points(:,ii==simplex);
                Mu   = mu(:,:,kk)*(par.^2);
                lmis = [lmis, (PoleMat <= -(eigtol + Mu)*eye(size(PoleMat))):lmiName];
            else
                lmis = [lmis, (PoleMat <= -eigtol*eye(size(PoleMat))):lmiName];
            end
            
        end
    end
end

