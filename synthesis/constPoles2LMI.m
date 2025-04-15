function [lmis,lmisChecked] = constPoles2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% *** for internal use ***
%
% CONSTPOLES2LMI returns a set of LMIs for Pole placement constraint
%
% Use: 
%   lmis = constPoles2LMI(pdG,synSet,lyapSet,ctrlSet,constraint,options)

% fbianchi - 2020-04-24
      
% optimization settings
eigtol = options.eigtol;

% number of systems in the description
nv = size(synSet.LPVatPoints,3);

% system order
ns = order(synSet.LPVatPoints(:,:,1));

% % number of terms
% nt = nsys(pdG);

% LMIs
lmisChecked = [];
lmis = [];
switch ctrlSet.fback         % type of feedback
    
    case 'output'

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
            
            mu = [];
            for kk = 1:length(constraint.region.L)
                
                M = constraint.region.M{kk};
                
                if any(any(M))
                    switch lyapType
                        
                        case {'aff','affdXY'}
                            mu_k = sdpvar(1,nt-1,'full');
                            mu = cat(3,mu,mu_k);
                            for jj = 2:nt
                                A = pdGlocal(jj).A;
                                Y = Yaff(:,:,jj);
                                X = Xaff(:,:,jj);
                                
                                AclXcl = blkdiag(A*Y,X*A);
                                Myx = kron(M,AclXcl) + kron(M',AclXcl');
                                
                                lmiName = sprintf('Poles: aff-XY (Region:%d, term: %d)',kk,jj);
                                lmis = [lmis, (Myx + mu_k(jj-1)*eye(size(Myx)) >= 0):lmiName];
                            end
                            lmiName = 'Poles: aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                            
                        case {'affX','affdX'}
                            mu_k = sdpvar(1,nt-1,'full');
                            mu = cat(3,mu,mu_k);
                            for jj = 2:nt
                                A = pdGlocal(jj).A;
                                X = Xaff(:,:,jj);
                                
                                Mx = kron(M,X*A) + kron(M',(X*A)');
                                lmiName = sprintf('Poles: aff-X (Region:%d, term: %d)',kk,jj);
                                lmis = [lmis, (Mx + mu_k(jj-1)*eye(size(Mx)) >= 0):lmiName];
                                
                            end
                            lmiName = 'Poles: aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                            
                        case {'affY','affdY'}
                            mu_k = sdpvar(1,nt-1,'full');
                            mu = cat(3,mu,mu_k);
                            for jj = 2:nt
                                A = pdGlocal(jj).A;
                                Y = Yaff(:,:,jj);
                                
                                My = kron(M,A*Y) + kron(M',(A*Y)');
                                lmiName = sprintf('Poles: aff-Y (Region:%d, term: %d)',kk,jj);
                                lmis = [lmis, (My + mu_k(jj-1)*eye(size(My)) >= 0):lmiName];
                            end
                            lmiName = 'Poles: aff mu > 0';
                            lmis = [lmis, (mu >= 0):lmiName];
                    end
                end
            end
            
            % Terms (one for each grid point)
            for ii = simplex
                
        
                % system matrices at par
%                 A = ssdata(synSet.LPVatPoints(:,:,ii));
                [A,b2,c2] = ssdata(synSet.LPVatPoints('meas','ctrl',ii));
                
                % controller vars at par
                [Aki,Bki,Cki,Dki,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
                
                % Lyapunov matrices at par
                [Y,X,dY,dX] = evalLyapFcn(lyapSet,ii,simplex);
                
                % Pole placement constraints
                for kk = 1:length(constraint.region.L)
                    
                    L = constraint.region.L{kk};
                    M = constraint.region.M{kk};
                    
                    if any(any(L))
                        LX = kron(L,[Y eye(ns);eye(ns) X]);
                    else
                        LX = 0;
                    end
                    if any(any(M))
                        AclXcl = [A*Y+b2*Cki,  A+b2*Dki*c2;
                                  Aki,         X*A+Bki*c2];
                        MAclXcl = kron(M,AclXcl) + kron(M',AclXcl');
                    else
                        MAclXcl = 0;
                    end
                    
                    % ensuring a given controller parameter dependence by
                    % relaxation
                    na = size(E,2);
                    if (ctrlSet.extFcn == 1) && (na > 0) && any(any(M)) && (ctrlSet.XArYbnd == 0)
                        It  = [zeros(na,2*na);eye(na) zeros(na)];
                        T   = kron(M,It) + kron(M.',It') + eye(2*na*size(M));
                        [U,S] = svd(T); 
                        T12 = U*sqrt(S);
                        extra   = kron(eye(size(M)),blkdiag(Y*F,X*E))*T12;
                        PoleMat = [LX + MAclXcl,  extra;...
                                   extra',       -eye(2*na*size(M))];
                    else
                        PoleMat = LX + MAclXcl;
                    end
                    
                    lmiName = sprintf('Poles: Region:%d, @(%d)',kk,ii);
                    if isnumeric(PoleMat)
                        aux.lmiValue = all(eig(PoleMat) <= 0);
                        aux.lmiName = lmiName;
                        lmisChecked = [lmisChecked, aux];
                        
                    else
                        if strcmp(lyapSet.type(1:3),'aff')
                            % par  = pdG.parset.points(:,ii);
                            par = pdGlocal.parset.points(:,ii==simplex);
                            Mu   = mu(:,:,kk)*(par.^2);
                            lmis = [lmis, (PoleMat <= -(eigtol + Mu)*eye(size(PoleMat))):lmiName];
                        else
                            lmis = [lmis, (PoleMat <= -eigtol*eye(size(PoleMat))):lmiName];
                        end
                    end
                    
                end
            end
        end
        
    case 'state'
        
            [~,b2] = ssdata(synSet.LPVatPoints(:,'ctrl',1));

            % extra constraints for affine cases
            mu = [];
            for kk = 1:length(constraint.region.L)
                
                M = constraint.region.M{kk};
                
                if any(any(M))
                    switch lyapSet.type
                        case 'pwa'
                            
                        case {'aff','affY','affdY'}
                            mu_k = sdpvar(1,nt-1,'full');
                            mu = cat(3,mu,mu_k);
                            for jj = 2:nsys(pdG)
                                A = pdG(jj).A;
                                Y = lyapSet.Y(:,:,jj);
                                My = kron(M,A*Y) + kron(M',(A*Y)');
                                lmiName = sprintf('Poles: aff-Y (Region:%d, term: %d)',kk,jj);
                                lmis = [lmis, (My + mu_k(jj-1)*eye(size(My)) >= 0):lmiName];
                            end
                    end
                end
            end
            lmis = [lmis, mu >= 0];
            
            % Terms (one for each grid point)
            for ii = 1:nv
                
                % system matrices at par
                A = ssdata(synSet.LPVatPoints(:,:,ii));
                
                % controller vars at par
                [~,~,~,Wi] = evalCtrlVars(ctrlSet,ii);%par,ii);
                
                % Lyapunov matrices at par
                Y = evalLyapFcn(lyapSet,ii);%par);
                
                % ----------------------
                % Pole placement constraints
                for kk = 1:length(constraint.region.L)
                    
                    L = constraint.region.L{kk};
                    M = constraint.region.M{kk};
                    
                    if any(any(L))
                        LX = kron(L,Y);
                    else 
                        LX = 0;
                    end
                    if any(any(M))
                        AclXcl = A*Y + b2*Wi;
                        MAclXcl = kron(M,AclXcl) + kron(M',AclXcl');
                    else
                        MAclXcl = 0;
                    end
                    PoleMat = LX + MAclXcl;
                    
                    lmiName = sprintf('Poles: Region:%d, @(%d)',kk,ii);
                    if strcmp(lyapSet.type(1:3),'aff')
                        par = pdG.parset.points(:,ii);
                        Mu = mu(:,:,kk)*(par.^2);
                        lmis = [lmis, (PoleMat <= -(eigtol + Mu)*eye(size(PoleMat))):lmiName];
                    else
                        lmis = [lmis, (PoleMat <= -eigtol*eye(size(PoleMat))):lmiName];
                    end
                   
                end
            end
            
    otherwise
        error('LPVSYN:constPoles2LMI:LMI',...
            'Type of feedback invalid')

end



