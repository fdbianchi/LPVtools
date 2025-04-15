function pdK = lpvctrl(pdG,synSet,lyapSet,ctrlSet,constraint)

% *** for internal use ***
%
% LPVCTRL computes a controller satisfying a performance constraints
%
% Use: 
%   pdK = lpvctrl(pdG,synSet,lyapSet,ctrlSet,constraint)
%
% Inputs:
%   - synSet:       structure with synthesis information (see synlpv)
%   - pdG:          augmented plant (pxss object)
%   - lyapSet:      structure with Lyapunov fcn information
%   - ctrlSet:      structure with controller parameter information
%   - constraint:   constraints
%
% Outputs:
%   - pdK: pass, ppss, pgss or pcss dependending on pdG and lyapSet

% fbianchi - 2020-05-26

% number of systems in the description
nv = nsys(pdG);
% system order
ns = order(pdG);       

% controller reconstruction
if strcmp(ctrlSet.fback,'state')
    % =====================================================================
    % state feedback controller
    
    switch ctrlSet.sch
        
        case 'gs'
            % case state feedback gain scheduling
            if strcmp(lyapSet.type,'cte')
                % Y = cte
                W = double(ctrlSet.D);
                Y = evalLyapFcn(lyapSet);
                dk = zeros(size(W));
                for ii = 1:nv
                    dk(:,:,ii) = W(:,:,ii)/Y;
                end
                
                switch ctrlSet.class%synSet.pdGclass
                    case 'pass'
                        pdK = pass([],[],[],dk,ctrlSet.parset,...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'pgss'
                        pdK = pgss([],[],[],dk,ctrlSet.parset,ctrlSet.parfcn,...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'ppss'
                        pdK = ppss([],[],[],dk,ctrlSet.parset,...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case {'ss','zpk','tf'}
                        pdK = ss([],[],[],dk,...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'psys-aff'
                        pdKaux = pass([],[],[],dk,ctrlSet.parset);
                        pdK = psys(pdKaux);
                        
                    case 'psys-pol'
                        pdKaux = ppss([],[],[],dk,ctrlSet.parset);
                        pdK = psys(pdKaux);

                    otherwise
                        error('LPVCTRL:inputError','Invalid pdG class')
                        
                end
                
            else
                % when Y is not cte, the controller must be computed online
                switch ctrlSet.interp
                    case 'pwa'
                        pdKhat = ppss(ctrlSet.A, ctrlSet.B,...
                            ctrlSet.C, ctrlSet.D,ctrlSet.parset,...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'aff'
                        
                        if isa(ctrlSet.parfcn,'function_handle')
                            pdKhat = pgss(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, ctrlSet.parset, ctrlSet.parfcn,...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        else
                            pdKhat = pass(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, ctrlSet.parset,...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        end
                        
                    otherwise
                        error('LPVCTRL:inputError','Invalid interpolation method')
                        
                end                 
                
                if isa(lyapSet.parfcnY,'function_handle')
                    pdY = pgss([],[],[],lyapSet.Y,...
                        pdG.parset,lyapSet.parfcnY);
                else
                    pdY = lyapSet.Y;
                end
                pdK = pcss(pdKhat, [], pdY);
               
            end
            
        case 'robust'
            % case state feedback robust
            if strcmp(lyapSet.type,'cte')
                % case state feedback gain scheduling
                W  = double(ctrlSet.D);
                Y  = evalLyapFcn(lyapSet);
                pdK = ss(W/Y,'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                
            else
                % when Y is not cte, the controller must be computed
                % online
                error('LPVCTRL:inputError',...
                    'when the controller is robust, Y must be constant')
            end
    end
    
else
    % =====================================================================
    % output feedback case
    switch synSet.type
        
        case 'basic'
            % -------------------------------------------------------------
            % controller change of variable
            %

            % parameter functions according to complexity reduction
            if (ctrlSet.extFcn) && strcmp(ctrlSet.class,'pgss')
                ctrlSet.parfcn = subfunc(ctrlSet.parfcn,ctrlSet.idxfcn(2:end)-1,nsys(pdG)-1);
            end
                
            % undoing the change of variables
            % Case: X and Y constant
            if strcmp(lyapSet.type,'cte')
                
                % closed-loop Lyapunov matrix reconstruction
                [Y,X] = evalLyapFcn(lyapSet);
                % Compute N & M
                [U,S,V] = svd(eye(ns) - X*Y);
                Si  = diag(1./diag(S))^(1/2);
                Ni  = Si*U';
                Mti = V*Si;
                
                % auxiliary controller variables
                akh = ctrlSet.A;
                bkh = ctrlSet.B;
                ckh = ctrlSet.C;
                dk = ctrlSet.D;
                nk = size(akh,3);
                ns = size(akh,1);
                [nu,ny] = size(dk(:,:,1));
                
                % controller variables
                ak = zeros(ns,ns,nk);
                bk = zeros(ns,ny,nk);
                ck = zeros(nu,ns,nk);
                
                % plant
                pdG22 = pdG('meas','ctrl');
                A = pdG22.A;
                B2 = pdG22.B;
                C2 = pdG22.C;
                
                if all(synSet.pdMeasCtrl == [0 0])
                    % b2 and c2 constant
                    
                    for ii = 1:nk
                        
                        a = A(:,:,ctrlSet.idxfcn(ii));
                        b2 = B2(:,:,1);
                        c2 = C2(:,:,1);
                        aki = akh(:,:,ii);
                        bki = bkh(:,:,ii);
                        cki = ckh(:,:,ii);
                        dki = dk(:,:,ii);
                        % Undo the change of variable
                        ak(:,:,ii) = Ni*(aki - X*(a + b2*dki*c2)*Y - bki*c2*Y - X*b2*cki)*Mti;
                        bk(:,:,ii) = Ni*(bkh(:,:,ii) - X*b2*dki);
                        ck(:,:,ii) = (cki - dki*c2*Y)*Mti;
                        
                    end
                    
                elseif (all(synSet.pdMeasCtrl == [1 0]) && (ctrlSet.pdOut == 0))
                    % c2, ck and dk constant
                    
                    for ii = 1:nk
                        
                        a = A(:,:,ctrlSet.idxfcn(ii));
                        b2 = B2(:,:,ii);
                        c2 = C2(:,:,1);
                        aki = akh(:,:,ii);
                        bki = bkh(:,:,ii);
                        cki = ckh(:,:,1);
                        dki = dk(:,:,1);
                        % Undo the change of variable
                        ak(:,:,ii) = Ni*(aki - X*(a + b2*dki*c2)*Y - bki*c2*Y - X*b2*cki)*Mti;
                        bk(:,:,ii) = Ni*(bkh(:,:,ii) - X*b2*dki);
                        if (ii > 1)
                            ck(:,:,ii) = (cki - dki*c2*Y)*Mti;
                        end
                        
                    end
                    
                elseif (all(synSet.pdMeasCtrl == [0 1]) && (ctrlSet.pdIn == 0))
                    % b2, bk and dk constant
                    
                    for ii = 1:nk
                        
                        a = A(:,:,ctrlSet.idxfcn(ii));
                        b2 = B2(:,:,1);
                        c2 = C2(:,:,ii);
                        aki = akh(:,:,ii);
                        bki = bkh(:,:,1);
                        cki = ckh(:,:,ii);
                        dki = dk(:,:,1);
                        % Undo the change of variable
                        ak(:,:,ii) = Ni*(aki - X*(a + b2*dki*c2)*Y - bki*c2*Y - X*b2*cki)*Mti;
                        if (ii > 1)
                            bk(:,:,ii) = Ni*(bkh(:,:,ii) - X*b2*dki);
                        end
                        ck(:,:,ii) = (cki - dki*c2*Y)*Mti;
                        
                    end
                    
                elseif strcmp(ctrlSet.class,'pgss')
                    % the rest of cases lead to multi-affine controller,
                    % therefore only for general plants
                    
                    % expanding the parameter functions with the cross
                    % product terms (e.g. b2(p)*ck(p))
                    if ((synSet.pdMeasCtrl(1) == 1) && (ctrlSet.pdOut == 1)) || ...
                            ((synSet.pdMeasCtrl(2) == 1) && (ctrlSet.pdIn == 1))
                        nt = 1 + nk^2;
                        ctrlSet.parfcn =@(p) [ctrlSet.parfcn(p);reshape(ctrlSet.parfcn(p)*ctrlSet.parfcn(p)',(nk-1)^2,1)];
                    else
                        nt = nk;
                    end
                    % controller variables
                    ak = zeros(ns,ns,nt);
                    bk = zeros(ns,ny,nt);
                    ck = zeros(nu,ns,nt);                    
                    
                    if ((synSet.pdMeasCtrl(2) == 1) && (ctrlSet.pdIn == 1))
                        % c2, bk and dk parameter dependent
                    
                    elseif ((synSet.pdMeasCtrl(1) == 1) && (ctrlSet.pdOut == 1))
                        % c2, bk and dk parameter dependent
                    
                    
                    
                    % constant and linear terms
                    for ii = 1:nk
                        
                        a = A(:,:,ctrlSet.idxfcn(ii));
                        b2 = B2(:,:,ii);
                        c2 = C2(:,:,ii);
                        aki = akh(:,:,ii);
                        bki = bkh(:,:,1);
                        cki = ckh(:,:,1);
                        dki = dk(:,:,1);
                        % Undo the change of variable
                        ak(:,:,ii) = Ni*(aki - X*(a + b2*dki*c2)*Y - bki*c2*Y - X*b2*cki)*Mti;
                        bk(:,:,ii) = Ni*(bkh(:,:,ii) - X*b2*dki);
                        ck(:,:,ii) = (cki - dki*c2*Y)*Mti;
                        
                    end
                   
                elseif strcmp(ctrlSet.class,'pgss') && (ctrlSet.dk == 0)
                    

                    
                    for ii = 1:nk

                        a = A(:,:,ctrlSet.idxfcn(ii));
                        b2 = B2(:,:,ii);
                        c2 = C2(:,:,ii);
                        aki = akh(:,:,ii);
                        bki = bkh(:,:,1);
                        cki = ckh(:,:,1);
                        % Undo the change of variable
                        ak(:,:,ii) = Ni*(aki - X*a*Y - bki*c2*Y - X*b2*cki)*Mti;
                        bk(:,:,ii) = Ni*bki;
                        ck(:,:,ii) = cki*Mti;
                        
                    end
                    
                    if ((synSet.pdMeasCtrl(1) == 1) && (ctrlSet.pdOut == 1))
                        
                        % b2 and ck parameter dependent
                        
                        for jj = 1:nk   % extra linear terms
                            ak(:,:,ii) = -Ni*X*b2(:,:,jj)*ckh(:,:,1)*Mti;
                        end
                        
                        idx = nchoosek(2:nk,2);
                        for jj = 1:size(idx,1)   % multi-linear terms
                            [~,b2] = ssdata(pdG('meas','ctrl',idx(jj,1)));
                            ak(:,:,ii) = - Ni*X*(b2(:,:,idx(jj,1))*ckh(:,:,idx(jj,2)) + ...
                                b2(:,:,idx(jj,2))*ckh(:,:,idx(jj,1)))*Mti;
                            ii = ii + 1;
                        end
                        
                        for jj = 2:nk   % quadratic terms
                            ak(:,:,ii) = -Ni*X*b2(:,:,jj)*ckh(:,:,jj)*Mti;
                            ii = ii + 1;
                        end
                        
                    end
                    
                    if ((synSet.pdMeasCtrl(2) == 1) && (ctrlSet.pdIn == 1))
                        
                        % c2 and bk parameter dependent
                        
                        for jj = 1:nk   % extra linear terms
                            ak(:,:,ii) = -Ni*bkh(:,:,1)*c2(:,:,jj)*Y*Mti;
                        end
                        
                        idx = nchoosek(2:nk,2);
                        for jj = 1:size(idx,1)   % multi-linear terms
                            ak(:,:,ii) = -Ni*(bkh(:,:,idx(jj,1))*c2(:,:,idx(jj,2)) +...
                                bkh(:,:,idx(jj,2))*c2(:,:,idx(jj,1)))*Y*Mti;
                            ii = ii + 1;
                        end
                        
                        for jj = 2:nk   % quadratic terms
                            ak(:,:,ii) = -Ni*bkh(:,:,jj)*c2(:,:,jj)*Y*Mti;
                            ii = ii + 1;
                        end
                    end
                    
                    
                else
                    error('LPVCTRL:inputError','It is not possible to compute the controller')
                end
                    
                    % controller model
                    switch ctrlSet.class%synSet.pdGclass
                        case 'pass'
                            pdK = pass(ak,bk,ck,dk,ctrlSet.parset,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                            
                        case 'pgss'
                            
                            pdK = pgss(ak,bk,ck,dk,ctrlSet.parset,parfcn,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                            
                        case 'ppss'
                            pdK = ppss(ak,bk,ck,dk,ctrlSet.parset,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                            
                        case {'ss','zpk','tf'}
                            pdK = ss(ak,bk,ck,dk,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                            
                        case 'psys-aff'
                            pdKaux = pass(ak,bk,ck,dk,ctrlSet.parset);
                            pdK = psys(pdKaux);
                            
                        case 'psys-pol'
                            pdKaux = ppss(ak,bk,ck,dk,ctrlSet.parset);
                            pdK = psys(pdKaux);
                            
                        otherwise
                            error('LPVCTRL:inputError','Invalid pdG class')
                            
                    end
                    
                
            else
                
                % when X or Y are not cte, the controller must be computed
                % online
                
                % auxiliary controller matrices
                switch ctrlSet.interp
                    
                    case 'pwa'
                        pdKhat = ppss(ctrlSet.A, ctrlSet.B,...
                            ctrlSet.C, ctrlSet.D,ctrlSet.parset,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'aff'
                        
                        if isa(ctrlSet.parfcn,'function_handle')
                            pdKhat = pgss(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, ctrlSet.parset, ctrlSet.parfcn,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        else
                            pdKhat = pass(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, ctrlSet.parset,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        end
                    otherwise
                        error('LPVCTRL:inputError','Invalid interpolation method')
                        
                end
                
                % plant matrices
                pdG22 = pdG('meas','ctrl');
                
                % Lyapunov functions
                if isa(lyapSet.parfcnY,'function_handle')
                    pdY = pgss([],[],[],lyapSet.Y,...
                        pdG.parset,lyapSet.parfcnY);

                elseif (size(lyapSet.Y,3) > 1)
                    pdY = ppss([],[],[],lyapSet.Y,...
                        pdG.parset);
                else
                    pdY = lyapSet.Y;
                end
                
                if isa(lyapSet.parfcnX,'function_handle')
                    pdX = pgss([],[],[],lyapSet.X,...
                        pdG.parset,lyapSet.parfcnX);

                elseif (size(lyapSet.X,3) > 1)
                    pdX = ppss([],[],[],lyapSet.X,...
                        pdG.parset);
                else
                    pdX = lyapSet.X;
                end
                
                % all info in PCSS object
                pdK = pcss(pdKhat, pdG22, pdY, pdX);
                
            end
            
%     case {'prAhinf','prAh2'}
%         % projected w.r.t. Ak
%         %
%         for ii=1:nv,
%             
%             % controller variables
%             Bki = double(vars.Bk(:,:,ii));
%             Cki = double(vars.Ck(:,:,ii));
%             dk  = double(vars.Dk(:,:,ii));
%             
%             % Compute Ak
%             [A,B1,b2,C1,c2,D11,d12,d21] = parsysdata(pdG(ii),ios);
% 
%             if strcmp(type(4:end),'hinf')
%                 gopt = double(bnd);
%                 delta = [-gopt*eye(ios.nw), (D11+d12*dk*d21)';D11+d12*dk*d21, -gopt*eye(ios.nz)];
%                 Aki = -(A+b2*dk*c2)' + ...
%                     [Xopt*B1+Bki*d21, (C1+d12*dk*c2)']/delta*[(B1+b2*dk*d21)';C1*Yopt+d12*Cki];
%             else
%                 Aki = -(A+b2*dk*c2)'-(Xopt*B1+Bki*d21)*(B1+b2*dk*d21)';
%             end
%             
%             % Undo the change of variable
%             ak = Ni*(Aki - Xopt*(A-b2*dk*c2)*Yopt - Bki*c2*Yopt - Xopt*b2*Cki)*Mti;
%             bk = Ni*(Bki - Xopt*b2*dk);
%             ck = (Cki - dk*c2*Yopt)*Mti;
%             
%             pdk(:,:,ii) = ss(ak,bk,ck,dk);
%             
%             
%         end
        
        case 'prFull'
            % ToDo: program my own function!
            % totally projected
            %
            g = constraint.bound;
%             g = sqrt(constraint.bound);
            orderK = -1;
            pdk = ss();
            % number of ctrl input and output
            ny = length(pdG.OutputGroup.meas);
            nu = length(pdG.InputGroup.ctrl);
            ioChannel.perf = constraint.outputs;
            ioChannel.dist = constraint.inputs;

            count = 0;
            while isempty(pdk)
                for ii = 1:nv
                    
                    [Y,X] = evalLyapFcn(lyapSet,ii);
                    
                    [a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(pdG(:,:,ii),ioChannel);
                    sys = ltisys(a,[b1 b2],[c1; c2],[d11 d12;d21 d22]);
                    
%                     [kv, gfin] = klmi(sys,[ny, nu],g,...
%                         Y,g*eye(ns),X,g*eye(ns),1e-4);
                    [kv, gfin] = klmi(sys,[ny, nu],g,...
                        g\Y,g*eye(ns),g*X,g*eye(ns),1e-4);
%                     [kv, gfin] = klmi(sys,[ny, nu],g,...
%                         Y,g*eye(ns),X,g\eye(ns),1e-4);
                    
                    nk = sinfo(kv);   % order of K
                    
                    if (~isempty(pdk) && nk ~= orderK)
                        pdk = ss();  g = 1.01*g;
                        break
                    else
                        orderK = nk;
                        pdk(:,:,ii) = mat2lti(kv);
                    end
                end
                count = count + 1;
            end
            
            pdK = ppss(pdk,pdG.parset,...
                        'OutputName',pdG.OutputName(pdG.OutputGroup.meas),...
                        'InputName',pdG.InputName(pdG.InputGroup.ctrl));
            
    end
end



