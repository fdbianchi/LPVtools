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
            [~,b2,c2] = ssdata(pdG('meas','ctrl',1));
            
            % controller parameter dependence
%             if isa(pdG,'pass')
%                 % Fixme: ???
% %                 parset = ctrlSet.parset(ctrlSet.idxfcn(2:end)-1);
%                 parset = ctrlSet.parset;
%                
%             elseif isa(pdG,'pgss')
%                 parset = ctrlSet.parset;
%                 % FixMe: there are some bug, just a dirt & quick solution
%                 %parfcn = subfunc(ctrlSet.parfcn,ctrlSet.idxfcn(2:end)-1);
%                 if (nsys(pdG) ~= ctrlSet.nk)
%                     parfcn = subfunc(ctrlSet.parfcn,ctrlSet.idxfcn(2:end)-1,nsys(pdG)-1);
%                 else
%                     parfcn = ctrlSet.parfcn;
%                 end
%             else
%                 parset = pdG.parset;
%             end
                
            if strcmp(lyapSet.type,'cte')
                % controller when X and Y are constant

                [Y,X] = evalLyapFcn(lyapSet);
                
                % Compute N & M
                [U,S,V] = svd(eye(ns) - X*Y);
%                 fprintf('\nMin svd(I-XY) %d\n',min(diag(S)));
                Si  = diag(1./diag(S))^(1/2);
                Ni  = Si*U';
                Mti = V*Si;
                
                ak = zeros(size(ctrlSet.A));
                bk = zeros(size(ctrlSet.B));
                ck = zeros(size(ctrlSet.C));
                dk = zeros(size(ctrlSet.D));
                nk = size(ak,3);
                for ii = 1:nk
                    % plant data for computation of the controller
                    A = pdG.A(:,:,ctrlSet.idxfcn(ii));
                    
                    % controller variables
                    Aki = ctrlSet.A(:,:,ii);
                    Bki = ctrlSet.B(:,:,ii);
                    Cki = ctrlSet.C(:,:,ii);
                    Dki = ctrlSet.D(:,:,ii);
                    
                    % Undo the change of variable
                    ak(:,:,ii) = Ni*(Aki - X*(A-b2*Dki*c2)*Y - Bki*c2*Y - X*b2*Cki)*Mti;
                    bk(:,:,ii) = Ni*(Bki - X*b2*Dki);
                    ck(:,:,ii) = (Cki - Dki*c2*Y)*Mti;
                    dk(:,:,ii) = Dki;

                end
                
                % controller model
                if (ctrlSet.extFcn)%(nsys(pdG) ~= ctrlSet.nk)
                    ctrlSet.parfcn = subfunc(ctrlSet.parfcn,ctrlSet.idxfcn(2:end)-1,nsys(pdG)-1);
                end
                switch ctrlSet.class%synSet.pdGclass
                    case 'pass'
                        pdK = pass(ak,bk,ck,dk,ctrlSet.parset,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'pgss'
                        
                        pdK = pgss(ak,bk,ck,dk,ctrlSet.parset,...
                            ctrlSet.parfcn,...
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
                            ctrlSet.C, ctrlSet.D,parset,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'aff'
                        
                        if isa(ctrlSet.parfcn,'function_handle')
                            pdKhat = pgss(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, parset, parfcn,...
                                'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                                'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        else
                            pdKhat = pass(ctrlSet.A, ctrlSet.B,...
                                ctrlSet.C, ctrlSet.D, parset,...
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



