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
                
            % Case: X and Y constant
            if strcmp(lyapSet.type,'cte')
                
                % auxiliary controller variables
                akh = ctrlSet.A;
                bkh = ctrlSet.B;
                ckh = ctrlSet.C;
                dk = ctrlSet.D;
                nk = size(akh,3);
                [nu,ny] = size(dk(:,:,1));
                
                % condition to compute a controller
                if ((synSet.pdMeasCtrl(1) == 1) && (ctrlSet.pdOut == 1)) || ...
                   ((synSet.pdMeasCtrl(2) == 1) && (ctrlSet.pdIn == 1)) || ...
                   (all(synSet.pdMeasCtrl) && ctrlSet.dk > 0)
                    % these cases will be multi-affine
                    if ~strcmp(ctrlSet.class,'pgss')
                        error('LPVCTRL:inputError',...
                            'It is not possible to compute the controller.\nThe controller should be multi-affine and the plant is affine/pwa')
                    elseif all([synSet.pdMeasCtrl ctrlSet.pdIn ctrlSet.pdOut]) && ctrlSet.dk > 0
                            error('LPVCTRL:inputError',...
                                'It is not possible to compute a multi-affine controller')
                    end
                    % expanding the parameter functions with the cross
                    % product terms (e.g. b2(p)*ck(p))
                    nt = 2*nk - 1 + nchoosek(nk-1,2);
                    ctrlSet.parfcn = multfunc(ctrlSet.parfcn,ctrlSet.parset);
                else
                    nt = nk;
                end
                
                % closed-loop Lyapunov matrix reconstruction
                [Y,X] = evalLyapFcn(lyapSet);
                % Compute N & M
                [U,S,V] = svd(eye(ns) - X*Y);
                Si = diag(1./diag(S))^(1/2);
                Ni = Si*U';
                Mti = V*Si;
                
                % controller variables
                ak = zeros(ns,ns,nt);
                bk = zeros(ns,ny,nt);
                ck = zeros(nu,ns,nt);
                
                % plant
                pdG22 = pdG('meas','ctrl');
                A = pdG22.A;
                B2 = pdG22.B;
                C2 = pdG22.C;
                
                % matrix products
                dkc2 = matProd(dk,C2,[],ctrlSet.interp);
                b2dkc2 = matProd(B2,dkc2,nt,ctrlSet.interp);
                bkc2 = matProd(bkh,C2,nt,ctrlSet.interp);
                b2ck = matProd(B2,ckh,nt,ctrlSet.interp);
                dkc2 = matProd(dk,C2,nt,ctrlSet.interp);
                b2dk = matProd(B2,dk,nt,ctrlSet.interp);
               
                % Undo the change of variable
                for ii = 1:nk
                    % affine terms
                    a = A(:,:,ctrlSet.idxfcn(ii));
                    
                    ak(:,:,ii) = Ni*(akh(:,:,ii) - X*(a - b2dkc2(:,:,ii))*Y -...
                                     bkc2(:,:,ii)*Y - X*b2ck(:,:,ii))*Mti;
                    bk(:,:,ii) = Ni*(bkh(:,:,ii) - X*b2dk(:,:,ii));
                    ck(:,:,ii) = (ckh(:,:,ii) - dkc2(:,:,ii)*Y)*Mti;
                end
                
                for jj = ii+1:nt
                    % rest of terms in case of pgss plants
                    ak(:,:,jj) = Ni*(X*b2dkc2(:,:,jj)*Y - bkc2(:,:,jj)*Y - X*b2ck(:,:,jj))*Mti;
                    bk(:,:,jj) = -Ni*X*b2dk(:,:,jj);
                    ck(:,:,jj) = -dkc2(:,:,jj)*Y*Mti;
                end
                
                % controller model
                switch ctrlSet.class%synSet.pdGclass
                    case 'pass'
                        pdK = pass(ak,bk,ck,dk,ctrlSet.parset,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'pgss'
                        
                        pdK = pgss(ak,bk,ck,dk,ctrlSet.parset,ctrlSet.parfcn,...
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
            
                switch ctrlSet.class%synSet.pdGclass
                       
                    case 'ppss'
                        pdK = ppss(pdk.A,pdk.B,pdk.C,pdk.D,ctrlSet.parset,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case {'ss','zpk','tf'}
                        pdK = ss(pdk.A,pdk.B,pdk.C,pdk.D,...
                            'InputName',pdG.OutputName(pdG.OutputGroup.meas),...
                            'OutputName',pdG.InputName(pdG.InputGroup.ctrl));
                        
                    case 'psys-pol'
                        pdKaux = ppss(pdk.A,pdk.B,pdk.C,pdk.D,ctrlSet.parset);
                        pdK = psys(pdKaux);
                        
                    otherwise
                        error('LPVCTRL:inputError','Invalid pdG class')
                        
                end                
            
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions

% computing the new parfcn when there are cross products
function fout = multfunc(fcn,parset)

if ~isa(fcn,'function_handle')
    error('SUBFUNC:inputerror','fcn must be an anonymous function')
end
    
try
    nv = length(fcn(parset.points(:,1)));
catch err
    error('PGSS:PGSS:inputError',...
        'Function definition does not correspond with the parameter set')
end
    
% string with old function
strFcn = func2str(fcn);

% get the argument
idx1 = find(strFcn=='@',1,'first');
idx2 = find(strFcn==')',1,'first');
strCtrlFcn = strFcn(idx1:idx2);

% elements in oldFcn
idx1 = find(strFcn=='[',1,'first') + 1;
idx2 = find(strFcn==']',1,'first') - 1;
strElems = strsplit(strFcn(idx1:idx2),{';',',',' '});
ne = length(strElems);

if (ne ~= nv)
    
    % "linear" terms
    strFcn = ['@(f)[' sprintf('f(%d);',1:nv)];
    % "multi-affine" terms
    idx = nchoosek(1:nv,2);
    for ii = 1:size(idx,1)
        strFcn = [strFcn, sprintf('f(%d).*f(%d);',idx(ii,:))];
    end
    % "quadratic" terms
    for jj = 1:nv
        strFcn = [strFcn, sprintf('f(%d).^2;',jj)];
    end
    strFcn(end) = ']';
    f = str2func(strFcn);
    fout =@(p) f(fcn(p));
    
else
    
    % "linear" terms
    fi = strElems;
    
    % "multi-affine" terms
    idx = nchoosek(1:ne,2);
    for ii = 1:size(idx,1)
        fi{ii+ne} = [strElems{idx(ii,1)} '*' strElems{idx(ii,2)}];
    end
    
    % "quadratic" terms
    for jj = 1:ne
        ii = ne + ii + 1;
        fi{ii} = [strElems{jj} '^2'];
    end
    
    strNewFcn = [strCtrlFcn ' ['];
    for ii = 1:length(fi)
        strNewFcn = [strNewFcn fi{ii} ';'];
    end
    strNewFcn(end) = ']';
    
    fout = str2func(strNewFcn);
end

% ------------------------------------------------------------------------
% parameter dependent matrices products
function C = matProd(A,B,nt,int)

tol = 1e-6;
if all(A(:,:,2:end) == 0,'all')
    % A is constant (case affine)
    na = 1;
else
    % case pwa
    nv = size(A,3);
    na = 1;
    for ii = 2:nv
        if (norm(A(:,:,ii) - A(:,:,1),1) > tol)
            na = nv;
        end
    end
end

if all(B(:,:,2:end) == 0,'all')
    % B is constant (case affine)
    nb = 1;
else
    % case pwa
    nv = size(B,3);
    nb = 1;
    for ii = 2:nv
        if (norm(B(:,:,ii) - B(:,:,1),1) > tol)
            nb = nv;
        end
    end
end

if (na > 1) && (nb > 1) && (na ~= nb)
    error('matProd:inputError','the third dimension of A and B must be equal')
end
    
nr = size(A,1);     % # of rows
nc = size(B,2);     % # of columns
m = max(na,nb);     % # of individual terms

if (na == 1) && (nb == 1)
    % result: constant matrix
    C(:,:,1) = A(:,:,1)*B(:,:,1);
    if (nargin > 2) 
        if strcmp(int,'aff')
            C = cat(3,C,zeros(nr,nc,nt-1));
        else
            if isempty(nt)
                C = repmat(C,1,1,m);
            else
                C = repmat(C,1,1,nt);
            end
        end
    end
    
elseif (na == 1) || (nb == 1)
    % result: affine
    C = zeros(nr,nc,m);
    if (na == 1)
        for ii = 1:m
            C(:,:,ii) = A(:,:,1)*B(:,:,ii);
        end
    else
        for ii = 1:m
            C(:,:,ii) = A(:,:,ii)*B(:,:,1);
        end
    end        
    if (nargin > 2) 
        if strcmp(int,'aff')
            C = cat(3,C,zeros(nr,nc,nt-m));
        end
    end
else
    % result: multi-affine
    nt = 2*m - 1 + nchoosek(m-1,2);   % # of terms
    C = zeros(nr,nc,nt);
    % constant term
    C(:,:,1) = A(:,:,1)*B(:,:,1);
    % affine terms
    for ii = 2:m
        C(:,:,ii) = A(:,:,1)*B(:,:,ii) + A(:,:,ii)*B(:,:,1);
    end
    % multi-affine terms
    idx = nchoosek(2:m,2);
    for jj = 1:size(idx,1)
        ii = ii + 1;
        C(:,:,ii) = A(:,:,idx(jj,1))*B(:,:,idx(jj,2)) + A(:,:,idx(jj,2))*B(:,:,idx(jj,1));
    end
    % quadratic terms
    for jj = 2:m
        ii = ii + 1;
        C(:,:,ii) = A(:,:,jj)*B(:,:,jj);
    end
    
end



