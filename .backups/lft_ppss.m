function pdGout = lft(pdG, pdK, points)

% LFT returns the closed loop interconnection
%
% Use:
%   pdGcl = LFT(pdG, pdK, [points])
%
% where:
%   pdG is the plant model (ppss object)
%   pdK is the controller (ppss or LTI object)
%
% when pdG is ppss object, pdK is a ppss or LTI object and the meas/ctrl 
% channel is not parameter dependent, the resulting model is also 
% a ppss object.
% 
% When pdK is a pcss object, pdGcl is a ppss in which the model 
% corresponding to each point in the parameter set is the lft 
% interconnection of the underlying LTI models for the plant and the 
% controller
%
% see also ppss, pcss

% fbianchi - 2024-01-26

if isa(pdK, 'ppss')
    
    if isequal(pdG.parset, pdK.parset)
        
        [no,ni,ns,~,nv] = size(pdG);
        [nu,ny,nk] = size(pdK);
        if (nk == 0)
            % state feedback
            ny = 0;
        end
        
        nw = ni - nu;
        nz = no - ny;
        
        if (nw == 0) || (nz == 0)
            error('the resulting model has no input or output')
        end
        
        % u->y transfer 
        S.type = '()';
        S.subs{1} = nz + (1:ny);
        S.subs{2} = nw + (1:nu);
        pdG_yu = subsref(pdG,S);
        [ioDiag,ioMsg] = ispd(pdG_yu);
        if any(ioDiag(2:3))
            error(ioMsg)
        end
        
        A = zeros(ns+nk,ns+nk,nv);
        B = zeros(ns+nk,nw,nv);
        C = zeros(nz,ns+nk,nv);
        D = zeros(nz,nw,nv);
        
        S.subs(2) = [];
        S.subs{1} = 1;
        G = subsref(pdG,S);
        [~,~,b2,~,c2,~,d12,d21] = parsysdata(G,[ny nu]);
        for ii = 1:nv
            
            S.subs{1} = ii;
            G = subsref(pdG,S);
            [a,b1,~,c1,~,d11] = parsysdata(G,[ny nu]);
            K = subsref(pdK,S);
            [ak,bk,ck,dk] = ssdata(K);
            
            if (nk == 0)
                % state feedback
                A(:,:,ii) = a + b2*dk;
                B(:,:,ii) = b1;
                C(:,:,ii) = c1 + d12*dk;
                D(:,:,ii) = d11;
            else
                A(:,:,ii) = [a+b2*dk*c2 b2*ck; bk*c2 ak];
                B(:,:,ii) = [b1+b2*dk*d21; bk*d21];
                C(:,:,ii) = [c1+d12*dk*c2 d12*ck];
                D(:,:,ii) = d11+d12*dk*d21;
            end
            
        end
        
        pdGout = ppss(A,B,C,D,pdG.parset,...
                     'StateName',[pdG.StateName; pdK.StateName],...
                     'InputName',pdG.InputName(1:nw),...
                     'OutputName',pdG.OutputName(1:nz));
        
        
    else
       error('PPSS:LFT','the parameter set of both model must be equal') 
    end
    
elseif (isa(pdK,'ss') || isa(pdK,'zpk') || isa(pdK,'tf'))% || isnumeric(pdK))    
    
    [no,ni,ns,~,nv] = size(pdG);
    [nu,ny] = size(pdK);
    nk = order(pdK);
    if (nk == 0)
        % state feedback
        ny = 0;
    end
    
    nw = ni - nu;
    nz = no - ny;
    
    if (nw == 0) || (nz == 0)
        error('PPSS:LFT','the resulting model has no input or output')
    end
    
    % u->y transfer
    S.type = '()';
    S.subs{1} = nz + (1:ny);
    S.subs{2} = nw + (1:nu);
    pdG_yu = subsref(pdG,S);
    [ioDiag,ioMsg] = ispd(pdG_yu);
    if all(ioDiag(2:3))
        error(ioMsg)
    end
    
    A = zeros(ns+nk,ns+nk,nv);
    B = zeros(ns+nk,nw,nv);
    C = zeros(nz,ns+nk,nv);
    D = zeros(nz,nw,nv);
    
    [ak,bk,ck,dk] = ssdata(pdK);
    
    S.subs(2) = [];
    for ii = 1:nv
        
        S.subs{1} = ii;
        G = subsref(pdG,S);
        [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(G,[ny nu]);
        
        if (nk == 0)
            % state feedback
            A(:,:,ii) = a + b2*dk;
            B(:,:,ii) = b1;
            C(:,:,ii) = c1 + d12*dk;
            D(:,:,ii) = d11;
        else
            A(:,:,ii) = [a+b2*dk*c2 b2*ck; bk*c2 ak];
            B(:,:,ii) = [b1+b2*dk*d21; bk*d21];
            C(:,:,ii) = [c1+d12*dk*c2 d12*ck];
            D(:,:,ii) = d11+d12*dk*d21;
        end
    end

    
    pdGout = ppss(A,B,C,D,pdG.parset,...
        'StateName',[pdG.StateName; pdK.StateName],...
        'InputName',pdG.InputName(1:nw),...
        'OutputName',pdG.OutputName(1:nz));

elseif isa(pdK,'pcss')
    
    [no,ni,ns,np,~] = size(pdG);
    [nu,ny,nk] = size(pdK);
    if (nk == 0)
        % state feedback
        ny = 0;
    end
    
    nw = ni - nu;
    nz = no - ny;
    
    if (nw == 0) || (nz == 0)
        error('PPSS:LFT','the resulting model has no input or output')
    end
    
    if (nargin == 3)
        if ~isnumeric(points) || ~(size(points,1) == np)
            error('PPSS:LFT',...
                  'The 3rd arg must be numeric matrix with %d rows',np)
        end
    else
        points = pdG.parset.points;
    end
    nv = size(points,2);
    
    A = zeros(ns+nk,ns+nk,nv);
    B = zeros(ns+nk,nw,nv);
    C = zeros(nz,ns+nk,nv);
    D = zeros(nz,nw,nv);
    
    G = ss(pdG, points);
    K = ss(pdK, points);
    
    for ii = 1:nv
        
        [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(G(:,:,ii),[ny nu]);
        
        [ak,bk,ck,dk] = ssdata(K(:,:,ii));
        
        if (nk == 0)
            % state feedback
            A(:,:,ii) = a + b2*dk;
            B(:,:,ii) = b1;
            C(:,:,ii) = c1 + d12*dk;
            D(:,:,ii) = d11;
        else
            A(:,:,ii) = [a+b2*dk*c2 b2*ck; bk*c2 ak];
            B(:,:,ii) = [b1+b2*dk*d21; bk*d21];
            C(:,:,ii) = [c1+d12*dk*c2 d12*ck];
            D(:,:,ii) = d11+d12*dk*d21;
        end
    end

    Pset = pset.Gral(points,pset,pdG.parset.rate);
    pdGout = ppss(A,B,C,D,Pset,...
        'StateName',[pdG.StateName; K.StateName],...
        'InputName',pdG.InputName(1:nw),...
        'OutputName',pdG.OutputName(1:nz));

else
    error('PPSS:LFT','Both model must be ppss object')
        
end
