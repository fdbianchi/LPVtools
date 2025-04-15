function pdGout = lft(pdG, pdK, points)

% LFT returns the closed loop interconnection
%
% Use:
%   pdGcl = LFT(pdG, pdK, [points])
%
% where:
%   pdG is the plant model (pgss object)
%   pdK is the controller (pgss or LTI object)
%
% when pdG is pgss object, pdK is a pgss or LTI object and the meas/ctrl 
% channel is not parameter dependent, the resulting model is also 
% a pgss object.
% 
% When pdK is a pcss object, pdGcl is a ppss in which the model 
% corresponding to each point in the parameter set is the lft 
% interconnection of the underlying LTI models for the plant and the 
% controller
%
% see also pgss, pcss

% fbianchi - 2024-01-26

if isa(pdG, 'pgss') && isa(pdK, 'pgss') && (nargin < 3)
    
    if isequal(pdG.parset, pdK.parset) && isequal(pdG.parfcn, pdK.parfcn)
        
        [no,ni,ns,~,nv] = size(pdG);
        [nu,ny,nk] = size(pdK);
        if (nk == 0)
            % state feedback
            ny = 0;
        end
        
        nw = ni - nu;
        nz = no - ny;
        
        if (nw == 0) || (nz == 0)
            error('PGSS:LFT','the resulting model has no input or output')
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
        
        S.subs(2) = [];
        S.subs{1} = 1;
        K = subsref(pdG,S);
        [~,~,b2,~,c2,~,d12,d21] = parsysdata(K,[ny nu]);
        for ii = 1:nv
            
            S.subs{1} = ii;
            K = subsref(pdG,S);
            [a,b1,~,c1,~,d11] = parsysdata(K,[ny nu]);
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
        
        pdGout = pgss(A,B,C,D,pdG.parset,pdG.parfcn,...
                     'StateName',[pdG.StateName; pdK.StateName],...
                     'InputName',pdG.InputName(1:nw),...
                     'OutputName',pdG.OutputName(1:nz));
        
        
    else
       error('PGSS:LFT',...
           'The parameter set and functions of both model must be equal')
       
    end
    
elseif (isa(pdG,'ss') || isa(pdG,'zpk') || isa(pdG,'tf')) && (nargin < 3)
    
    [no,ni,ns,~,nv] = size(pdK);
    [nu,ny] = size(pdG);
    nk = order(pdG);
    if (nk == 0)
        % state feedback
        ny = 0;
    end
    
    nw = ni - nu;
    nz = no - ny;
    
    if (nw == 0) || (nz == 0)
        error('PGSS:LFT','the resulting model has no input or output')
    end
    
    % u->y transfer
    S.type = '()';
    S.subs{1} = nz + (1:ny);
    S.subs{2} = nw + (1:nu);
    pdK_yu = subsref(pdK,S);
    [ioDiag,ioMsg] = ispd(pdK_yu);
    if all(ioDiag(2:3))
        error(ioMsg)
    end
    
    A = zeros(ns+nk,ns+nk,nv);
    B = zeros(ns+nk,nw,nv);
    C = zeros(nz,ns+nk,nv);
    D = zeros(nz,nw,nv);
    
    S.subs(2) = [];
    S.subs{1} = 1;
    K = subsref(pdK,S);
    [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(K,[ny nu]);

    [ak,bk,ck,dk] = ssdata(pdG);
    
    if (nk == 0)
        % state feedback
        A(:,:,1) = a + b2*dk;
        B(:,:,1) = b1;
        C(:,:,1) = c1 + d12*dk;
        D(:,:,1) = d11;
    else
        A(:,:,1) = [a + b2*dk*c2 b2*ck; bk*c2 ak];
        B(:,:,1) = [b1 + b2*dk*d21; bk*d21];
        C(:,:,1) = [c1 + d12*dk*c2 d12*ck];
        D(:,:,1) = d11 + d12*dk*d21;
    end
    
    for ii = 2:nv
        
        S.subs{1} = ii;
        K = subsref(pdK,S);
        
        [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(K,[ny nu]);
        
        if (nk == 0)
            % state feedback
            A(:,:,ii) = a + b2*dk;
            B(:,:,ii) = b1;
            C(:,:,ii) = c1 + d12*dk;
            D(:,:,ii) = d11;
        else
            A(:,:,ii) = [a + b2*dk*c2 b2*ck; bk*c2 zeros(nk,nk)];
            B(:,:,ii) = [b1 + b2*dk*d21; bk*d21];
            C(:,:,ii) = [c1 + d12*dk*c2 d12*ck];
            D(:,:,ii) = d11 + d12*dk*d21;
        end
        
    end
    
    pdGout = pgss(A,B,C,D,pdK.parset,pdK.parfcn,...
                'StateName',[pdK.StateName; pdD.StateName],...
                'InputName',pdK.InputName(1:nw),...
                'OutputName',pdK.OutputName(1:nz));

elseif (isa(pdK,'ss') || isa(pdK,'zpk') || isa(pdK,'tf')) && (nargin < 3)
    
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
        error('PGSS:LFT','the resulting model has no input or output')
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
    
    S.subs(2) = [];
    S.subs{1} = 1;
    K = subsref(pdG,S);
    [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(K,[ny nu]);

    [ak,bk,ck,dk] = ssdata(pdK);
    
    if (nk == 0)
        % state feedback
        A(:,:,1) = a + b2*dk;
        B(:,:,1) = b1;
        C(:,:,1) = c1 + d12*dk;
        D(:,:,1) = d11;
    else
        A(:,:,1) = [a + b2*dk*c2 b2*ck; bk*c2 ak];
        B(:,:,1) = [b1 + b2*dk*d21; bk*d21];
        C(:,:,1) = [c1 + d12*dk*c2 d12*ck];
        D(:,:,1) = d11 + d12*dk*d21;
    end
    
    for ii = 2:nv
        
        S.subs{1} = ii;
        K = subsref(pdG,S);
        
        [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(K,[ny nu]);
        
        if (nk == 0)
            % state feedback
            A(:,:,ii) = a + b2*dk;
            B(:,:,ii) = b1;
            C(:,:,ii) = c1 + d12*dk;
            D(:,:,ii) = d11;
        else
            A(:,:,ii) = [a + b2*dk*c2 b2*ck; bk*c2 zeros(nk,nk)];
            B(:,:,ii) = [b1 + b2*dk*d21; bk*d21];
            C(:,:,ii) = [c1 + d12*dk*c2 d12*ck];
            D(:,:,ii) = d11 + d12*dk*d21;
        end
        
    end
    
    pdGout = pgss(A,B,C,D,pdG.parset,pdG.parfcn,...
                'StateName',[pdG.StateName; pdK.StateName],...
                'InputName',pdG.InputName(1:nw),...
                'OutputName',pdG.OutputName(1:nz));

elseif isa(pdK,'pcss') || (nargin == 3)
    
    [no,ni,ns,np] = size(pdG);
    [nu,ny,nk] = size(pdK);
    if (nk == 0)
        % state feedback
        ny = 0;
    end
    
    nw = ni - nu;
    nz = no - ny;
    
    if (nw == 0) || (nz == 0)
        error('PGSS:LFT','the resulting model has no input or output')
    end
    
    if (nargin == 3)
        if ~isnumeric(points) || ~(size(points,1) == np)
            error('PGSS:LFT',...
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
    
    K = ss(pdG, points);
    K = ss(pdK, points);
    
    for ii = 1:nv
        
        [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(K(:,:,ii),[ny nu]);
        
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

    Pset = pset.Gral(points,pdG.parset.rate);
    pdGout = ppss(A,B,C,D,Pset,...
        'StateName',[pdG.StateName; K.StateName],...
        'InputName',pdG.InputName(1:nw),...
        'OutputName',pdG.OutputName(1:nz));

else
    error('PGSS:LFT','Both model must be pgss object')
       
end
