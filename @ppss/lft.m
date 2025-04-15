function pdGout = lft(pdG, pdK, points)

% LFT returns the closed loop interconnection
%
% Use:
%   pdGcl = LFT(pdG, pdK, [points])
%
% where:
%   pdG is the plant model (ppss object)
%   pdK is the controller (ppss, pcss or LTI object)
%
% when pdK is a ppss and the meas/ctrl channel is not parameter dependent, 
% the resulting model is also a ppss object. The meas or ctrl channel can be 
% parameter dependent if the controller matrices B and D or C and D are 
% constant.
% 
% Otherwise, LFT will produce a 3D ss object in which the LTI models in the
% description corresponds to the LFT interconnection of the LTI models for 
% the plant and the controller evaluated at the frozen points stored in the
% parameter set. That is,
%   pdGcl = LFT(pdG, pdK)
% will be the same as
%   G = ss(pdG); K = ss(pdK);
%   pdGcl = ss(LFT(G,K));
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
        [ioDiagG,ioMsgG] = ispd(pdG_yu);
        [ioDiagK,ioMsgK] = ispd(pdK);
        if ((ioDiagG(2) && ioDiagK(3)) || (ioDiagG(3) && ioDiagK(2)) || all(ioDiagG(2:3)))
            warning('PPSS:LFT','pdG: %s\npdK: %s\n%s',ioMsgG,ioMsgK,...
                    'The resulting model will be ss 3D-object')
            G = ss(pdG);
            K = ss(pdK);
            pdGout = lft(G,K);
            return
        end
        
        A = zeros(ns+nk,ns+nk,nv);
        B = zeros(ns+nk,nw,nv);
        C = zeros(nz,ns+nk,nv);
        D = zeros(nz,nw,nv);
        
        S.subs(2) = [];
        for ii = 1:nv
            
            S.subs{1} = ii;
            G = subsref(pdG,S);
            [a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(G,[ny nu]);
            if any(d22,'all')
                error('PPSS:LFT','d22 must be zero')
            end
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
    
    A = zeros(ns+nk,ns+nk,nv);
    B = zeros(ns+nk,nw,nv);
    C = zeros(nz,ns+nk,nv);
    D = zeros(nz,nw,nv);
    
    [ak,bk,ck,dk] = ssdata(pdK);
    if all(ioDiag(2:3)) && any(dk > 0)
        error(ioMsg)
        error('PPSS:LFT','pdG: %s\nand dk is not zero',ioMsg)
    end
    
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
    
    warning('pdK is pcss, then the resulting model is a 3D LTI array');
    
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

%     Pset = pset.Gral(points,pset,pdG.parset.rate);
%     pdGout = ppss(A,B,C,D,Pset,...
    pdGout = ss(A,B,C,D,...
                'StateName',[pdG.StateName; K.StateName],...
                'InputName',pdG.InputName(1:nw),...
                'OutputName',pdG.OutputName(1:nz));

else
    error('PPSS:LFT','Both model must be ppss object')
        
end
