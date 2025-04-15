function pdGout = lft(pdG, pdK, points)

% LFT returns the closed loop interconnection
%
% Use:
%   pdGcl = LFT(pdG, pdK, [points])
%
% where:
%   pdG is the plant model (pass object)
%   pdK is the controller (pass or LTI object)
%
% when pdG is pass object, pdK is a pass or LTI object and the meas/ctrl 
% channel is not parameter dependent, the resulting model is also 
% a pass object.
% 
% When pdK is a pcss object, pdGcl is a ppss in which the model 
% corresponding to each point in the parameter set is the lft 
% interconnection of the underlying LTI models for the plant and the 
% controller
%
% see also pass, pcss


% fbianchi - 2024-01-26


if isa(pdK, 'pass')
    
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
            error('PASS:LFT','the resulting model has no input or output')
        end
        
        % u->y transfer 
        S.type = '()';
        S.subs{1} = nz + (1:ny);
        S.subs{2} = nw + (1:nu);
        pdG_yu = subsref(pdG,S);
%         [ioDiag,ioMsg] = ispd(pdG_yu);
%         if any(ioDiag(2:3))
%             error(ioMsg)
%         end
        [ioDiagG,ioMsgG] = ispd(pdG_yu);
        [ioDiagK,ioMsgK] = ispd(pdK);
        if ~((ioDiagG(1) && ioDiagK(2)) || (ioDiagG(1) && ioDiagK(2)) || all(ioDiagG))
            error('PASS:LFT','pdG: %s\npdK: %s',ioMsgG,ioMsgK)
        end
        
        n = ns + nk;
        A = zeros(n,n,nv);
        B = zeros(n,nw,nv);
        C = zeros(nz,n,nv);
        D = zeros(nz,nw,nv);
        
        % plant matrices
        a = pdG.A;
        b1 = pdG.B(:,1:nw,:);     b2 = pdG.B(:,nw+1:end,:);
        c1 = pdG.C(1:nz,:,:);     c2 = pdG.C(nz+1:end,:,:);
        d11 = pdG.D(1:nz,1:nw,:); d12 = pdG.D(1:nz,nw+1:end,:);
        d21 = pdG.D(nz+1:end,1:nw,:); 
        d22 = pdG.D(nz+1:end,nw+1:end,:);

        if all(d22,'all')
            error('PASS:LFT','d22 must be zero')
        end
        
        % controller matrices
        ak = pdK.A;             bk = pdK.B;
        ck = pdK.C;             dk = pdK.D;
        
        % closed-loop matrices
        S21 = [c2 d21];
        S12 = [b2; d12];
        M1 = multMat(bk,S21);
        M2 = multMat(S12,ck);
        M3 = multMat(S12,multMat(dk,S21));
        
        A = [a + M3(1:ns,1:ns,:), M2(1:ns,:,:);...
             M1(:,1:ns,:),        ak];
        B = [b1 + M3(1:ns,ns+1:end,:); M1(:,ns+1:end,:)];
        C = [c1 + M3(ns+1:end,1:ns,:) M2(ns+1:end,:,:)];
        D = d11 + M3(ns+1:end,ns+1:end,:);
        
%         S.subs(2) = [];
% %         S.subs{1} = 1;
% %         G = subsref(pdG,S);
% %         [~,~,b2,~,c2,~,d12,d21] = parsysdata(G,[ny nu]);
%         for ii = 1:nv
%             
%             S.subs{1} = ii;
%             G = subsref(pdG,S);
%             [a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(G,[ny nu]);
%             if any(d22,'all')
%                 error('PASS:LFT','d22 must be zero')
%             end
%             K = subsref(pdK,S);
%             [ak,bk,ck,dk] = ssdata(K);
%             
%             if (nk == 0)
%                 % state feedback
%                 A(:,:,ii) = a + b2*dk;
%                 B(:,:,ii) = b1;
%                 C(:,:,ii) = c1 + d12*dk;
%                 D(:,:,ii) = d11;
%             else
%                 A(:,:,ii) = [a+b2*dk*c2 b2*ck; bk*c2 ak];
%                 B(:,:,ii) = [b1+b2*dk*d21; bk*d21];
%                 C(:,:,ii) = [c1+d12*dk*c2 d12*ck];
%                 D(:,:,ii) = d11+d12*dk*d21;
%             end
%             
%         end
        
        pdGout = pass(A,B,C,D,pdG.parset,...
                     'StateName',[pdG.StateName; pdK.StateName],...
                     'InputName',pdG.InputName(1:nw),...
                     'OutputName',pdG.OutputName(1:nz));
        
        
    else
       error('PASS:LFT','the parameter set of both model must be equal') 
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
        error('PASS:LFT','the resulting model has no input or output')
    end

    % controller matrices
    [ak,bk,ck,dk] = ssdata(pdK);
    
    % u->y transfer
    S.type = '()';
    S.subs{1} = nz + (1:ny);
    S.subs{2} = nw + (1:nu);
    pdG_yu = subsref(pdG,S);
    [ioDiag,ioMsg] = ispd(pdG_yu);
    if all(ioDiag(2:3)) && any(dk)
        error('PASS:LFT','the resulting model will not be affine\nbecause %s',ioMsg)
    end
    
    A = zeros(ns+nk,ns+nk,nv);
    B = zeros(ns+nk,nw,nv);
    C = zeros(nz,ns+nk,nv);
    D = zeros(nz,nw,nv);
    
    S.subs(2) = [];
%     S.subs{1} = 1;
%     G = subsref(pdG,S);
%     [a,b1,b2,c1,c2,d11,d12,d21] = parsysdata(G,[ny nu]);
% 
%     if (nk == 0)
%         % state feedback
%         A(:,:,1) = a + b2*dk;
%         B(:,:,1) = b1;
%         C(:,:,1) = c1 + d12*dk;
%         D(:,:,1) = d11;
%     else
%         A(:,:,1) = [a + b2*dk*c2 b2*ck; bk*c2 ak];
%         B(:,:,1) = [b1 + b2*dk*d21; bk*d21];
%         C(:,:,1) = [c1 + d12*dk*c2 d12*ck];
%         D(:,:,1) = d11 + d12*dk*d21;
%     end
    
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
            A(:,:,ii) = [a + b2*dk*c2 b2*ck; bk*c2 zeros(nk,nk)];
            B(:,:,ii) = [b1 + b2*dk*d21; bk*d21];
            C(:,:,ii) = [c1 + d12*dk*c2 d12*ck];
            D(:,:,ii) = d11 + d12*dk*d21;
        end
        
    end
    
    pdGout = pass(A,B,C,D,pdG.parset,...
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
        error('PASS:LFT','the resulting model has no input or output')
    end
    
    if (nargin == 3)
        if ~isnumeric(points) || ~(size(points,1) == np)
            error('PASS:LFT',...
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

    Pset = pset.Gral(points,pdG.parset.rate);
    pdGout = ppss(A,B,C,D,Pset,...
        'StateName',[pdG.StateName; K.StateName],...
        'InputName',pdG.InputName(1:nw),...
        'OutputName',pdG.OutputName(1:nz));
    
else
    error('PASS:LFT','PDK must be a pass, LTI or pcss object')
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = multMat(A,B)

nr = size(A,1);
nc = size(B,2);
nt = size(A,3);

if all(A(:,:,2:end) == 0,'all')
    % A is constant (case affine)
    na = 1;
else
    na = size(A,3);
end

if all(B(:,:,2:end) == 0,'all')
    % B is constant (case affine)
    nb = 1;
else
    nb = size(B,3);
end

C = zeros(nr,nc,nt);
if (na == 1)
    % A is constant
    for ii = 1:nb
       C(:,:,ii) = A(:,:,1)*B(:,:,ii);
    end
else
    % B is constant
    for ii = 1:na
       C(:,:,ii) = A(:,:,ii)*B(:,:,1);
    end
end    
    

