function Go = horzcat(varargin)

% HORZCAT returns a ppss model corresponding to the horizontal concatenation
%
% Use:
%   Gout = HORZCAT(G1,G2,G3,...,Gn) or
%   Gout = [G1,G2,G3,...,Gn]
%
% where:
%   G1, G2, ..., Gn are ppss objects or ss, tf, or zpk objects
%
% in case of more than one ppss object, all of them must have the same
% parameter set
%
% See also ppss, append, series

% fbianchi - 2025-01-16


if (nargin == 2)
    
    G1 = varargin{1};
    G2 = varargin{2};

    if isnumeric(G1)
        G1 = ss(G1);
    end
    if isnumeric(G2)
        G2 = ss(G2);
    end
    
    % model dimensions
    [no1,ni1] = iosize(G1); ns1 = order(G1);
    [no2,ni2] = iosize(G2); ns2 = order(G2);
    
    if (no1 ~= no2)
        error('PGSS:horzcat:inputError','G1 or G2 must have the same number of outputs')
    end
    
    if isa(G1,'ppss') && isa(G2,'ppss')
        
        bool = isequal(G1.parset,G2.parset);
        nv1 = nsys(G1);
        nv2 = nsys(G2);
        
        if bool && (nv1 == nv2)
            Ao(ns1+ns2,ns1+ns2,nv1) = 0;
            Bo(ns1+ns2,ni1+ni2,nv1) = 0;
            Co(no1,ns1+ns2,nv1) = 0;
            Do(no1,ni1+ni2,nv1) = 0;
            for ii = 1:nsys(G1)
                Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), G2.A(:,:,ii));
                Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), G2.B(:,:,ii));
                Co(:,:,ii) = [G1.C(:,:,ii), G2.C(:,:,ii)];
                Do(:,:,ii) = [G1.D(:,:,ii), G2.D(:,:,ii)];
            end
            
            pv  = G1.parset;
        else
            error('PPSS:horzcat:inputError','G1 or G2 must have the same parameter set and the same number of vertices')
        end
        
        
    elseif (isa(G1,'ss') || isa(G1,'zpk') || isa(G1,'tf') || isnumeric(G1))...
            && isa(G2,'ppss')
        
        G1 = ss(G1);
        [A1,B1,C1,D1] = ssdata(G1);
        
        % model dimensions
        ns1 = order(G1);    [no1,ni1] = iosize(G1);
        ns2 = order(G2);    [no2,ni2] = iosize(G2);
        nv  = nsys(G2);
        
        Ao = zeros(ns1+ns2,ns1+ns2,nv);
        Bo = zeros(ns1+ns2,ni1+ni2,nv);
        Co = zeros(no1,ns1+ns2,nv);
        Do = zeros(no1,ni1+ni2,nv);
        
        for ii = 1:nv
            Ao(:,:,ii) = blkdiag(A1, G2.A(:,:,ii));
            Bo(:,:,ii) = blkdiag(B1, G2.B(:,:,ii));
            Co(:,:,ii) = [C1, G2.C(:,:,ii)];
            Do(:,:,ii) = [D1, G2.D(:,:,ii)];
        end
        
        pv = G2.parset;
        
    elseif (isa(G2,'ss') || isa(G2,'zpk') || isa(G2,'tf') || isnumeric(G2))...
            && isa(G1,'ppss')
        
        G2 = ss(G2);
        [A2,B2,C2,D2] = ssdata(G2);
        
        % model dimensions
        ns1 = order(G1);    [no1,ni1] = iosize(G1);
        ns2 = order(G2);    [no2,ni2] = iosize(G2);
        nv  = nsys(G1);
        
        Ao = zeros(ns1+ns2,ns1+ns2,nv);
        Bo = zeros(ns1+ns2,ni1+ni2,nv);
        Co = zeros(no1,ns1+ns2,nv);
        Do = zeros(no1,ni1+ni2,nv);
        
        for ii = 1:nv
            Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), A2);
            Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), B2);
            Co(:,:,ii) = [G1.C(:,:,ii), C2];
            Do(:,:,ii) = [G1.D(:,:,ii), D2];
        end
        
        pv = G1.parset;
        
    else
        error('PPSS:horzcat','G1 or G2 must be LTI systems or ppss objects')
    end
    
    Go = ppss(Ao, Bo, Co, Do, pv,...
            'StateName',[G1.StateName; G2.StateName],...
            'InputName',[G1.InputName; G2.InputName],...
            'OutputName',[G1.OutputName]);
else
    Go = varargin{1};
    for ii = 2:nargin
        G2 = varargin{ii};
        Go = horzcat(Go,G2);
    end
end
    