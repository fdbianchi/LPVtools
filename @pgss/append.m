function Go = append(varargin)

% APPEND returns a pgss model with the appended inputs and outputs
%
% Use:
%   Gout = APPEND(G1,G2,G3,...,Gn)
%
% where:
%   G1, G2, ..., Gn are pgss objects or ss, tf, or zpk objects
%
% in case of more than one pgss object, all of them must have the same 
% parset and functions
%
% See also pgss, append, series

% fbianchi - 2024-12-05

if (nargin == 2)
    
    % case with only two models
    G1 = varargin{1};
    G2 = varargin{2};
    
    % converting numeric matrices into ss objects
    if isnumeric(G1)
        G1 = ss(G1);
    end
    if isnumeric(G2)
        G2 = ss(G2);
    end
    
    % model dimensions
    [no1,ni1] = iosize(G1); ns1 = order(G1);
    [no2,ni2] = iosize(G2); ns2 = order(G2);
    
    if isa(G1,'pgss') && isa(G2,'pgss')
        % both models LPV
        
        bool1 = isequal(G1.parset,G2.parset);
        bool2 = strcmp(func2str(G1.parfcn),func2str(G2.parfcn));
        
        if (bool1 && bool2)
            
            nv = nsys(G1);
            Ao(ns1+ns2,ns1+ns2,nv) = 0;
            Bo(ns1+ns2,ni1+ni2,nv) = 0;
            Co(no1+no2,ns1+ns2,nv) = 0;
            Do(no1+no2,ni1+ni2,nv) = 0;
            for ii = 1:nv
                Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), G2.A(:,:,ii));
                Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), G2.B(:,:,ii));
                Co(:,:,ii) = blkdiag(G1.C(:,:,ii), G2.C(:,:,ii));
                Do(:,:,ii) = blkdiag(G1.D(:,:,ii), G2.D(:,:,ii));
            end
            
            pv = G1.parset;
            fcn = G1.parfcn;
            
        else
            error('PGSS:append:inputError','G1 or G2 must have the same parameter set and functions')
        end
        
    elseif (isa(G1,'ss') || isa(G1,'zpk') || isa(G1,'tf')) && isa(G2,'pgss')
        
        [A1,B1,C1,D1] = ssdata(ss(G1));
        
        Ao(:,:,1) = blkdiag(A1, G2.A(:,:,1));
        Bo(:,:,1) = blkdiag(B1, G2.B(:,:,1));
        Co(:,:,1) = blkdiag(C1, G2.C(:,:,1));
        Do(:,:,1) = blkdiag(D1, G2.D(:,:,1));
        for ii = 2:nsys(G2)
            Ao(:,:,ii) = blkdiag(zeros(ns1), G2.A(:,:,ii));
            Bo(:,:,ii) = blkdiag(zeros(ns1,ni1), G2.B(:,:,ii));
            Co(:,:,ii) = blkdiag(zeros(no1,ns1), G2.C(:,:,ii));
            Do(:,:,ii) = blkdiag(zeros(no1,ni1), G2.D(:,:,ii));
        end
        
        pv = G2.parset;
        fcn = G2.parfcn;
        
    elseif (isa(G2,'ss') || isa(G2,'zpk') || isa(G2,'tf')) && isa(G1,'pgss')
        
        [A2,B2,C2,D2] = ssdata(ss(G2));
        
        Ao(:,:,1) = blkdiag(G1.A(:,:,1), A2);
        Bo(:,:,1) = blkdiag(G1.B(:,:,1), B2);
        Co(:,:,1) = blkdiag(G1.C(:,:,1), C2);
        Do(:,:,1) = blkdiag(G1.D(:,:,1), D2);
        for ii = 2:nsys(G1)
            Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), zeros(ns2));
            Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), zeros(ns2,ni2));
            Co(:,:,ii) = blkdiag(G1.C(:,:,ii), zeros(no2,ns2));
            Do(:,:,ii) = blkdiag(G1.D(:,:,ii), zeros(no2,ni2));
        end
        
        pv = G1.parset;
        fcn = G1.parfcn;
        
    else
        error('PGSS:append:inputError','G1 and G2 must be an LTI system or pgss objects')
        
    end
    
    Go = pgss(Ao, Bo, Co, Do, pv, fcn,...
            'StateName',[G1.StateName; G2.StateName],...
            'InputName',[G1.InputName; G2.InputName],...
            'OutputName',[G1.OutputName; G2.OutputName]);
    
else
    
    % case with more than two models
    Go = varargin{1};
    
    for ii = 2:nargin
        G2 = varargin{ii};
        Go = append(Go,G2);
    end
    
end