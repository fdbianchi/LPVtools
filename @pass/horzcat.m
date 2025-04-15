function Go = horzcat(varargin)

% HORZCAT returns a pass model corresponding to the horizontal concatenation
%
% Use:
%   Gout = HORZCAT(G1,G2,G3,...,Gn) or
%   Gout = [G1,G2,G3,...,Gn]
%
% where:
%   G1, G2, ..., Gn are pass objects or ss, tf, or zpk objects
%
% in case of more than one pass object, all of them must have the same parset
%
% See also pass, append, series

% fbianchi - 2025-01-16

if (nargin == 2)
    % with only two models
    
    G1 = varargin{1};
    G2 = varargin{2};

    % converting numeric arrays into ss objects
    if isnumeric(G1)
        G1 = ss(G1);
    end
    if isnumeric(G2)
        G2 = ss(G2);
    end
    
    % model dimensions
    [no1,ni1] = size(G1); ns1 = order(G1); 
    [no2,ni2] = size(G2); ns2 = order(G2); 
    if (no1 ~= no2)
        error('PASS:horzcat:inputError','G1 or G2 must have the same number of outputs')
    end
    
    if isa(G1,'pass') && isa(G2,'pass')
        
        % case: two lpv models
        
        bool = isequal(G1.parset,G2.parset);
        if bool
            nv = nsys(G1);
            Ao(ns1+ns2,ns1+ns2,nv) = 0;
            Bo(ns1+ns2,ni1+ni2,nv) = 0;
            Co(no1,ns1+ns2,nv) = 0;
            Do(no1,ni1+ni2,nv) = 0;
            for ii = 1:nsys(G1)
                Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), G2.A(:,:,ii));
                Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), G2.B(:,:,ii));
                Co(:,:,ii) = [G1.C(:,:,ii), G2.C(:,:,ii)];
                Do(:,:,ii) = [G1.D(:,:,ii), G2.D(:,:,ii)];
            end
            
            pv  = G1.parset;
        else
            error('PASS:horzcat:inputError','G1 or G2 must have the same parameter set')
        end
        
    elseif (isa(G1,'ss') || isa(G1,'zpk') || isa(G1,'tf')) && isa(G2,'pass')
        
        % case: G2 lti
        
        [A1,B1,C1,D1] = ssdata(ss(G1));
        
        Ao(:,:,1) = blkdiag(A1, G2.A(:,:,1));
        Bo(:,:,1) = blkdiag(B1, G2.B(:,:,1));
        Co(:,:,1) = [C1, G2.C(:,:,1)];
        Do(:,:,1) = [D1, G2.D(:,:,1)];
        for ii = 2:nsys(G2)
            Ao(:,:,ii) = blkdiag(zeros(ns1), G2.A(:,:,ii));
            Bo(:,:,ii) = blkdiag(zeros(ns1,ni1), G2.B(:,:,ii));
            Co(:,:,ii) = [zeros(no1,ns1), G2.C(:,:,ii)];
            Do(:,:,ii) = [zeros(no1,ni1), G2.D(:,:,ii)];
        end
        
        pv  = G2.parset;
        
    elseif (isa(G2,'ss') || isa(G2,'zpk') || isa(G2,'tf')) && isa(G1,'pass')
        
        % case: G1 lti

        [A2,B2,C2,D2] = ssdata(ss(G2));
        
        Ao(:,:,1) = blkdiag(G1.A(:,:,1), A2);
        Bo(:,:,1) = blkdiag(G1.B(:,:,1), B2);
        Co(:,:,1) = [G1.C(:,:,1), C2];
        Do(:,:,1) = [G1.D(:,:,1), D2];
        for ii = 2:nsys(G1)
            Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), zeros(ns2));
            Bo(:,:,ii) = blkdiag(G1.B(:,:,ii), zeros(ns2,ni2));
            Co(:,:,ii) = [G1.C(:,:,ii), zeros(no2,ns2)];
            Do(:,:,ii) = [G1.D(:,:,ii), zeros(no2,ni2)];
        end
        
        pv  = G1.parset;
        
    else
        error('PASS:horzcat:inputError','G1 and G2 must be LTI models or pass objects')
        
    end
    
    % resulting models
    Go = pass(Ao, Bo, Co, Do, pv,...
            'StateName',[G1.StateName; G2.StateName],...
            'InputName',[G1.InputName; G2.InputName],...
            'OutputName',[G1.OutputName]);
    
else
    
    % for more than two models
    Go = varargin{1};
    for ii = 2:nargin
        G2 = varargin{ii};
        Go = horzcat(Go,G2);
    end
    
end