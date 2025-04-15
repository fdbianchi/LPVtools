function Go = series(G1,G2,out1,inp2)

% SERIES connects LPV or LTI models in series
%
% Use:
%   G3 = SERIES(G1,G2,out1,inp2)
%
% where:
%   G1 or G2 are pgss. ss, tf, or zpk objects
%   out1,inp2 indices of indicating the inputs and outputs connected in 
%   series, see lti.series
%
% See also pgss, series

% fbianchi - 2021-04-06
% fbianchi - 2024-12-09 - rev. 


if (nargin == 2)
    
    % with only models as input arguments
    
    % converting numeric arrays into lti objects
    if isnumeric(G1)
        G1 = ss(G1);
    end
    if isnumeric(G2)
        G2 = ss(G2);
    end
    
    % model dimensions
    [no1,ni1] = size(G1); ns1 = order(G1); 
    [no2,ni2] = size(G2); ns2 = order(G2); 
    
    if (ni2 ~= no1)
        error('PGSS:series:inputError','The number of inputs and outpus is inconsistent')
        
    else
        
        if isa(G1,'pgss') && isa(G2,'pgss')
            
            % case two lpv model
            
            bool1 = isequal(G1.parset,G2.parset);
            bool2 = strcmp(func2str(G1.parfcn),func2str(G2.parfcn));
            % checking if parameter dependence channels to ensure a
            % resulting pgss object
            dgn1 = ispd(G1);    dgn2 = ispd(G2);
            
            if (bool1 && bool2) && ~(dgn1(3) && dgn2(2))
                
                nv = nsys(G1);
                if dgn1(3)
                    % output channel of G1 parameter dependent
                    G2.B = repmat(G2.B(:,:,1),1,1,nv);
                    G2.D = repmat(G2.D(:,:,1),1,1,nv);
                else
                    % input channel of G2 parameter dependent
                    G1.C = repmat(G1.C(:,:,1),1,1,nv);
                    G1.D = repmat(G1.D(:,:,1),1,1,nv);
                end                    
                
                Ao(ns1+ns2,ns1+ns2,nv) = 0;
                Bo(ns1+ns2,ni1,nv) = 0;
                Co(no2,ns1+ns2,nv) = 0;
                Do(no2,ni1,nv) = 0;
                for ii = 1:nv
                    Ao(:,:,ii) = [G1.A(:,:,ii), zeros(ns1,ns2)
                                  G2.B(:,:,ii)*G1.C(:,:,ii), G2.A(:,:,ii)];
                    Bo(:,:,ii) = [G1.B(:,:,ii); G2.B(:,:,ii)*G1.D(:,:,ii)];
                    Co(:,:,ii) = [G2.D(:,:,ii)*G1.C(:,:,ii), G2.C(:,:,ii)];
                    Do(:,:,ii) = G2.D(:,:,ii)*G1.D(:,:,ii);
                end
                
                pv = G1.parset;
                fcn = G1.parfcn;

            else
                error('PGSS:series:inputError',...
                    'G1 or G2 must have the same parameter set, the same functions and\n the input channel of G2 or the output channel of G1 must be parameter independent')
            end
        
        
        elseif (isa(G1,'ss') || isa(G1,'zpk') || isa(G1,'tf')) && isa(G2,'pgss')
            
            % G1 is not lpv
            G1 = ss(G1);
            [A1,B1,C1,D1] = ssdata(G1);
            
            Ao(:,:,1) = [A1, zeros(ns1,ns2)
                         G2.B(:,:,1)*C1, G2.A(:,:,1)];
            Bo(:,:,1) = [B1; G2.B(:,:,1)*D1];
            Co(:,:,1) = [G2.D(:,:,1)*C1, G2.C(:,:,1)];
            Do(:,:,1) = G2.D(:,:,1)*D1;
            for ii = 2:nsys(G2)
                Ao(:,:,ii) = [zeros(ns1,ns1+ns2);
                              G2.B(:,:,ii)*C1, G2.A(:,:,ii)];
                Bo(:,:,ii) = [zeros(ns1,ni1); G2.B(:,:,ii)*D1];
                Co(:,:,ii) = [G2.D(:,:,ii)*C1, G2.C(:,:,ii)];
                Do(:,:,ii) = G2.D(:,:,ii)*D1;
            end
            
            pv = G2.parset;
            fcn = G2.parfcn;
            
        elseif (isa(G2,'ss') || isa(G2,'zpk') || isa(G2,'tf')) && isa(G1,'pgss')
            
            % G2 is not lpv
            G2 = ss(G2);
            [A2,B2,C2,D2] = ssdata(G2);
            
            Ao(:,:,1) = [G1.A(:,:,1), zeros(ns1,ns2);
                         B2*G1.C(:,:,1), A2];
            Bo(:,:,1) = [G1.B(:,:,1); B2*G1.D(:,:,1)];
            Co(:,:,1) = [D2*G1.C(:,:,1), C2];
            Do(:,:,1) = D2*G1.D(:,:,1);
            for ii = 2:nsys(G1)
                Ao(:,:,ii) = [G1.A(:,:,ii), zeros(ns1,ns2);
                              B2*G1.C(:,:,ii), zeros(ns2)];
                Bo(:,:,ii) = [G1.B(:,:,ii); B2*G1.D(:,:,ii)];
                Co(:,:,ii) = [D2*G1.C(:,:,ii), zeros(no2,ns2)];
                Do(:,:,ii) = D2*G1.D(:,:,ii);
            end
            
            pv = G1.parset;
            fcn = G1.parfcn;
            
        else
            error('PGSS:Series:inputError','G1 or G2 must be LTI models or pgss objects')
        end
        
        Go = pgss(Ao, Bo, Co, Do, pv,fcn,...
                'StateName',[G1.StateName; G2.StateName],...
                'InputName',G1.InputName,...
                'OutputName',G2.OutputName);
    end
    
elseif (nargin == 4)
    
    % case with input and output specifications
    
    % Validate indices
    [no1,~] = iosize(G1);
    [~,ni2] = iosize(G2);
    lo = length(out1);
    li = length(inp2);
    if (li ~= lo)
        error('PGSS:Series:inputError','The length of OUT2 and INP2 must be the same')
    elseif (li > ni2)
        error('PGSS:Series:inputError','The length of INP1 must be lower than the number of inputs of G2')
    elseif (lo > no1)
        error('PGSS:Series:inputError','The length of OUT2 must be lower than the number of outputs of G1')
    elseif any(inp2 <= 0) || any(inp2 > ni2)
        error('PGSS:Series:inputError','INP2 out of range')
    elseif any(out1 <= 0) || any(out1 > no1)
        error('PGSS:Series:inputError','OUT2 out of range')
    end
    
    % Build series interconnection
    if isa(G1,'pgss')
        s.type = '()';
        s.subs = {out1,':',':'};
        G1 = subsref(G1,s);
    else
        G1 = G1(out1,:);
    end
    if isa(G2,'pgss')
        s.type = '()';
        s.subs = {':',inp2,':'};
        G2 = subsref(G2,s);
    else
        G2 = G2(:,inp2);
    end
    % resulting model
    Go = G2*G1;
    
else
    error('PGSS:Series:inputError','The number of input arguments must be 2 or 4');
end


