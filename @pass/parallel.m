function Go = parallel(G1,G2,inp1,inp2,out1,out2)

% PARALLEL connects LPV or LTI models in parallel
%
% Use:
%   G3 = PARALLEL(G1,G2,inp1,inp2,out1,out2)
%
% where:
%   G1, G2 are pass, ss, tf, or zpk objects
%   inp1,inp2: input indices of the resulting model
%   out1,out2: outputs to be summed
%   see lti.parallel
%
% See also pass, parallel, series

% fbianchi - 2021-03-31
% fbianchi - 2024-12-09 - rev

if (nargin == 2)
    
    % case with only models as input arguments
    
    % converting numeric arrays into ss
    if isnumeric(G1)
        G1 = ss(G1);
    end
    if isnumeric(G2)
        G2 = ss(G2);
    end
    
    % model dimensions
    [no1,ni1] = size(G1); ns1 = order(G1); 
    [no2,ni2] = size(G2); ns2 = order(G2); 
    if (no1 ~= no2) || (ni1 ~= ni2)
        error('PASS:parallel:inputError','All model must have the same number of inputs and output')

    else
        
        if isa(G1,'pass') && isa(G2,'pass')
            
            % case both models are LPV
            
            bool = isequal(G1.parset,G2.parset);
            if bool
                nv = nsys(G1);
                Ao(ns1+ns2,ns1+ns2,nv) = 0;
                Bo(ns1+ns2,ni1,nv) = 0;
                Co(no1,ns1+ns2,nv) = 0;
                Do(no1,ni1,nv) = 0;
                for ii = 1:nsys(G1)
                    Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), G2.A(:,:,ii));
                    Bo(:,:,ii) = [G1.B(:,:,ii); G2.B(:,:,ii)];
                    Co(:,:,ii) = [G1.C(:,:,ii), G2.C(:,:,ii)];
                    Do(:,:,ii) = G1.D(:,:,ii) + G2.D(:,:,ii);
                end
                
                pv  = G1.parset;
            else
                error('PASS:parallel:inputError','G1 or G2 must have the same parameter set')
            end
            
        elseif (isa(G1,'ss') || isa(G1,'zpk') || isa(G1,'tf')) && isa(G2,'pass')
            
            % case: G1 lti
            [A1,B1,C1,D1] = ssdata(ss(G1));
            
            Ao(:,:,1) = blkdiag(A1, G2.A(:,:,1));
            Bo(:,:,1) = [B1; G2.B(:,:,1)];
            Co(:,:,1) = [C1, G2.C(:,:,1)];
            Do(:,:,1) = G2.D(:,:,1) + D1;
            for ii = 2:nsys(G2)
                Ao(:,:,ii) = blkdiag(zeros(size(A1)), G2.A(:,:,ii));
                Bo(:,:,ii) = [zeros(size(B1)); G2.B(:,:,ii)];
                Co(:,:,ii) = [zeros(size(C1)), G2.C(:,:,ii)];
                Do(:,:,ii) = G2.D(:,:,ii);
            end
            
            pv  = G2.parset;
            
        elseif (isa(G2,'ss') || isa(G2,'zpk') || isa(G2,'tf')) && isa(G1,'pass')
            
            % case: G2 lti
            [A2,B2,C2,D2] = ssdata(ss(G2));
            
            Ao(:,:,1) = blkdiag(G1.A(:,:,1), A2);
            Bo(:,:,1) = [G1.B(:,:,1); B2];
            Co(:,:,1) = [G1.C(:,:,1), C2];
            Do(:,:,1) = D2 + G1.D(:,:,1);
            for ii = 2:nsys(G1)
                Ao(:,:,ii) = blkdiag(G1.A(:,:,ii), zeros(size(A2)));
                Bo(:,:,ii) = [G1.B(:,:,ii); zeros(size(B2))];
                Co(:,:,ii) = [G1.C(:,:,ii), zeros(size(C2))];
                Do(:,:,ii) = G1.D(:,:,ii);
            end
            
            pv  = G1.parset;
            
        else
            error('PASS:Parallel:inputError','G1 or G2 must be LTI models or pass objects')
            
        end
        
        % resulting model
        Go = pass(Ao, Bo, Co, Do, pv,...
                'StateName',[G1.StateName; G2.StateName],...
                'InputName',G1.InputName,...
                'OutputName',G2.OutputName);
    end
    
elseif (nargin == 6)
    
    % with input and output specifications

    [no1,ni1] = size(G1);
    [no2,ni2] = size(G2);
    
    % input validations
    li1 = length(inp1); inp1 = reshape(inp1,1,li1);
    li2 = length(inp2); inp2 = reshape(inp2,1,li2);
    lo1 = length(out1); out1 = reshape(out1,1,lo1);
    lo2 = length(out2); out2 = reshape(out2,1,lo2);
    if (li1 ~= li2)
        error('PASS:Parallel:inputError','IN1 and IN2 must have the same length');
    elseif (lo1 ~= lo2)
        error('PASS:Parallel:inputError','OUT1 and OUT2 must have the same length');
    elseif any(inp1 <= 0) || any(inp1 > ni1)
        error('PASS:Parallel:inputError','IN1 must be valid model inputs');
    elseif any(inp2 <= 0) || any(inp2 > ni2)
        error('PASS:Parallel:inputError','IN2 must be valid model inputs');
    elseif any(out1 <= 0) || any(out1 > no1)
        error('PASS:Parallel:inputError','OUT1 must be valid model outputs');
    elseif any(out2 <= 0) || any(out2 > no2)
        error('PASS:Parallel:inputError','OUT2 must be valid model outputs');
    end
    
    % Build parallel interconnection
    iv1 = 1:ni1;   iv1(inp1) = [];
    iz1 = 1:no1;   iz1(out1) = [];
    iv2 = 1:ni2;   iv2(inp2) = [];
    iz2 = 1:no2;   iz2(out2) = [];
    if isa(G1,'pass')
        s.type = '()';
        s.subs = {[iz1,out1],[iv1,inp1],':'};
        G1 = subsref(G1,s);
    else
        G1 = G1([iz1,out1],[iv1,inp1]);
    end
    if isa(G2,'pass')
        s.type = '()';
        s.subs = {[out2,iz2],[inp2,iv2],':'};
        G2 = subsref(G2,s);
    else
        G2 = G2([out2,iz2],[inp2,iv2]);
    end
    % Reevaluate sizes in case of repeated indices
    [no1,ni1] = iosize(G1);
    [no2,ni2] = iosize(G2);
    % resulting model
    Go = append(G1,zeros(no2-lo2,ni2-li2)) + append(zeros(no1-lo1,ni1-li1),G2);

else
    error('PASS:Parallel:inputError','The number of input arguments must be 2 or 6');
end