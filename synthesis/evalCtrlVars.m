function [ak,bk,ck,dk,E,F] = evalCtrlVars(ctrlSet,idx,pdG)

% *** for internal use ***
%
% EVALCTRLVARS returns a set of sdpvar obj corresponding to the controller
% variables
%
% Use: 
%   [ak,bk,ck,dk,E,F] = evalCtrlVars(ctrlSet,idx,pdG)
%
% where:
%   - CtrlSet:  struct with controller variables
%   - idx:      optimization variable index (for ppss/PWA)
%   - pdG:      plant model, only used to find E and F in case of
%               controller simplifications

% fbianchi - 2020-04-24

E = []; F = [];
par = ctrlSet.parset.points(:,idx);

switch ctrlSet.interp

    case 'cte'
        % case robust
        ak = ctrlSet.A;
        bk = ctrlSet.B;
        ck = ctrlSet.C;
        dk = ctrlSet.D;
    
    case 'pwa'
        % case polytopic
        if strcmp(ctrlSet.fback,'state')
            ak = [];
            bk = [];
            ck = [];
            dk = ctrlSet.D(:,:,idx);
        else
            ak = ctrlSet.A(:,:,idx);
            bk = ctrlSet.B(:,:,idx);
            ck = ctrlSet.C(:,:,idx);
            dk = ctrlSet.D(:,:,idx);
        end
    
    case 'aff'
        % affine o general dependency
        if isa(ctrlSet.parfcn,'function_handle')
            fpar = [1; ctrlSet.parfcn(par)];
        else
            fpar = [1; par];
        end
    
        % in case of given controller parameter dependence, only the terms in
        % CtrlSet.idxfcn are considered
        Iset = ctrlSet.idxfcn;
        if (ctrlSet.extFcn == 1)
            Jset = setdiff(1:length(fpar),ctrlSet.idxfcn);
            % part of A not included in Ak
            Ar = 0;
            for ii = Jset
                Ar = Ar + pdG.A(:,:,ii)*fpar(ii);
            end
            ra = rank(Ar);
            if (ra > 0)
                % Ar factorization
                [U,Sv,V] = svds(Ar,ra);
                E = U*sqrt(Sv); F = V*sqrt(Sv);
            end
        end
        
        if ~isempty(ctrlSet.A)
            ak = 0;
            bk = 0;
            ck = 0;
            dk = 0;
            for ii = 1:length(Iset)
                ak = ak + ctrlSet.A(:,:,ii)*fpar(Iset(ii));
                bk = bk + ctrlSet.B(:,:,ii)*fpar(Iset(ii));
                ck = ck + ctrlSet.C(:,:,ii)*fpar(Iset(ii));
                dk = dk + ctrlSet.D(:,:,ii)*fpar(Iset(ii));
            end
            
        elseif ~isempty(ctrlSet.B)
            ak = [];
            bk = 0;
            ck = 0;
            dk = 0;
            for ii = 1:length(Iset)
                bk = bk + ctrlSet.B(:,:,ii)*fpar(Iset(ii));
                ck = ck + ctrlSet.D(:,:,ii)*fpar(Iset(ii));
                dk = dk + ctrlSet.C(:,:,ii)*fpar(Iset(ii));
            end
            
        else
            ak = [];
            bk = [];
            ck = [];
            dk = 0;
            for ii = 1:length(Iset)
                dk = dk + ctrlSet.D(:,:,ii)*fpar(Iset(ii));
            end
            
        end
        
    otherwise
        error('EVALCTRLVARS','Invalid controller matrix interpolation');
        
end

