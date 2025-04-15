function [x,dx,Xaff] = evalLyapFcn_a(lyapSet, idx, simplex)

% *** for internal use ***
%
% EVALLYAPFCN_A returns a set of sdpvar obj corresponding to the Lyapunov
% functions
%
% Use: 
%   [x,dx,Xaff] = evalLyapFcn_a(lyapSet, par, idx)
%
% where:
%   - lyapSet:  struct with Lyapunov function variables
%   - par:      parameter value
%   - idx:      optimization variable index (for ppss/PWA)

% fbianchi - 2020-04-24


Xaff = [];
dx = {0};

switch lyapSet.type(1:3)
    
    case 'cte'
        % constant
        x  = lyapSet.X;
        
    case 'pwa'
        % PWA case
        if (size(lyapSet.X,3) > 1)
            x = lyapSet.X(:,:,idx);
        else
            x = lyapSet.X;
        end
        
        if (nargout > 1)
            % local affine transformation
            
            ns = size(x,1);
            
            % local affine descprition
            nv = length(simplex);
            Points = lyapSet.parset.points;
            
            % central point
            pc = Points(:,simplex)*ones(nv,1)/nv;
            Delta = [ones(1,nv);
                Points(:,simplex) - pc*ones(1,nv)];
            
            if (size(lyapSet.X,3) > 1)
                aux = reshape(lyapSet.X(:,:,simplex),ns*ns,nv);
                Xaff = reshape(aux/Delta,ns,ns,nv);
                
                if strcmp(lyapSet.type,'pwadX') && ~isempty(lyapSet.parset.rate) && any(any(lyapSet.parset.rate))
                    % removing inf parameter rates
                    aux = lyapSet.parset.rate;
                    idx = isinf(aux(:,1)) | isinf(aux(:,2));
                    aux(idx,:) = 0;
                    % hypercube corresponding to the rate
                    dgrd = pgrid(aux);
                    for ii = 1:size(dgrd,2)
                        auxdx = 0;
                        for jj = 1:size(dgrd,1)
                            auxdx = auxdx + Xaff(:,:,jj+1)*dgrd(jj,ii);
                        end
                        dx{ii} = auxdx;
                    end
                end
            end
        end
        
    otherwise
        % affine or general
        
        par = lyapSet.points(:,idx);
        
        x = lyapSet.X(:,:,1);
        if isa(lyapSet.parfcn,'function_handle')
            fpar = lyapSet.parfcn(par);
            for ii = 1:length(fpar)
                x = x + lyapSet.X(:,:,ii+1)*fpar(ii);
            end
        end
        % compute derivatives for each vertex in rate
        if isa(lyapSet.dparfcn,'function_handle') && ~isempty(lyapSet.parset.rate)
            % removing inf parameter rates
            aux = lyapSet.parset.rate;
            idx = isinf(aux(:,1)) | isinf(aux(:,2));
            aux(idx,:) = 0;
            % hypercube corresponding to the rate
            dgrd = pgrid(aux);
            for ii = 1:size(dgrd,2)
                try
                    dfpar = lyapSet.dparfcn(par)*dgrd(:,ii);
                catch ME
                    error('EVALLYAPFCN:inputError',...
                        'dparfcnX is incompatible with the parameter rate definition')
                end
                auxdx = 0;
                for jj = 1:length(dfpar)
                    auxdx = auxdx + lyapSet.X(:,:,jj+1)*dfpar(jj);
                end
                dx{ii} = auxdx;
            end
        end
end

















% switch lyapSet.type
%     
%     case 'cte'
%         % constant
%         x  = lyapSet.X;
%         dx = {0};
%         
%     case 'pwa'
%         % ToDo: check case PWA
% %         y = lyapSet.Y(:,:,idx);
%         x = lyapSet.X(:,:,idx);
%         dx = {0};
% %         dy = {0};
%         
%     otherwise
%         
% %         par = lyapSet.parset.points(:,idx);
%         par = lyapSet.points(:,idx);
% 
%         % output feedback
%         x = lyapSet.X(:,:,1);
%         if isa(lyapSet.parfcn,'function_handle')
%             fpar = lyapSet.parfcn(par);
%             for ii = 1:length(fpar)
%                 x = x + lyapSet.X(:,:,ii+1)*fpar(ii);
%             end
%         end
%         % compute derivatives for each vertex in rate
%         if isa(lyapSet.dparfcn,'function_handle')
%             dgrd = pgrid(lyapSet.parset.rate);  % hypercube corresponding
%             % to the rate
%             for ii = 1:size(dgrd,2)
%                 dfpar = lyapSet.dparfcn(par)*dgrd(:,ii);
%                 auxd = 0;
%                 for jj = 1:length(dfpar)
%                     auxd = auxd + lyapSet.X(:,:,jj+1)*dfpar(jj);
%                 end
%                 dx{ii} = auxd;
%             end
%         else
%             dx = {0};
%         end
%         
% end

