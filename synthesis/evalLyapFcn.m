function [y,x,dy,dx,Yaff,Xaff] = evalLyapFcn(lyapSet, idx, simplex)

% *** for internal use ***
%
% EVALLYAPFCN returns a set of sdpvar obj corresponding to the Lyapunov
% functions
%
% Use: 
%   [y,x,dy,dx] = evalLyapFcn(lyapSet, idx, simplex)
%
% where:
%   - lyapSet:  struct with Lyapunov function variables
%   - idx:      optimization variable index (for ppss/PWA)
%               point index in affine case
%   - simplex:  vertices of the simplex in which variable corresponding to
%               the index idx is located. Only used for PWA and Affine 
%               parameter dependent Lyapunov functions

% fbianchi - 2020-04-24

Yaff = [];
Xaff = [];
dy = {0};
dx = {0};

switch lyapSet.type(1:3)
    
    case 'cte'
        % constant
        y = lyapSet.Y;
        x = lyapSet.X;
        
    case 'pwa'
        % PWA case
        if (size(lyapSet.Y,3) > 1)
            y = lyapSet.Y(:,:,idx);
        else
            y = lyapSet.Y;
        end
        if (size(lyapSet.X,3) > 1)
            x = lyapSet.X(:,:,idx);
        else
            x = lyapSet.X;
        end
        
        if (nargout > 2)
            
            ns = size(y,1);
            
            % local affine description 
            nv = length(simplex);
            Points = lyapSet.parset.points;
            
            % central point
            pc = Points(:,simplex)*ones(nv,1)/nv;
            Delta = [ones(1,nv);
                     Points(:,simplex) - pc*ones(1,nv)];

            if (size(lyapSet.Y,3) > 1)
                aux = reshape(lyapSet.Y(:,:,simplex),ns*ns,nv);
                Yaff = reshape(aux/Delta,ns,ns,nv);
                
                if any(strcmp(lyapSet.type,'pwadY'))
                    dgrd = pgrid(lyapSet.parset.rate);  % hypercube corresponding
                    % to the rate
                    for ii = 1:size(dgrd,2)
                        auxdy = 0;
                        for jj = 1:size(dgrd,1)
                            if ~isinf(dgrd(jj,ii))
                                auxdy = auxdy + Yaff(:,:,jj+1)*dgrd(jj,ii);
                            end
                        end
                        % dy(:,:,ii) = auxdy; % not supported
                        dy{ii} = auxdy;
                    end
                end
            end
            
            if (size(lyapSet.X,3) > 1)
                aux = reshape(lyapSet.X(:,:,simplex),ns*ns,nv);
                Xaff = reshape(aux/Delta,ns,ns,nv);
                
                if any(strcmp(lyapSet.type,'pwadX'))
                    dgrd = pgrid(lyapSet.parset.rate);  % hypercube corresponding
                    % to the rate
                    for ii = 1:size(dgrd,2)
                        auxdx = 0;
                        for jj = 1:size(dgrd,1)
                            if ~isinf(dgrd(jj,ii))
                                auxdx = auxdx + Xaff(:,:,jj+1)*dgrd(jj,ii);
                            end
                        end
                        % dx(:,:,ii) = auxdx; % not supported
                        dx{ii} = auxdx;
                    end
                end
            end
        end
        
        
        
    otherwise
        
        par = lyapSet.parset.points(:,idx);

        % affine or general
        if isempty(lyapSet.X)
            % state feedback
            x  = [];
        else
            % output feedback
            x = lyapSet.X(:,:,1);
            if isa(lyapSet.parfcnX,'function_handle')
                fparX = lyapSet.parfcnX(par);
                for ii = 1:length(fparX)
                    x = x + lyapSet.X(:,:,ii+1)*fparX(ii);
                end
            end
            % compute derivatives for each vertex in rate
            if isa(lyapSet.dparfcnX,'function_handle') && any(any(lyapSet.parset.rate))
                % removing inf parameter rates
                aux = lyapSet.parset.rate;
                idx = isinf(aux(:,1)) | isinf(aux(:,2));
                aux(idx,:) = 0;
                % hypercube corresponding to the rate
                dgrd = pgrid(aux);
                for ii = 1:size(dgrd,2)
                    try
                        dfparX = lyapSet.dparfcnX(par)*dgrd(:,ii);
                    catch ME
                        error('EVALLYAPFCN:inputError',...
                              'dparfcnX is incompatible with the parameter rate definition')
                    end
                    auxdx = 0;
                    for jj = 1:length(dfparX)
                        auxdx = auxdx + lyapSet.X(:,:,jj+1)*dfparX(jj);
                    end
                    dx{ii} = auxdx;
                end
            end
        end
        
        y = lyapSet.Y(:,:,1);
        if isa(lyapSet.parfcnY,'function_handle')
            fparY = lyapSet.parfcnY(par);
            for ii = 1:length(fparY)
                y = y + lyapSet.Y(:,:,ii+1)*fparY(ii);
            end
        end
        % compute derivatives for each vertex in rate
        if isa(lyapSet.dparfcnY,'function_handle') && any(any(lyapSet.parset.rate))
            % removing inf parameter rates
            aux = lyapSet.parset.rate;
            idx = isinf(aux(:,1)) | isinf(aux(:,2));
            aux(idx,:) = 0;
            % hypercube corresponding to the rate
            dgrd = pgrid(aux);%,2*(~idx) + 1*(idx));%lyapSet.parset.rate);
            for ii = 1:size(dgrd,2)
                    try
                        dfparY = lyapSet.dparfcnY(par)*dgrd(:,ii);
                    catch ME
                        error('EVALLYAPFCN:inputError',...
                              'dparfcnY is incompatible with the parameter rate definition')
                    end
                auxdy = 0;
                for jj = 1:length(dfparY)
                    auxdy = auxdy + lyapSet.Y(:,:,jj+1)*dfparY(jj);
                end
                % dy(:,:,ii) = auxdy; % not supported
                dy{ii} = auxdy;
            end
        end
end

