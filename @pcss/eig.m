function varargout = eig(obj,pval)

% EIG computes the poles of an lpv model at frozen parameter values
%
% Use:
%   EIG(pdG,[pval])
%   e = EIG(pdG,[pval])
%
% where
%   - pdG: lpv model
%   - pval: a matrix of np x nv with the parameter values in which the
%           lpv model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - p:    vector with the eigenvalues

% fbianchi - 2020-07-02

% if pval is not given, the dc gain is computed for the point in pset
if (nargin < 2)
    % all points in parset
    pval = obj.ctrller.parset.points;
end
    
% the pcss object is transformed into multidimensional ss obj
[G,msg] = subs(obj,pval);
if isempty(G)
    error(msg)
end
    
% computing the eigenvalues
e  = eig(G);
nv = size(pval,2);

% displaying values when no output argument
if (nargout == 0)
    ns = order(obj);
    e  = reshape(e,ns,nv);
    for ii =  1:nv
        if (nargin < 2)
            strpar = sprintf('%g, ',obj.ctrller.parset.points(:,ii));
        else
            strpar = sprintf('%g, ',pval(:,ii));
        end
        
        fprintf('Eigenvalues at [%s]\n',strpar(1:end-2));
        
        for jj = 1:length(e(:,ii))
            if isreal(e(jj,ii))
                fprintf('\t %+9d\n',e(jj,ii))
            else
                fprintf('\t %+9d %+9di\n',real(e(jj,ii)),imag(e(jj,ii)))
            end
        end

        fprintf('\n')
    end
    varargout = {e};
    
else
    varargout = {e};
end
