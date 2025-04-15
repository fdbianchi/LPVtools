function varargout = eig(obj,pval)

% EIG computes the poles of the LTI models corresponding to an LPV model
%   evaluated at frozen parameter values
%
% Use:
%   EIG(pdG,[pval])
%   e = EIG(pdG,[pval])
%
% where
%   - pdG:  LPV model
%   - pval: a matrix of (np x nv) with the parameter values at which the
%           LPV model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - e:    vector with the eigenvalues
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31

if (nargin < 2)
    % all points in parset
    pval = obj.parset.points;
end
    
% the p_ss object is transformed into multidimensional ss obj
[G,msg] = subs(obj,pval);
if isempty(G)
    error(msg)
end
    
% computing the eigenvalues
e = eig(G);
nv = size(pval,2);

% displaying values when no output argument
if (nargout == 0)
    ns = order(obj);
    e  = reshape(e,ns,nv);
    for ii =  1:nv
        if (nargin < 2)
            strpar = sprintf('%g, ',obj.parset.points(:,ii));
        else
            strpar = sprintf('%g, ',pval(:,ii));
        end
        
        fprintf('Eigenvalues at [%s]\n',strpar(1:end-2));
        
        for jj = 1:length(e(:,ii))
            if isreal(e(jj,ii))
                fprintf('\t %+.5g\n',e(jj,ii))
            else
                if (sign(imag(e(jj,ii))) == -1)
                    sgn = '-';
                else
                    sgn = '+';
                end
                fprintf('\t %+.5g %s %.5gi\n',real(e(jj,ii)),sgn,abs(imag(e(jj,ii))))
            end
        end

        fprintf('\n')
    end
    varargout = {e};
    
else
    varargout = {e};
end
