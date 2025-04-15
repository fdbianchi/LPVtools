function varargout = tzero(obj,pval)

% TZERO computes the zeros of the LTI model corresponding to an LPV model 
%   evaluated at frozen parameter values
%
% Use:
%   TZERO(pdG,[pval])
%   z = TZERO(pdG,[pval])
%
% where
%   - pdG:  LPV model
%   - pval: a matrix of (np x nv) with the parameter values in which the
%           LPV model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - z:    vector of transmission zeros 
%
% See also pass, ppss, pgss, pcss, tzero

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

% compute zeros for each LTI model
nv = size(pval,2);
for ii = 1:nv
    z{ii} = tzero(G(:,:,ii));
end

% displaying values when no output argument
if (nargout == 0)
    for ii = 1:nv
        if (nargin < 2)
            strpar = sprintf('%g, ',obj.parset.points(:,ii));
        else
            strpar = sprintf('%g, ',pval(:,ii));
        end
        
        fprintf('Transmission zeros at [%s]\n',strpar(1:end-2));
        
        if isempty(z{ii})
            fprintf('\tNo transmission zeros\n')
        else
            if all(isreal(z{ii}))
                fprintf('\t %8.4f\n',z{ii})
                
            else
                fprintf('\t %8.4f %+8.4fi\n',[real(z{ii}),imag(z{ii})]')
            end
            
            fprintf('\n')
        end
    end
    varargout = {};
    
else
    varargout = {z};
end 
