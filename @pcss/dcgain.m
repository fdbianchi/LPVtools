function varargout = dcgain(obj,pval)

% DCGAIN computes the DC gain of an lpv model at frozen parameter values
%
% Use:
%   DCGAIN(pdG,[pval])
%   g = DCGAIN(pdG,[pval])
%
% where 
%   - pdG: lpv model
%   - pval: a matrix of np x nv with the parameter values in which the
%           lpv model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - g:    vector with the DC gains
%
% See also pass, ppss, pgss, pcss, dcgain

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
    
% computing the DC gain
g  = dcgain(G);
nv = size(pval,2);

% displaying values when no output argument
if (nargout == 0)
    [no,ni] = iosize(G);
    for ii =  1:nv
        if (nargin < 2)
            strpar = sprintf('%g, ',obj.ctrller.parset.points(:,ii));
        else
            strpar = sprintf('%g, ',pval(:,ii));
        end
        
        fprintf('DC gain at [%s]\n',strpar(1:end-2));
        for jj = 1:no
            for kk = 1:ni
                fprintf('\t%8.4f',g(jj,kk,ii))
            end
            fprintf('\n')
        end
        fprintf('\n')
    end
    varargout = {};
    
else
    varargout = {g};
end
 