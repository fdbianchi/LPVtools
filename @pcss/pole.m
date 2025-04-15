function varargout = pole(varargin)

% POLE computes the poles of the LTI models corresponding to an LPV model
%   evaluated at frozen parameter values
%
% Use:
%   POLE(pdG,[pval])
%   p = POLE(pdG,[pval])
%
% where
%   - pdG:  lpv system
%   - pval: a matrix of np x nv with the parameter values at which the
%           lpv model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - p:    vector with the eigenvalues
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2020-07-02

if (nargout == 0)
    eig(varargin{:});
    
else
    varargout{1} = eig(varargin{:});
    
end
