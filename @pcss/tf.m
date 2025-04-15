function sys = tf(pdG,pval)

% TF converts a pcss object into a 3-D tf object
%
% Use:
%   G = TF(pdG)
%   G = TF(pdG,pval)
%
% where 
%   - pdG:  lpv model
%   - pval: a matrix of np x nv with the parameter values at which the
%           lpv model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - G:    transfer functions corresponding to the frozen parameter values
%           (tf object)
%
% TF(pdG) evaluates an LPV model pdG at frozen parameter values 
%   and produces a transfer functions G
%
% See also pass, ppss, pgss, pcss, tf

% fbianchi - 2020-07-02

if (nargin > 1)
    % the system is evaluated at pval
    sys = tf(ss(pdG,pval));

else
    % the system is evaluated at the points in the parameter set
    sys = tf(ss(pdG));
end

