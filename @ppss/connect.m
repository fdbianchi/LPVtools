function pdGo = connect(varargin)

% CONNECT connects the block diagram elements based on i/o names
%
% Use:
%   sysc = CONNECT(sys1,...,sysN,inputs,outputs)
%
% where ONE system can be a ppss object.
%
% See also ppss, connect

% fbianchi - 2021-03-31

% -------------------------------------------------------------------------
% check if there is more than 1 lpv system to interconnect
types = cellfun(@class, varargin, 'uni', false);
typePass = ismember(types,'ppss');  
if (sum(typePass) > 1)
    error('Only 1 ppss object can be interconnected');
end

% ppss plant
pdG = varargin{typePass};

% to use matlab functions first the LPV model is transformed into ss
% object and after using connect transformed again into ppss object
%
% convertion into ss object
varargin{typePass} = ss(pdG);

% to avoid simplifications in each system
opt = connectOptions('Simplify',false);
varargin{end+1} = opt;

% interconnection
sys_tmp = connect(varargin{:});

% convertion into ppss object
warning('off','Control:ltiobject:RepeatedChannelNames')
pdGo = ppss(sys_tmp,pdG.parset);
warning('on','Control:ltiobject:RepeatedChannelNames')
