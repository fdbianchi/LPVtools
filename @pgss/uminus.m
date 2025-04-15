function Go = uminus(Gi)

% UMINUS unary minus for LPV models.
%
%   MM = UMINUS(M) is invoked by MM = -M.
%
%
% See also pgss

% fbianchi - 2020-02-20

Go = pgss(Gi.A, Gi.B, -Gi.C, -Gi.D, Gi.parset, Gi.parfcn,...
        'StateName',Gi.StateName,...
        'InputName',Gi.InputName,...
        'OutputName',Gi.OutputName,...
        'InputGroup',Gi.InputGroup,...
        'OutputGroup',Gi.OutputGroup);

