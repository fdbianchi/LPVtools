classdef synConst < matlab.mixin.Heterogeneous
    
    % SYNCONST is an heterogeneous super-class to contain several constraints
    %
    % See also synCont.Gain, synCont.Poles, synCont.GainH2
    
    % fbianchi - 2021-06-30
    

   methods (Static, Sealed, Access = protected)
      function defaultObject = getDefaultScalarElement
         defaultObject = synConst.Gain;
      end
   end
   
   methods (Sealed)
        
%         function varargout = char(obj)
%             
%             % CHAR method for synConst class
%             
%             if (nargout < 2)
%                 varargout{1} = [];
%             else
%                 varargout = {[],[],[]};
%             end
%             
%             varargout = 'SynConst object';
%             
%         end
        
        function disp(obj)
            
            nro_const = length(obj);
            
            if (nro_const == 1)
                disp(char(obj))
            else
            for ii = 1:length(obj)
%                 str = sprintf(
                str = char(obj(ii));
                fprintf('%2.0f) %s',ii,str(2:end));
                
%                 if ismember(class(obj(ii)),{'synConst.Gain','synConst.GainH2',...
%                         'synConst.GainH2g','synConst.Poles'})
%                     disp(char(obj(ii)))
%                 else
%                     disp('Empty synConst object');
%                 end
            end
            end
        end
    end
end