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
       
       
       function bool = isempty(obj)
           
           nro_const = length(obj);
           
           if (nro_const == 0)
               bool = true;
               
           elseif (nro_const == 1)
               
               if ismember(class(obj),{'synConst.Gain','synConst.GainH2',...
                       'synConst.GainH2g'})
                   
                   bool = isempty(obj.inputs) || isempty(obj.outputs);
                   
               elseif isa(obj,'synConst.Poles')
                   
                   bool = (obj.MinDecay == 0) && (obj.MinDamping == 0) && isinf(obj.MaxFreq);
               end
               
           else
               bool = false;
               
           end
       end
       
       
       function disp(obj)
           
           nro_const = length(obj);
           
           if (nro_const == 1)
               disp(char(obj))
           else
               for ii = 1:length(obj)
                   str = char(obj(ii));
                   fprintf('%2.0f) %s',ii,str(2:end));
               end
           end
       end
   end

end