function ns = order(obj)

% ORDER(pdG) returns the order of the lpv model pdG

% fbianchi - 2020-06-28

ns = size(obj.ctrller.A(:,:,1),1);
      
