function C = matProd(A,B,nt,int)

if all(A(:,:,2:end) == 0,'all')
    % A is constant (case affine)
    na = 1;
else
    auxA = num2cell(A,[1 2]);
    if isequal(auxA{:})
        % A is constant (case pwa)
        na = 1;
    else
        na = size(A,3);
    end
end

if all(B(:,:,2:end) == 0,'all')
    % B is constant (case affine)
    nb = 1;
else
    auxB = num2cell(B,[1 2]);
    if isequal(auxB{:})
        % B is constant (case pwa)
        nb = 1;
    else
        nb = size(B,3);
    end
end

if (na > 1) && (nb > 1) && (na ~= nb)
    error('matProd:inputError','the third dimension of A and B must be equal')
end
    
nr = size(A,1);     % # of rows
nc = size(B,2);     % # of columns
m = max(na,nb);     % # of individual terms

if (na == 1) && (nb == 1)
    % result: constant matrix
    C(:,:,1) = A(:,:,1)*B(:,:,1);
    if (nargin > 2) 
        if strcmp(int,'aff')
            C = cat(3,C,zeros(nr,nc,nt-1));
        else
            C = repmat(C,1,1,nt);
        end
    end
    
elseif (na == 1) || (nb == 1)
    % result: affine
    C = zeros(nr,nc,m);
    if (na == 1)
        for ii = 1:m
            C(:,:,ii) = A(:,:,1)*B(:,:,ii);
        end
    else
        for ii = 1:m
            C(:,:,ii) = A(:,:,ii)*B(:,:,1);
        end
    end        

else
    % result: multi-affine
    nt = 2*m - 1 + nchoosek(m-1,2);   % # of terms
    C = zeros(nr,nc,nt);
    % constant term
    C(:,:,1) = A(:,:,1)*B(:,:,1);
    % affine terms
    for ii = 2:m
        C(:,:,ii) = A(:,:,1)*B(:,:,ii) + A(:,:,ii)*B(:,:,1);
    end
    % multi-affine terms
    idx = nchoosek(2:m,2);
    for jj = 1:size(idx,1)
        ii = ii + 1;
        C(:,:,ii) = A(:,:,idx(jj,1))*B(:,:,idx(jj,2)) + A(:,:,idx(jj,2))*B(:,:,idx(jj,1));
    end
    % quadratic terms
    for jj = 2:m
        ii = ii + 1;
        C(:,:,ii) = A(:,:,jj)*B(:,:,jj);
    end
    
end


