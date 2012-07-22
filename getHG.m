function [ H, G ] = getHG( n, R )
%BUILDH Creates the parity check matrix H reading the prototype from file for
% a codeword length n and rate R

Zsize = [27, 54, 81];           % Square submatrices available size

col = 24;       % Number of submatrices in a row

if R == 1/2
    row = 12;       % Number of submatrices in a column
    rStr = '12';
elseif R == 2/3
    row= 8;
    rStr = '23';
elseif R == 3/4
    row = 6;
    rStr = '34';
else
    row = 4;
    rStr = '56';
end

% This block will be needed only if the support for different
% codeword size will be implemented 
if n == 648         
    Z = Zsize(1);    
elseif n == 1296
    Z = Zsize(2);
else
    Z = Zsize(3);
end

if(exist(['Matrix/',num2str(n),'_',rStr,'.mat'],'file'))
    load(['Matrix/',num2str(n),'_',rStr],'H');
    load(['Matrix/',num2str(n),'_',rStr],'G');
else
    H = zeros(row*Z,col*Z);

    % Prototype matrix
    % If the element (i,j) is >= 0 the corresponding submatrix will be a column circular
    % shift of the identy matrix by the number of position indicated by the element
    % otherwise the submatrix will be a zero matrix

    protoH = load(['protoH/',num2str(n),'_',rStr]);

    for i = 1:row
        for j = 1:col
            if protoH(i,j) >= 0
                A = eye(Z);
                k = protoH(i,j);
                H((i-1)*Z+1:(i-1)*Z+Z,(j-1)*Z+1:(j-1)*Z+Z) = circshift(A,[0 k]);
            end        
        end
    end

    
    H = gf(H);
    k = n*R;
    B = H(:,1:n-k);
    C = H(:,n-k+1:n);
    C1 = inv(C);
    G = [eye(k);C1*B];
    
    CHECK = H*G;
    CHECK = CHECK.x;
    if sum(CHECK(:))
        disp('G is not correct');
        return
    end
    
    mkdir('Matrix');
    save(['Matrix/',num2str(n),'_',rStr],'H');
    save(['Matrix/',num2str(n),'_',rStr],'G','-append');
end
G = double(G.x);
H = double(H.x);

end

