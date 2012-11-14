function [ H, G, Z ] = buildHG( n, R )
%BUILDHG Creates the parity check matrix H and the generator matrix G
%(reading the prototype from file) for codeword length n and rate R

[H, Z] = buildH(n, R);

if R == 1/2
    rStr = '12';
elseif R == 2/3
    rStr = '23';
elseif R == 3/4
    rStr = '34';
else
    rStr = '56';
end

if(exist(['matrix/',num2str(n),'_',rStr,'.mat'],'file'))
    load(['matrix/',num2str(n),'_',rStr],'G');
else
    gfH = gf(H);
    k = n*R;
    B = gfH(:,1:n-k);
    C = gfH(:,n-k+1:n);
    C1 = inv(C);            % Very slow operation
    G = [eye(k);C1*B];

    CHECK = gfH*G;
    if any(CHECK(:))
        disp('G is not correct');
        return
    end

    if ~exist('matrix','dir')
        mkdir('matrix');
    end
    save(['matrix/',num2str(n),'_',rStr],'G','-append');
end
G = double(G.x);

end

