function [ u_output, ber, fer ] = ldpcTxSystemFast( u_input, R, gammaDB, backSubstitution, varargin)
%LDPCTXSYSTEMFAST Simulates a 802.11 transmission system with ldpc encoding

if length(varargin)
    iterations = varargin{1};
else
    iterations = 25;
end

mu = length(u_input);            % Input length

% Codeword length available in the standard: 648, 1296, 1944
% 648, 1296 are used only when there are few input bits
% 1944 used in all other cases (will be used here)

n = 1944;   % Codeword length
k = n*R;    % Payload length

%%% PARITY CHECK MATRIX %%%%

if backSubstitution
    [H,Z] = buildH(n,R);
else
    [H,G,Z] = buildHG(n,R);
end

% To improve performance instead of using cell, a special matrix is
% used. The first element of each row contains the number of useful element
% in each row starting from index 2
A = zeros(n,23);          % Row i contains the indexes of the check nodes in which the variable i is involved
B = zeros(n-k,13);        % Row i contains the indexes of the variable nodes involved in the check i

for vn=1:n
    row = find(H(:,vn));
    len = length(row);
    A(vn,1) = len;
    A(vn,2:len+1) = row;
end
for cn=1:(n-k)
    row = find(H(cn,:));
    len = length(row);
    B(cn,1) = len;
    B(cn,2:len+1) = row;
end


%%%%%%% ENCODER %%%%%%%

nCW = ceil(mu/k);                   % Number of codewords
Npad = nCW*k - mu;                  % Number of padding bits
u_in = [ u_input zeros(1,Npad) ];

if backSubstitution         % Encoding by backsubstitution
    c = zeros(n,nCW);
    H1 = H(:,1:k);
    H2 = H(:,(k+1):(k+Z));    
    
    for i=1:nCW
        payload = u_in((i-1)*k+1:i*k)';        
        s = H1 * payload;
        s = reshape(s,Z,(n-k)/Z);
        p1 = mod(sum(s,2),2);
        temp = H2 * p1;
        temp = reshape(temp,Z,(n-k)/Z);
        stilda = mod(s + temp,2);
        prevp2 = 0;
        ll = length(stilda(1,:)) - 1;
        p2 = zeros(ll*Z,1);
        for j=1:ll
            prevp2 = mod(stilda(:,j) + prevp2,2);
            p2((j-1)*Z+1:j*Z) = prevp2;            
        end
        c(:,i) = [ payload; p1; p2];
    end    
else        % Enconding by matrix product
    c = mod(G*reshape(u_in,k,length(u_in)/k),2);    
end


%%%%%%%%% P/S %%%%%%%%%

d = c(:)';


%%% BPAM MODULATOR %%%%

sTx = d;
sTx(sTx==0) = -1;    % Map 0 to -1 and produce the transmitted signal


%%%%%%% CHANNEL %%%%%%%

r = awgn(sTx,gammaDB);          % Received signal


%%%%%%%%% S/P %%%%%%%%%

r = reshape(r,n,length(r)/n);


%%% MESSAGE PASSING DECODER %%%

sigmaw2 = 1/(10^(gammaDB/10));      % Noise variance
        
[u_out, checkOK] = mexfastdecoder(k,n,nCW,sigmaw2,A,B,H,2*r/sigmaw2,iterations);

% Remove the padding bits
u_output = u_out(1:mu);

% Compute bit error rate and frame error rate
ber = sum(u_output~=u_input)/mu;
fer = 1 - checkOK/nCW;
end


