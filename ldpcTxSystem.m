function [ u_output ] = ldpcTxSystem( u_input, R, gammaDB)
%CONVOLUTIONALTXSYSTEM Simulates a 802.11 transmission system with ldpc encoding

mu = length(u_input);            % Input length

% Codeword length available 648, 1296, 1944 in the standard
% 648, 1296 are used only when there are few input bits
% 1944 will be used here

n = 1944;
k = n*R;

[H,G] = getHG(n,R);

sigmaw2 = 1/(10^(gammaDB/10));

%%%%%%% ENCODER %%%%%%%

nCW = ceil(mu/k);               % Number of codewords
Nshrt = nCW*k - mu;             % Number of padding bits
u_in = [ u_input zeros(1,Nshrt) ];
c = zeros(n,nCW);
G = double(G.x);
H = double(H.x);
for i=1:nCW
    payload = u_in((i-1)*k+1:i*k)';
    c(:,i) = mod(G*payload,2);
end

%sum(mod(H*y(:,i),2))

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

for i=1:nCW                 % For each codeword    
    llr = zeros(n,1);
    llr_est = llr;          % Each element is the LLR given from the sum of all the LLR sent 
                            % to the variable i from the corresponding check nodes 
    llr_int = llr;
    M = zeros(n-k,n);

    for j=1:length(llr)
        llr_int(j) = -2*r(j,i)/sigmaw2;      % Intrinsic information LLR
    end
    
    
    for ii=1:30               
        llr = llr_int + llr_est;        
        
        % Get the estimated output from the llr
        yCap = zeros(n,1);      
        yCap(llr<0) = 1;
        
        
        % If the estimated codeword satisfies
        % all the check rules break the cycle        
        if(sum(mod(H*yCap,2))==0)                         
            break;
        end
    
        % For each check node compute the messages to the variable nodes        
        for cn=1:(n-k)
            idx = find(H(cn,:));   % Find the indexes of the variable nodes involved in this check
            for vn=idx   
                %M(cn,vn) = cnMess(vn,cn,idx,llr);                                
                inLLR = llr(idx,1);
                inLLR(idx==vn) = [];
                M(cn,vn) = prod(sign(inLLR))*(lntanh(sum(lntanh(abs(inLLR)))));
            end            
        end    
        
        % Bound the LLR from the check nodes to the variable nodes to prevent numerical instability
        M(M>40) = 200;
        M(M<-40) = -200;
        
        % For each variable nodes compute the extrinsing information
        for vn=1:n
            idx = find(H(:,vn));  
            llr_est(vn) = sum(M(idx,vn));
            %llr_est(vn) = vnMess(vn,idx,M);        
        end
       
    end
    % Finally extract the payload from the codeword
    % quite easy since the code is in systematic form
    u_out((i-1)*k+1:i*k) = yCap(1:k,1)';    
    
end
% Remove the padding bits
u_output = u_out(1:mu);
end

