function [ u_output ] = ldpcTxSystem( u_input, R, gammaDB, mexEnabled, backSubstitution)
%CONVOLUTIONALTXSYSTEM Simulates a 802.11 transmission system with ldpc encoding

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


A = cell(1,n);          % Element i contains the indexes of the check nodes in which the variable i is involved
B = cell(1,n-k);        % Element i contains the indexes of the variable nodes involved in the check i

for vn=1:n
    A{vn} = find(H(:,vn))';
end
for cn=1:(n-k)
    B{cn} = find(H(cn,:));
end



%%%%%%% ENCODER %%%%%%%

nCW = ceil(mu/k);                   % Number of codewords
Npad = nCW*k - mu;                  % Number of padding bits
u_in = [ u_input zeros(1,Npad) ];

if backSubstitution
    c = zeros(n,nCW);
    H1 = H(:,1:12*Z);
    H2_1 = H(:,12*Z+1:13*Z);    
    tic;
    for i=1:nCW
        payload = u_in((i-1)*k+1:i*k)';        
        s = H1 * payload;
        s = reshape(s,Z,(n-k)/Z);
        p1 = mod(sum(s,2),2);
        temp = H2_1 * p1;
        temp = reshape(temp,Z,(n-k)/Z);
        stilda = mod(s + temp,2);
        prevp2 = 0;
        p2 = zeros(11*Z,1);
        for j=1:11
            prevp2 = mod(stilda(:,j) + prevp2,2);
            p2((j-1)*Z+1:j*Z) = prevp2;            
        end
        c(:,i) = [ payload; p1; p2];
    end
    cc = toc;
else    
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
a = 0;
b = 0;
d = 0;

u_out = zeros(1,length(u_in));
sigmaw2 = 1/(10^(gammaDB/10));

L = zeros(1,n);         % Each element is the LLR of the corresponding variable node 
    
    % Matrix M: Element (i,j) contains the message from the variable node j to the check node i
    % Matrix E: Element (i,j) contains the message from the check node i to the variable node j
    
for i=1:nCW                 % For each codeword            
                    
    M = ones(n-k,1) * (-2*r(:,i)/sigmaw2)';    % Initialize the LLR with the intrinsic information                
        
    for ii=1:30
         
        if mexEnabled    
            %tic;
            E = cnMess(M,B);                    
            a = a + toc;
        else
            % For each check node compute the messages to the variable nodes        
            for cn=1:(n-k)
                for vn=B{cn}                   
                    inM = M(cn,B{cn}(B{cn}~=vn)); % Messages from all the relevant variables nodes except vn                    
                    E(cn,vn) = prod(sign(inM))*(lntanh(sum(lntanh(abs(inM)))));                
                end            
            end 
        end
                
        tic;
        % For each variable nodes compute the LLR as a sum of intrinsinc and extrinsing information
        for vn=1:n            
            L(vn) = sum(E(A{vn},vn)) - 2*r(vn,i)/sigmaw2;             
        end                                      
        d = d + toc;
        c;    
        % Get the estimated output from the llr
        yCap = zeros(n,1);      
        yCap(L<0) = 1;
        
        
        % If the estimated codeword satisfies
        % all the check rules break the cycle        
        if(sum(mod(H*yCap,2))==0)              
            break;
        else         
            if mexEnabled
                %tic;
                M = vnMess(E,A,(2*r(:,i)/sigmaw2)');
                b = b + toc;
            else            
                for vn=1:n
                    for cn=A{vn}                    
                        M(cn,vn) = sum(E(A{vn}(A{vn}~=cn),vn)) - 2*r(vn,i)/sigmaw2;
                    end
                end 
            end
        end                                     
        
    end

    % Finally extract the payload from the codeword
    % quite easy since the code is in systema form
    u_out((i-1)*k+1:i*k) = yCap(1:k,1);    
end
a
b
d
% Remove the padding bits
u_output = u_out(1:mu);
end


