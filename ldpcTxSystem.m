function [ u_output ] = ldpcTxSystem( u_input, R, gammaDB, mexEnabled)
%CONVOLUTIONALTXSYSTEM Simulates a 802.11 transmission system with ldpc encoding

mexEnabled = 1;     % Enable the C++ version of the decoder
gpuEnabled = 0;
backSub = 1;        % Enable encoding by back substitution

mu = length(u_input);            % Input length

% Codeword length available 648, 1296, 1944 in the standard
% 648, 1296 are used only when there are few input bits
% 1944 will be used here

n = 1944;   % Codeword length
k = n*R;    % Payload length


if backSub
    [H,Z] = buildH(n,R);
else
    [H,G,Z] = buildHG(n,R);
end
%F = H;
%H = sparse(H);

A = cell(1,n);          % Element i contains the indexes of the check nodes in which the variable i is involved
B = cell(1,n-k);        % Element i contains the indexes of the variable nodes involved in the check i

for vn=1:n
    A{vn} = find(H(:,vn))';
end
for cn=1:(n-k)
    B{cn} = find(H(cn,:));
end


sigmaw2 = 1/(10^(gammaDB/10));

%%%%%%% ENCODER %%%%%%%

nCW = ceil(mu/k);                   % Number of codewords
Npad = nCW*k - mu;                  % Number of padding bits
u_in = [ u_input zeros(1,Npad) ];
c = zeros(n,nCW);


if backSub
    H1 = H(:,1:12*Z);
    H2_1 = H(:,12*Z+1:13*Z);    
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
else
    for i=1:nCW
        payload = u_in((i-1)*k+1:i*k)';
        c(:,i) = mod(G*payload,2);
    end
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

%PbitTh = qfunc(sqrt(1/sigmaw2))
%PbitEff = sum(((sign(r)+1)/2)~=c)/n

%%% MESSAGE PASSING DECODER %%%
         

for i=1:nCW                 % For each codeword    
    
    L = zeros(1,n);         % Each element is the LLR of the corresponding variable node 
        
    E = zeros(n-k,n);       % Element (i,j) contains the message from the check node i to 
                            % the variable node j
    M = zeros(n-k,n);       % Element (i,j) contains the message from the variable node i to 
                            % the check node j
%mexEnabled=0;                            
% for xx=1:2
%     if xx == 2    
%         mexEnabled=1;
%     end
%     
    for j=1:n        
        M(:,j) = -2*r(j,i)/sigmaw2;      % Initialize the LLR with the intrinsic information
    end        
    aa = 0;
    bb = 0;
    cc = 0;
    
    %l = fec.ldpcdec(H);
    %a = decode(l,r(:,i)');
    
    for ii=1:30
         tic
        if mexEnabled
            E = cnMess(k,n,M,A,B);
        elseif gpuEnabled
            % For each check node compute the messages to the variable nodes        
            for cn=1:(n-k)
                for vn=B{cn}                   
                    inM = M(cn,B{cn}(B{cn}~=vn)); % Messages from all the relevant variables nodes except vn                    
                    if(isnan(inM))
                        disp('NaN');
                    end
                    z = prod(tanh(inM/2));
                    E(cn,vn) = log((1+z)/(1-z));                
                end            
            end 
        else
            % For each check node compute the messages to the variable nodes        
            for cn=1:(n-k)
                for vn=B{cn}                   
                    inM = M(cn,B{cn}(B{cn}~=vn)); % Messages from all the relevant variables nodes except vn                    
                    if(isnan(inM))
                        disp('NaN');
                    end
                    E(cn,vn) = prod(sign(inM))*(lntanh(sum(lntanh(abs(inM)))));                
                end            
            end 
        end
        aa = aa +toc;
        tic
             % For each variable nodes compute the LLR as a sum of intrinsinc and extrinsing information
        for vn=1:n            
            L(vn) = sum(E(A{vn},vn)) - 2*r(vn,i)/sigmaw2;             
        end                       


        bb = bb + toc;
        tic
        % Bound the LLR from the check nodes to the variable nodes to prevent numerical instability
        %L(L>200) = 200;
        %L(L<-200) = -200;
            
        % Get the estimated output from the llr
        yCap = zeros(n,1);      
        yCap(L<0) = 1;
        
        
        % If the estimated codeword satisfies
        % all the check rules break the cycle        
        if(sum(mod(H*yCap,2))==0)  
            disp('Check ok');
            break;
        else
            for vn=1:n
                for cn=A{vn}
                    idx = A{vn}(A{vn}~=cn);
                    M(cn,vn) = sum(E(idx,vn)) - 2*r(vn,i)/sigmaw2;
                end
            end
        end                                     
        cc = cc+ toc;
    end
    aa/ii;
    bb/ii;
    cc/ii;
    % Finally extract the payload from the codeword
    % quite easy since the code is in systematic form
    u_out((i-1)*k+1:i*k) = yCap(1:k,1)';    
%end
%sum(aa - bb)
end
% Remove the padding bits




u_output = u_out(1:mu);
end


