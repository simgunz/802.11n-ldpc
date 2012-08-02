function [ u_output ] = uncodedTxSystem( u, gammaDB )
%CONVOLUTIONALTXSYSTEM Simulates a full transmission system with uncoded input

%%% BPAM MODULATOR %%%%

sTx = u;
sTx(sTx==0) = -1;    % Map 0 to -1 and produce the transmitted signal


%%%%%%% CHANNEL %%%%%%%

r = awgn(sTx,gammaDB);      % Received signal


%%%%% DEMODULTOR %%%%%%

u_output = sign(r);             % Decision
u_output(u_output==-1) = 0;     % Invese mapping

end

