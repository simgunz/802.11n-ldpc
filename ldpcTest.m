clc;
clear all;
close all;

%% TUNABLE PARAMETERS %%

mu = 10^6;         % Input length
iter = 1 ;         % Number of simulations (10000 for good results)
EbN0step = 0.25;    % Set to 0.5 to speed things up

mexEnabled = 1;             % Enable the C++ version of the decoder
backSubstitution = 1;       % Enable encoding by back substitution

R=1/2;             % Available rates 1/2, 2/3, 3/4, 5/6

%% PRESET %%
EbN0dB = 1:EbN0step:3;

preset = 0;
% 0: Iteration test
% 1: Rate 1/2 accurate test
% 2: Rate 1/2 longer input, less accurate test
% 3: Rate comparison less accurate test
% 4: Manual
switch(preset)
    case 0:
        mu = 10^6;
        R = 1/2;
        iter = 20;
        ldpcIter = 1:50;
        EbN0 = 2.5;        
    case 1:
        mu = 10^6;
        R = 1/2;
        iter = 100;
        break;
    case 2:
        mu = 10^7;
        R = 1/2;
        iter = 10;
        break;
    case 3:
        mu = 10^6;
        R = [1/2, 2/3, 3/4, 5/6];
        iter = 24;
        break;
end

%% SIMULATION %%


gammaDB = EbN0dB + 10*log(2*R);

ber_ldpc = zeros(iter,length(gamma),length(R));
fer_ldpc = zeros(iter,length(gamma),length(R));

tic
if preset == 0
    for j=1:length(ldpcIter)
        parfor i=1:iter
            u_input = round(rand(1,mu));       % Random input sequence
            [u_output, ber_ldpc(i,j), fer_ldpc(i,j)] = ldpcTxSystem( u_input, R, gammaDB, mexEnabled, backSubstitution, ldpcIter(j));
            %u_output2 = ldpcTxSystemWrong( u_input, R, gammaDB(k) );
            %Pbit2(k) = Pbit2(k) + sum(u_input ~= u_output2);
        end    
    end
else
    for k=1:length(R)
        for j=1:length(gamma)
            parfor i=1:iter
                u_input = round(rand(1,mu));       % Random input sequence
                [u_output, ber_ldpc(i,j,k), fer_ldpc(i,j,k)] = ldpcTxSystem( u_input, R, gammaDB(k), mexEnabled, backSubstitution);
                %u_output2 = ldpcTxSystemWrong( u_input, R, gammaDB(k) );
                %Pbit2(k) = Pbit2(k) + sum(u_input ~= u_output2);
            end    
        end
    end
end
        
time = toc

ber_ldpc = sum(ber_ldpc)/iter;
fer_ldpc = sum(fer_ldpc)/iter;

%% PLOT %%

%% SAVE DATA %%
save('output/workspace');
