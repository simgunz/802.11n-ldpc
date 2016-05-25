clc;clear all;close all;

%% TUNABLE PARAMETERS %%

backSubstitution = 1;       % Enable encoding by back substitution
parallelComputation = 0;    % Do not enable this if the matrices are not been created!

% Parameters presets:
% 0: BER vs Iterations test
% 1: Rate 1/2 accurate BER and FER test
% 2: Rate 5/6 accurate BER and FER test
% 3: Rates comparison, accurate
% 4: Manual
preset = 4;

% Manual mode parameters
mu = 10^5;                % Input length
R = 1/2;                  % Code rate. Available rates 1/2, 2/3, 3/4, 5/6
iter = 1;                 % Number of simulations (100 for good results)
EbN0step = 0.25;          % Set to 0.5 to speed things up
EbN0dB = 1:EbN0step:2.25; % Eb/N0 values, select a number of values multiple of 8
                          % to exploit all the 8 core of the cluster blade


%% PRESET %%

switch preset
    case 0,
        mu = 10^5;
        R = 1/2;
        iter = 1000;
        ldpcIter = 1:50;
        EbN0dB = 1.5;
    case 1,
        mu = 10^7;
        R = 1/2;
        iter = 10;          
        EbN0dB = 0.75:EbN0step:2.5;                       
    case 2,
        mu = 10^7;
        R = 5/6;
        iter = 10;
        EbN0dB = 2.5:EbN0step:4.5;
    case 3,                 
        mu = 10^7;
        R = [1/2, 2/3, 3/4, 5/6];        
        iter = 10;
        EbN0dB = 0.75:EbN0step:4.5;        
end


%% SIMULATION %%

if(parallelComputation && isempty(gcp('nocreate')))
    parpool;   % Enable parallel computation
end

fprintf('Running simulation...');
tic
if ~preset
    all_ber_ldpc = zeros(length(ldpcIter),iter);
    all_fer_ldpc = zeros(length(ldpcIter),iter);
    gammaDB = EbN0dB + 10*log10(2*R);
    parfor i=1:length(ldpcIter)
        for j=1:iter
            u_input = round(rand(1,mu));       % Random input sequence
            [u_output, all_ber_ldpc(i,j), all_fer_ldpc(i,j)] = ldpcTxSystemFast( u_input, R, gammaDB, backSubstitution, ldpcIter(i));
        end    
        disp(ldpcIter(i));      % To get a feedback on the state of the simulation
    end    
    ber_ldpc = sum(all_ber_ldpc,2)/iter;
    fer_ldpc = sum(all_fer_ldpc,2)/iter;            
else
    all_ber_ldpc = zeros(length(R),length(EbN0dB),iter);
    all_fer_ldpc = zeros(length(R),length(EbN0dB),iter);
    for i=1:length(R)
        gammaDB = EbN0dB + 10*log10(2*R(i));
        parfor j=1:length(gammaDB)
            for k=1:iter                
                u_input = round(rand(1,mu));       % Random input sequence
                [u_output, all_ber_ldpc(i,j,k), all_fer_ldpc(i,j,k)] = ldpcTxSystemFast( u_input, R(i), gammaDB(j), backSubstitution);                
            end    
            disp(j);            % To get a feedback on the state of the simulation
        end
    end    
    ber_ldpc = sum(all_ber_ldpc,3)/iter;
    fer_ldpc = sum(all_fer_ldpc,3)/iter;
end
time = toc        % Get the simulation time

if(parallelComputation && ~isempty(gcp('nocreate')))
    delete(gcp('nocreate'));
end

%% SAVE DATA %%
if ~exist('output','dir')
    mkdir('output');
end
save('output/workspace');


%% PLOT %%

if ~preset
    plotLDPC(ldpcIter, ber_ldpc, fer_ldpc, 1);
else
    plotLDPC(EbN0dB, ber_ldpc, fer_ldpc);
end
