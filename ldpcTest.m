clc;
clear all;
close all;

%%% TUNABLE PARAMETERS %%%

mu = 800;        % Input length
iter = 1;          % Number of simulations (10000 for good results)
EbN0step = 0.5;    % Set to 0.5 to speed things up

R=1/2;             % Available rates 1/2, 2/3, 3/4, 5/6

%%%%%% SIMULATION %%%%%%%%

%EbN0 = 1:EbN0step:9;
EbN0 = 1:EbN0step:3.5;
%EbN0 = 2;
EbN0dB = 10*log10(EbN0);

gamma = 2*R*EbN0;
gammaDB = 10*log10(gamma);

Pbit = zeros(1,length(gamma));

tic
for k=1:length(gamma)
    for i=1:iter
        u_input = round(rand(1,mu));       % Random input sequence
        u_output = ldpcTxSystem( u_input, R, gammaDB(k) );
        Pbit(k) = Pbit(k) + sum(u_input ~= u_output);
    end
    Pbit(k) = Pbit(k)/(mu*iter);
end
toc

h = figure;
semilogy(EbN0dB,Pbit,'r');
%hold;
%line([0,10],[1e-5,1e-5],'Color','r');

legend('LDPC Simulated BER');
xlabel('Eb/N0 [dB]');
ylabel('Pbit');
%mkdir('output');
%save('output/workspace');
%saveas(h,'output/figure');
%saveas(h,'output/figure','pdf');