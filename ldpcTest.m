clc;
clear all;
close all;

%%% TUNABLE PARAMETERS %%%

mu = 1000;         % Input length
iter = 3 ;         % Number of simulations (10000 for good results)
EbN0step = 0.5;    % Set to 0.5 to speed things up

R=1/2;             % Available rates 1/2, 2/3, 3/4, 5/6

%%%%%% SIMULATION %%%%%%%%

%EbN0 = 1:EbN0step:9;
EbN0dB = 1:EbN0step:8;
%EbN0dB = fliplr(EbN0dB);
EbN0 = 10.^(EbN0dB/10);
%EbN0dB = 10*log10(EbN0);

gamma = 2*R*EbN0;
gammaUnc = 2*EbN0;
gammaDB = 10*log10(gamma);
gammaUncDB = 10*log10(gammaUnc);

Pbit1 = ones(1,length(gamma));
Pbit2 = ones(1,length(gamma));
Pbit3 = ones(1,length(gamma));

h = figure;

pause(1);
tic
for k=1:length(gamma)
    for i=1:iter
        u_input = round(rand(1,mu));       % Random input sequence
        u_output1 = ldpcTxSystem( u_input, R, gammaDB(k) );
        Pbit1(k) = Pbit1(k) + sum(u_input ~= u_output1);
        u_output2 = ldpcTxSystemWrong( u_input, R, gammaDB(k) );
        Pbit2(k) = Pbit2(k) + sum(u_input ~= u_output2);
        u_output3 = uncodedTxSystem( u_input, gammaUncDB(k) );
        Pbit3(k) = Pbit3(k) + sum(u_input ~= u_output3);
    end
    Pbit1(k) = Pbit1(k)/(mu*iter);
    if(Pbit1(k)==0)
        Pbit1(k)=10^-7;
    end
    Pbit2(k) = Pbit2(k)/(mu*iter);
    if(Pbit2(k)==0)
        Pbit2(k)=10^-7;
    end
    Pbit3(k) = Pbit3(k)/(mu*iter);
    if(Pbit3(k)==0)
        Pbit3(k)=10^-7;
    end    
    
    semilogy(EbN0dB,Pbit1,'--rs');
    hold on;
    plot(EbN0dB,Pbit2,'--gs');
    plot(EbN0dB,Pbit3,'--bs');
    hold off;
    mkdir('output');
    save('output/workspace');
    saveas(h,'output/figure');
    saveas(h,'output/figure','pdf');
    pause(1);
end
time = toc


%hold;
%line([0,10],[1e-5,1e-5],'Color','r');

legend('LDPC Simulated BER','Bad LDPC Simulated BER','Uncoded BER');
xlabel('Eb/N0 [dB]');
ylabel('Pbit');
save('output/workspace');
saveas(h,'output/figure');
saveas(h,'output/figure','pdf');