function [ ] = plotLDPC( x, ber_ldpc, fer_ldpc, varargin)
%PLOTLDPC Summary of this function goes here
%   Detailed explanation goes here

berColor = {'--sb','--+r','--dg','--om'};
ferColor = {'-sb','-+r','-dg','-om'};

f = figure;

if length(varargin)
    iteration = x;    
    semilogy(iteration,fer_ldpc,'-rd'); 
    hold on;
    plot(iteration,ber_ldpc,'--gs'); 
    hold off;
    xlabel('Number of iteration');
    ylabel('FER');
    saveas(f,'output/FERvsITER');
    saveas(f,'output/FERvsITER','pdf');
else
    EbN0dB = x;        
    semilogy(EbN0dB,ones(1,length(ber_ldpc(1,:))));
    hold on;
    for i=1:length(ber_ldpc(:,1))
        plot(EbN0dB,ber_ldpc(i,:),berColor{i}); 
        plot(EbN0dB,fer_ldpc(i,:),ferColor{i}); 
    end
    hold off;    
    legend('LDPC Simulated BER','LDPC Simulated FER');
    xlabel('Eb/N0 [dB]');
    ylabel('BER/FER');
    saveas(f,'output/BEReFERvsEbN0');
    saveas(f,'output/BEReFERvsEbN0','pdf');
end

end

