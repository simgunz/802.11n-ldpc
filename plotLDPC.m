function [ ] = plotLDPC( x, ber_ldpc, fer_ldpc, varargin)

berColor = {'--sb','--+r','--dg','--om'};
ferColor = {'-sb','-+r','-dg','-om'};
leg = {'LDPC Simulated BER, Rate 1/2','LDPC Simulated FER, Rate 1/2','LDPC Simulated BER, Rate 2/3','LDPC Simulated FER, Rate 2/3','LDPC Simulated BER, Rate 3/4','LDPC Simulated FER, Rate 3/4','LDPC Simulated BER, Rate 5/6','LDPC Simulated FER, Rate 5/6'};

f = figure;

if length(varargin)     % Use variable input length to select the kind of graph
    iteration = x;    
    semilogy(iteration,fer_ldpc,'-rd'); 
    hold on;
    plot(iteration,ber_ldpc,'--gs'); 
    hold off;
    grid on;
    xlabel('Number of iteration');
    ylabel('FER');
    saveas(f,'output/FERvsITER');
    saveas(f,'output/FERvsITER','pdf');
else
    EbN0dB = x;        
    semilogy(EbN0dB(1),1,'r');
    hold on;
    for i=1:length(ber_ldpc(:,1))
        h1(i)=plot(EbN0dB,ber_ldpc(i,:),berColor{i}); 
        h2(i)=plot(EbN0dB,fer_ldpc(i,:),ferColor{i}); 
    end
    hold off; 
    h = [h1; h2];
    h = h(:);    
    leg = leg(1:length(h));
    legend(h,leg);
    xlabel('E_b/N_0 [dB]');
    ylabel('BER/FER');
    grid on;
    saveas(f,'output/BEReFERvsEbN0');
    saveas(f,'output/BEReFERvsEbN0','pdf');
end

end

