 %% Clear 
    clear
    close all 
    clc
%% Sweep parameters
    SNR_dBV     = 0:1:10;  % vector of SNR values in dB
    SNR=10.^(SNR_dBV/10);
    SNR_dBVL    = length(SNR_dBV);   % length of SNR vector
    nMonteCarlo = 5;               % number of Monte Carlo to find the value of each point in the figure
    ofdmIn.Nt   = 2;                 % number of transmit antennas
    ofdmIn.Nr   = 4;                 % number of recieve antennas
    ofdmIn.ifDisplayResults    = 0;  % turn off the display
    % other parameters of ofdm can also be set. see help of MIMO_OFDM_LSE_CHAN_EST
%% Outputs
    MSE_CHAN_SIM = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in simulation
    MSE_CHAN_THR = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in theory
    MSE_CHAN_BER = zeros(nMonteCarlo,SNR_dBVL);     % BER of the MIMO OFDM with LSE channel estimation
%% Loop
cnt3 = 0;
M=[64 256];
for jj=1:length(M)
for cnt1 = 1 : SNR_dBVL
    chanIn.SNR_dB = SNR_dBV(cnt1); % load the SNR value    
    for cnt2 = 1 : nMonteCarlo
        [ofdm chan] = MIMO_OFDM(ofdmIn,chanIn,M(jj));
        MSE_CHAN_SIM(cnt2,cnt1) = chan.MSE_Simulation;
        MSE_CHAN_THR(cnt2,cnt1) = chan.MSE_Theory;
        MSE_CHAN_BER(cnt2,cnt1,jj) = ofdm.BER;
        
        % update the loop counter
        cnt3 = cnt3 + 1;
        disp([num2str(round(cnt3/(SNR_dBVL*nMonteCarlo)*1000)/10),' is done...'])        
    end   
end
end    


SNR_dBV     = 0:1:10;  % vector of SNR values in dB
    SNR=10.^(SNR_dBV/10);
    SNR_dBVL    = length(SNR_dBV);   % length of SNR vector
    nMonteCarlo = 5;               % number of Monte Carlo to find the value of each point in the figure
    ofdmIn.Nt   = 2;                 % number of transmit antennas
    ofdmIn.Nr   = 4;                 % number of recieve antennas
    ofdmIn.ifDisplayResults    = 0;  % turn off the display
    % other parameters of ofdm can also be set. see help of MIMO_OFDM_LSE_CHAN_EST
%% Outputs
    MSE_CHAN_SIM1 = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in simulation
    MSE_CHAN_THR1 = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in theory
    MSE_CHAN_BER1 = zeros(nMonteCarlo,SNR_dBVL);     % BER of the MIMO OFDM with LSE channel estimation
%% Loop
cnt3 = 0;
M=[64 256];
for jj=1:length(M)
for cnt1 = 1 : SNR_dBVL
    chanIn.SNR_dB = SNR_dBV(cnt1); % load the SNR value    
    for cnt2 = 1 : nMonteCarlo
        [ofdm chan] = MIMO_OFDM1(ofdmIn,chanIn,M(jj));
        MSE_CHAN_SIM1(cnt2,cnt1) = chan.MSE_Simulation;
        MSE_CHAN_THR1(cnt2,cnt1) = chan.MSE_Theory;
        MSE_CHAN_BER1(cnt2,cnt1,jj) = ofdm.BER;
        
        % update the loop counter
        cnt3 = cnt3 + 1;
        disp([num2str(round(cnt3/(SNR_dBVL*nMonteCarlo)*1000)/10),' is done...'])        
    end   
end
end    







    semilogy(SNR_dBV,mean(MSE_CHAN_BER(:,:,1)),'rx-')
    hold on
    semilogy(SNR_dBV,mean(MSE_CHAN_BER(:,:,2)),'cd-')
    semilogy(SNR_dBV,mean(MSE_CHAN_BER1(:,:,1)),'b')
    semilogy(SNR_dBV,mean(MSE_CHAN_BER1(:,:,2)),'gp-')
    xlabel('SNR [dB]')
    ylabel('Bit error rate (BER)')
    grid on
    %legend('16 QAM','64 QAM','128 QAM')
    legend('64 subcarrier ','256 subcarrier')
    %legend('64 subcarrier ')
    title('ZF VS MMSE')
    