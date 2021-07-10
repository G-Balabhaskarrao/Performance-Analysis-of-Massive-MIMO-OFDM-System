
function [ofdm chan] = MIMO_OFDM(ofdmIn,chanIn,M)
  
%% Parameters
    % system parameters (independent)
    ofdm.Nb      = 1e3;                 % number of blocks
    ofdm.Nt      = 2;                   % number of transmit antenna    
    ofdm.Nr      = 2;                   % number of receive antenna
    ofdm.K       = M;                 % number of subcarriers    
    ofdm.G       = 1/4;                 % Guard interval percentage    
    ofdm.Mod     = 16;                   % QPSK Modulation    
    ofdm.PSpace  = 1;                   % pilot space between two pilots
    
    
    % channel parameters
    %chan.SNR_dB  = 15;                  % signal to noise ratio
    chan.L       = 6;                   % number of channel taps between each transmit-receive antenna
    
    % control parameters
    ofdm.ifDemodulateData = 1;          % (1,0) if 1, the code demodulates the transmitted data via LS algorithm, and calculates the BER
    ofdm.ifDisplayResults = 1;          % (1,0) if 1, displays the results in the command window
    
    % inserting input data to the default data
    if nargin > 3
        error('Only two arguments can be set as inputs')
    elseif nargin == 3
        % updating the set parameters
        S = fieldnames(ofdmIn);
        for nS = 1:length(S)
             ofdm.(S{nS}) = ofdmIn.(S{nS});
        end
        S = fieldnames(chanIn);
        for nS = 1:length(S)
             chan.(S{nS}) = chanIn.(S{nS});
        end
    elseif nargin == 1
        S = fieldnames(ofdmIn);
        for nS = 1:length(S)
             ofdm.(S{nS}) = ofdmIn.(S{nS});
        end
    end
    % dependent parameters
    ofdm.PPos    = 1:(ofdm.PSpace+1):ofdm.K;    % OFDM pilot positionss
    ofdm.PL      = length(ofdm.PPos);           % Length of pilot subcarriers
    ofdm.DPos    = setxor(1:ofdm.K,ofdm.PPos);  % OFDM data positions
    ofdm.DL      = length(ofdm.DPos);           % Length of data subcarriers
    ofdm.BER     = 0;                           % set the BER to zero
    chan.sigma   = sqrt(10^(-0.1*chan.SNR_dB)); % noise power
    
    % normalization of the energy for the constelation        
        temp         = 0:ofdm.Mod-1;           % possible symbols
        temp         = qammod(temp,ofdm.Mod);  % modulated symbols
        temp         = abs(temp).^2;           % power of each point in the constellation
        temp         = mean(temp);             % average energy of the constellation
        ofdm.ModNorm = 1/sqrt(temp);           % normaliztion factor
    
 
%% Data generation
    % symbol generation
    ofdm.d      = randi(ofdm.Mod,ofdm.DL,ofdm.Nb,ofdm.Nt)-1;   % generation of a DL by nB by Nt matrix of data symbols

%% data Modulation
    ofdm.dMod   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    % memory allocation for the ofdm blocks transmitted from each Tx antenna
    if ofdm.DL > 0
        for nt = 1 : ofdm.Nt
            ofdm.dMod(ofdm.DPos,:,nt) = ofdm.ModNorm*qammod(ofdm.d(:,:,nt),ofdm.Mod);
        end
    end
%% Pilot insertion
    for nt = 1 : ofdm.Nt
        ofdm.dMod(ofdm.PPos,:,nt) = repmat(exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL),1,ofdm.Nb);
    end
    % checking the power of the transmit signal (it has to be 1 after normalization)
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
%% IFFT operation    
    ofdm.ifft   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    % memory allocation for the ofdm blocks transmitted from each Tx antenna after ifft
    for nt = 1 : ofdm.Nt
        ofdm.ifft(:,:,nt) = sqrt(ofdm.K)*ifft(ofdm.dMod(:,:,nt),ofdm.K);
    end
%% Cyclic perfix
    % copy the end of signal to the begining of signal
    ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:,:);ofdm.ifft];
%% Channel
    % for each block we generate a rayleigh fading MIMO channel which is
    % fixed over a block
   chan.Coeff = 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)); 
   
  % chan.Coeff = sqrt(chan.sigma/(chan.sigma+1)) + (sqrt(1/(chan.sigma+1)) )* 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)) ;   

%% Channel pass and filter
    if ofdm.K*ofdm.G < chan.L+1
        error('Guard interval is shorter than channel length, and the system does not function properly')
    end
    ofdm.Y = zeros(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr);
    for nb = 1 : ofdm.Nb
        for nt=1:ofdm.Nt
            for nr=1:ofdm.Nr
                ofdm.Y(:,nb,nr) = ofdm.Y(:,nb,nr) + filter(squeeze(chan.Coeff(nt,nr,:,nb)),1,ofdm.ifftG(:,nb,nt));
            end
        end
    end
    % adding noise
    ofdm.Y = ofdm.Y + chan.sigma*1/sqrt(2)*(         randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)+...
                                            sqrt(-1)*randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)     );
%% Cyclic prefix removal
    ofdm.fftG = ofdm.Y(ofdm.K*ofdm.G+1:ofdm.K*(1+ofdm.G),:,:);
%% FFT operation
    ofdm.fft  = zeros(ofdm.K,ofdm.Nb,ofdm.Nr);
    for nr = 1 : ofdm.Nr
        ofdm.fft(:,:,nr)  = 1/sqrt(ofdm.K)*fft(ofdm.fftG(:,:,nr),ofdm.K);
    end
%% Channel estimation
    % building the first L columns of the fft matrix
    F = dftmtx(ofdm.K);
    F = F(:,1:chan.L);
    % Memory allocation for the estimated channel coefficients
    chan.CoeffEst = zeros(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb);
    for nb = 1 : ofdm.Nb
        for nr = 1 : ofdm.Nr
            % Building matrix A (see the paper)
            chan.A = zeros(ofdm.PL,chan.L*ofdm.Nt);
            for nt = 1 : ofdm.Nt
                chan.A(:,(1:chan.L)+(nt-1)*chan.L) = diag(ofdm.dMod(ofdm.PPos,nb,nt))*F(ofdm.PPos,:);
            end
            ChanEst = pinv(chan.A)*ofdm.fft(ofdm.PPos,nb,nr);
            for nt = 1 : ofdm.Nt
                chan.CoeffEst(nt,nr,:,nb) = ChanEst((1:chan.L)+(nt-1)*chan.L);
            end
        end        
    end
%% MSE (Mean Square error calculation)
    chan.MSE_Simulation = var(chan.Coeff(:)-chan.CoeffEst(:));
   chan.MSE_Theory     = chan.sigma^2/ofdm.PL;
    if ofdm.ifDisplayResults
        disp(['MSE of channel estimation (theory)     is : ',num2str(chan.MSE_Theory)])
       disp(['MSE of channel estimation (simulation) is : ',num2str(chan.MSE_Simulation)])
    end
%% Demodulation
    if ofdm.ifDemodulateData == 1 && ofdm.DL > 0
        % Building channel coefficients in frequency domain
        chan.CoeffEstFreq = zeros(ofdm.K,ofdm.Nt,ofdm.Nr,ofdm.Nb);
        for nb = 1 : ofdm.Nb
            for nr = 1 : ofdm.Nr
                for nt = 1 : ofdm.Nt
                    chan.CoeffEstFreq(:,nt,nr,nb) = F*squeeze(chan.CoeffEst(nt,nr,:,nb));
                end
            end
        end
        
        ofdm.H = chan.CoeffEstFreq(:,:,:,:);
        % demodulation2
        ofdm.dDemod = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nb = 1 : ofdm.Nb
            for dl = 1 : ofdm.DL
                ofdm.dDemod(dl,nb,:) = pinv(reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr)).' ...  % ZF DETECTOR
                                                                                *squeeze(ofdm.fft(ofdm.DPos(dl),nb,:));          
                                         
                %ofdm.dDemod(dl,nb,:) = inv(pagetranspose(ofdm.H) .* ofdm.H + (2*chan.sigma^2*eye(2)))*pagetranspose(ofdm.H)...
                %                                                                               *squeeze(ofdm.fft(ofdm.DPos(dl),nb,:));  
               
                %ofdm.dDemod(dl,nb,:) = inv((reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr).')' * reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr).' + (chan.sigma^2*eye(2))) *(reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr).')'...
                 %                                                                                                                                                                                                                                        *squeeze(ofdm.fft(ofdm.DPos(dl),nb,:));
                                                                                                                                                                                                                                                 
            end
        end
        % detection
        ofdm.dEst = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nt = 1 : ofdm.Nt
            ofdm.dEst(:,:,nt) = qamdemod(1/ofdm.ModNorm * ofdm.dDemod(:,:,nt),ofdm.Mod);
        end
        % BER calculation
        [~,ofdm.BER]  = biterr(ofdm.d(:),ofdm.dEst(:),log2(ofdm.Mod));
        if ofdm.ifDisplayResults
            disp(['BER is = ',num2str(ofdm.BER)])
        end
    end

