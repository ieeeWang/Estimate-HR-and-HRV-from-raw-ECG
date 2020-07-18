% helper functions for analyzing ECG waveforms & extract HR HRV features
% for the case presentation of imec interview
% Lei Wang, ieeewangl@gmail.com

function MyUtilECG = MyUtilECG
    MyUtilECG.plot_ECG_Ann = @plot_ECG_Ann;
    MyUtilECG.ReLU = @ReLU;
    MyUtilECG.preprocessing = @preprocessing;
    MyUtilECG.findQRS = @findQRS;
    MyUtilECG.HRVfindR = @HRVfindR;
    MyUtilECG.remove_artifact_peaks = @remove_artifact_peaks;
    MyUtilECG.metric_peak_dtec = @metric_peak_dtec;
    MyUtilECG.getHRV_featureSet = @getHRV_featureSet;
    MyUtilECG.pading_nan= @pading_nan;
end


function plot_ECG_Ann(ecg, ann, fs)
% INPUT:
%     ecg, ECG waveform data
%     ann, the annotation points (in sec)
%     fs, sampling freq

N=length(ecg);
t=[1:N]/fs; %time period(total sample/Fs )
% ann = round(ann*fs);

figure
    plot(t, ecg); hold on
    plot(t(ann),ecg(ann),'ro'); grid on
    legend('ECG','detected peaks')
end


function X = ReLU(X)
% A ReLU layer performs a threshold operation to each element of the input, 
% where any value less than zero is set to zero.
% This operation is equivalent to ReLU for deeplearning.

loc = find(X<0);
X(loc) = 0;
end


function [ecg_filter, ecg_smooth] = preprocessing(ECG, plt)
%% setting
ecg = ECG;
fs = 256;
N=length(ecg);
t=[1:N]/fs; %time period(total sample/Fs )

%% filter by designfilt, for using zero-phase shift filtering
% .A) IIR highpass 1 Hz 
d1 =designfilt('highpassiir', 'FilterOrder', 6, 'HalfPowerFrequency', 1,...
    'SampleRate', fs);
% .B) IIR notch_50Hz (order=2) % remove 50 Hz power line noise
d2 =designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 49.5,...
    'HalfPowerFrequency2', 50.5, 'SampleRate', fs);

% zero-phase filter
tmp = filtfilt(d1, ecg);
ecg_filter = filtfilt(d2, tmp);

% using average filter to remove glitches to increase the performance of 
% peak detection 
% .C) 5-point moving average
ecg_smooth =smooth(ecg_filter); 


%% save data
% savename='ECG_75a';
% save(savename, 'ecg_smooth')

%% plot raw and filtered signals
if plt == 1
    figure % time domain 
    axt_1 = subplot(311);
    plot(t,ecg); grid on
    ylim([-500 500])
    title('Raw ECG')             
    xlabel('time (s)')
    ylabel('amplitude')
    axt_2 = subplot(312);
    plot(t, ecg_filter,'r'); grid on
    ylim([-500 500])
    title('Filtered ECG')             
    xlabel('time (s)')
    ylabel('amplitude')
    axt_3 = subplot(313);
    plot(t, ecg_smooth,'r'); grid on
    ylim([-500 500])
    title('Smooth ECG')             
    xlabel('time (s)')
    ylabel('amplitude') 
    linkaxes([axt_1 axt_2 axt_3],'xy')

    figure % frequency domain
    axf_1 = subplot(211);
    [P1, f] = plot_spectrum(ecg, fs, 0); 
    plot(f,P1); grid on
    ylim([0 3])
    title('Raw ECG')             
    xlabel('f (Hz)')
    ylabel('|Y(f)|')
    axf_2 = subplot(212);
    [P2, f] = plot_spectrum(ecg_filter, fs, 0); 
    plot(f, P2,'r'); grid on
    title('Filtered ECG')             
    xlabel('f (Hz)')
    ylabel('|Y(f)|')
    linkaxes([axf_1 axf_2],'x')
end
    
end

function [P1, f] = plot_spectrum(X, fs, plt) 
% plot the Amplitude Spectrum of X
% f resolusion =  1/t, t is duration (sec) of X
% DFT here is computed without setting a NFFT, ensuring maximum frequency
% resolution.
% OUTPUT
% P1 - |Y(f)|
% P2 - |Y(f)|^2


Y = fft(X); % same with Y = fft(X, lenth(X)); 
L = length(X);
P_abs = abs(Y/L);
P1 = P_abs(1 : (floor(L/2)+1)); % first half including the median point
P1(2:end-1) = 2*P1(2:end-1);
f = fs/2*(0:(2/L):1); % same with f = Fs*(0:(L/2))/L;
% f resolusion = fs/L = 1/t
P2 = P1.^2; 


if plt==1 % plot
    figure
        plot(f, P1) 
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|Y(f)|')
        grid on
end

end

function locs = findQRS(ecg, peak2detect, plt)
% INPUT:
%     ecg, filtered ECG signals
%     peak2detect, 'R' or 'S' at this moment, more to add later.
%     plt, (=1) plot details
% OUTPUT:
%     locs, times (sec) of target peaks    
%     index (in raw ECG) = round(locs*fs);

%% setting
fs = 256;
N=length(ecg);
t=[1:N]/fs; %time period(total sample/Fs )

%% remove artifacts from ECG before peak detection -- (not good)
% tw =1; % 1 sec to ensure most artifacts removed
% threshold = 300; % based on the data
% [ecg_smooth, ~] = zeros_artifact(ecg, fs, tw, threshold);

ecg_smooth = ecg;

%% Detection of R peaks of ECG
% (1) do MODWT
DWTlevel = 5; % 5-levels DWT is enough, 5th level: [4 8] Hz 
wt = modwt(ecg_smooth,DWTlevel);
wtrec = zeros(size(wt));
% wtrec(4:5,:) = wt(4:5,:); % for slower ECG
wtrec(3:4,:) = wt(3:4,:); % for slower & faster ECG
y_rec = imodwt(wtrec,'sym4');

%% (2) Use the squared absolute values of the signal approximation built from 
% the wavelet coefficients, work on the square magnitudes of the original 
% data, you find the capability of the wavelet transform to isolate the R 
% peaks makes the detection problem much easier. Working on the raw data 
% can cause misidentifications such as when the squared S-wave peak exceeds 
% the R-wave peak.

% option A: hard code to fixed value
% peak2peak = 100; 
% option B: adjustable value according measure session
peak2peak = get_loc_peak2peak(ecg_smooth, fs);
fprintf(['peak2peak amplitude: ', num2str(peak2peak), '\n'])
% used to mormalize re-constructed amplitude near 1
y_nom = 2*y_rec/peak2peak;

switch peak2detect
    case 'S' 
        % (A) use squre of re-constructed y.
        y2 = (abs(y_nom)).^2;
    case 'R'
        % (B) mask (i.e.,=0) the mimus elements of re-constructed y.
        MyUtil = MyUtilECG;% use defined function set: MyUtilECG
        y2 = MyUtil.ReLU(y_nom);
        %y2 = y2.^2;
end
% The following parameters learnt from MIT ECG dataset
[qrspeaks,locs] = findpeaks(y2,t,'MinPeakHeight', 0.35,...
    'MinPeakDistance',0.150);

if plt==1
    figure
    ax1 = subplot(211);
    plot(t, ecg_smooth)
    hold on
    plot(t, y_rec)
    grid on
    legend('ECG smooth','Wavelet Reconstruction')     
    title('Wavelet sym4 reconstruction ECG')             
    xlabel('time (s)')
    ylabel('Amp')

    ax2 = subplot(212);
    plot(t,y2)
    hold on
    plot(locs,qrspeaks,'ro')
    grid on
    linkaxes([ax1 ax2],'x')
    xlabel('time (s)')
    ylabel('Normalized Amp')
    title('R Peaks Localized by Wavelet Transform')
end

end



function Ann_sec = HRVfindR(sig_waveform, plt)
% INPUT:
%     sig_waveform, filtered ECG signals
%     plt, (=1) plot details
% OUTPUT:
%     Ann_sec, times (sec) of target peaks    
%     index (in raw ECG) = round(locs*fs);
% this function requres HRV-master toolbox in MATLAB path

%% settings for automated beat detection
load('qrs_settings.mat')
s = 1;    % s=1 for human ECG settings
Fs = 256; % set your sampling frequency
Beat_min = qrs_settings.Beat_min(s);
Beat_max = qrs_settings.Beat_max(s);
wl_tma = ceil(qrs_settings.wl_tma(s)*Fs);
wl_we  = ceil(qrs_settings.wl_we(s,:).*Fs);
d_fs = Fs;

% Heart beat detection
Ann = [];
seg = ceil(length(sig_waveform)/(300*Fs));
if seg>2        
    for i=0:seg
        sig_waveform_tmp = sig_waveform(max(300*Fs*i-10*Fs,1):min(300*Fs*(i+1),length(sig_waveform)));
        if sum(isnan(sig_waveform_tmp)) ~= length(sig_waveform_tmp)
            Ann_tmp = singleqrs(sig_waveform_tmp,Fs,'downsampling',d_fs,'Beat_min',Beat_min,'Beat_max',Beat_max,'wl_tma',wl_tma,'wl_we',wl_we); 
            Ann = [Ann; Ann_tmp+max(300*Fs*i-10*Fs,1)];
        end
    end
    Ann = Ann(Ann>0 & Ann<=length(sig_waveform));
    Ann = unique(sort(Ann));
    Ann(diff(Ann)<.05*Fs)=[];
else
    Ann = singleqrs(sig_waveform,Fs,'downsampling',d_fs,'Beat_min',Beat_min,'Beat_max',Beat_max,'wl_tma',wl_tma,'wl_we',wl_we);
end 
Ann_sec = Ann/Fs;

% use defined function set: MyUtilECG
if plt==1   
    MyUtilECG1 = MyUtilECG;
    MyUtilECG1.plot_ECG_Ann(sig_waveform, Ann, Fs)
end

end


function locs = remove_artifact_peaks(locs, ecg, fs)
% INPUT:
%     ecg_smooth, filtered ECG signals
%     locs, times (sec) of detected peaks 
% OUTPUT:
%     locs, times (sec) of detected peaks after removing artifacts   

tw = 2; % sliding window in sec
artifacts_threshold = 600;
% index_arti, 0: non-artifacts, 1: artifacts
[~, index_arti] = zeros_artifact(ecg, fs, tw, artifacts_threshold);
index_arti = [index_arti,1]; % in case it reach the end
for i = 1:length(locs)
    index = floor(locs(i)/tw)+1;
    if index_arti(index)==1
        locs(i) = nan;
    end
end
% remove nan
locs(isnan(locs))=[];
end
 

function [ecg, index_arti] = zeros_artifact(ecg, fs, tw, threshold)
% INPUT:
%     ecg, filtered ECG signals
%     tw, sliding window (sec)
% OUTPUT:
%     index_artifacts, 0: non-artifacts, 1: artifacts
%     ecg, ecg with artifact parts setted to 0

N=length(ecg);
ws = fs*tw;
n = floor(N/ws);
index_arti = zeros(1, n);
% sliding window no overlap
for i=1:n
    segm = ecg((i-1)*ws+1:i*ws);
    p2p = max(segm) - min(segm);
    if p2p > threshold
        index_arti(i) = 1;
        ecg((i-1)*ws+1:i*ws) = 0;
    end
end
end
  

function p2p_mean = get_loc_peak2peak(ecg, fs)
N=length(ecg);
ws = fs*2;
n = floor(N/ws);
p2p = [];
% sliding window no overlap
for i=1:n
    segm = ecg((i-1)*ws+1:i*ws);
    p2p(i) = max(segm) - min(segm);
end

p2p_mean = median(p2p);

end



function [recall, precition] = metric_peak_dtec(QRS_sec, locs2)
% evalute the peak detectiton perforamne, comapare with QRS_sec (groud truth)
% both inputs in sec
dis_Threshold = 0.008; % 8 ms, within two sampling points, 2/256 = 0.0078
% dis_Threshold = 0.012; % 3/256 = 0.0117 s
N_ann = length(QRS_sec);
dis = [];
for i = 1: N_ann
    tmp = QRS_sec(i)-locs2;
    dis(i) = min(abs(tmp));
end
posti_dtec = sum(dis<dis_Threshold);
% figure
%     plot(dis); grid on

recall = posti_dtec/N_ann;
precition = posti_dtec/length(locs2);
end


function HRV_featureSet = getHRV_featureSet(RR_filt, fs)
% INPUT:
%     RR_filt, filtered relative RR intervals
%     fs, sampling ratio of ECG
% OUTPUT:
%     HRV_featureSet, commonly used HRV featureSet in both time and freq 
% this function requres HRV-master toolbox in MATLAB path

RR_loc = RR_filt;
RR_loc(isnan(RR_loc))=[];

HR_globa  = HRV.HR(RR_loc,0);
SDNN  = HRV.SDNN(RR_loc,0)*1000;
RMSSD = HRV.RMSSD(RR_loc,0)*1000;
pNN50 = HRV.pNN50(RR_loc,0)*100;
rrHRV = HRV.rrHRV(RR_loc,0);

[pLF,pHF,LFHFratio,VLF,LF,HF] = HRV.fft_val(RR_loc,0,fs);

HRV_featureSet={};
HRV_featureSet.time =[HR_globa,SDNN,RMSSD,pNN50,rrHRV];
HRV_featureSet.freq =[pLF,pHF,LFHFratio,VLF,LF,HF];

end


function C1 = pading_nan(HR_loc, HR_loc_L_max)
L=HR_loc_L_max;
C1 = HR_loc(:);
C1 = [C1; nan(L-length(HR_loc),1)];

end





