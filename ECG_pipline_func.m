function output = ECG_pipline_func(filepath, QRSmethod, plt_out)
% this function is based on ECG_pipline.m, to do the same thing in one line
% Here add try-catch to support reading .x file where QRS anntation is absent
ECG = xlsread(filepath,'ecg','A:A');

QRS_Ann_exit = 1; % be default
try
    QRS = xlsread(filepath,'qrs','A:A');
catch
    warning('QRS Annotation does not exit, but it does not matter.');
    QRS_Ann_exit = 0;
end

fs = 256;
plt = 0; % not plot details
%% preprocessing
MyUtil = MyUtilECG;% use defined function set: MyUtilECG
[~, ecg_smooth] = MyUtil.preprocessing(ECG, plt);

%% Detection of R (or S) peaks of ECG
switch QRSmethod
    case 'wavelet'   
        % (option A) my implemented algrithm using MATLAB toolbox functions
        peak2find = 'R'; % 'R' or 'S'
        locs = MyUtil.findQRS(ecg_smooth, peak2find, plt);
    case 'toolbox' 
        % (option B) use algrithms from the open-source toolbox HRV-master
        locs = MyUtil.HRVfindR(ecg_smooth, plt);
end
% mask the detected peaks where siganls are too noisy
locs2 = MyUtil.remove_artifact_peaks(locs, ecg_smooth, fs);

%% check the detected peaks on ECG
N=length(ecg_smooth);
t=[1:N]/fs; %time period(total sample/Fs )
ann = round(locs*fs); % index in raw ECG
ann2 = round(locs2*fs); % index in raw ECG

if QRS_Ann_exit == 1
    QRS_sec = QRS/fs; % sec
    [recall, precition] = MyUtil.metric_peak_dtec(QRS_sec, locs);
    [recall2, precition2] = MyUtil.metric_peak_dtec(QRS_sec, locs2);
end

if plt_out == 1
    figure
    ax_1 = subplot(211);
    plot(t, ecg_smooth); hold on
    plot(t(ann),ecg_smooth(ann),'ro'); grid on
    switch QRS_Ann_exit
        case 1
            plot(t(QRS),ecg_smooth(QRS),'k*')
            legend('ECG smooth','detected peaks', 'R peaks (expert)')
            title(['Without post-processing, recall: ', sprintf('%.2f%%',100*recall),...
            ' precision: ', sprintf('%.2f%%',100*precition)])  
        case 0
            legend('ECG smooth','detected peaks')
             title('Without post-processing') 
    end
    xlabel('time (s)')
%     ylim([-200 200])
    
    ax_2 = subplot(212);
    plot(t, ecg_smooth); hold on
    plot(t(ann2),ecg_smooth(ann2),'ro'); grid on 
    switch QRS_Ann_exit
        case 1   
            plot(t(QRS),ecg_smooth(QRS),'k*')
            legend('ECG smooth','detected peaks', 'R peaks (expert)')
            title(['After post-processing, recall: ', sprintf('%.2f%%',100*recall2),...
            ' precision: ', sprintf('%.2f%%',100*precition2)])  
        case 0
             legend('ECG smooth','detected peaks')
             title('Without post-processing') 
    end     
%     ylim([-200 200])
    xlabel('time (s)') 
    linkaxes([ax_1 ax_2],'xy')
end
%% analyze RR intervals & extract HR HRV features
% get RR intervals and filtering of artifacts
RR = diff(locs2);
RR_filt = HRV.RRfilter(RR,20);

% Computation of HRV measures
HRV_featureSet = MyUtil.getHRV_featureSet(RR_filt, fs);
HRVt = HRV_featureSet.time;
fprintf('global HR is %3.2f bpm \n', HRVt(1))

% compute HR using every n RRs e.g., 10 or 20
nRR = 10; % <<<<<<<<<<<<<<<<<<<<
nRR_filt = RR_filt; % non-successive RRs
nRR_filt(isnan(nRR_filt))=[];
HR_loc = HRV.HR(nRR_filt,nRR);

% compute HR_loc for annotation R peaks
if QRS_Ann_exit == 1
    RR_Ann = diff(QRS_sec); % Annotation
    RR_filt_Ann = HRV.RRfilter(RR_Ann,20);
    RR_filt_Ann(isnan(RR_filt_Ann))=[];
    HR_loc_Ann = HRV.HR(RR_filt_Ann,nRR);
    % Computation of HRV measures
    HRVSet_Ann = MyUtil.getHRV_featureSet(RR_filt_Ann, fs);   
end


if plt_out == 1
    figure
        subplot(211)% visulize RR artifacts filtering
        plot(RR,'--'); hold on
        plot(RR_filt,'--.'); grid on
        xlabel('Index of RR intervals'); 
        ylabel('RR interval (s)')
        legend('RR raw','RR filtered')
        title('RR intervals before and after filter')
        subplot(212)
        plot(HR_loc,'--'); hold on
        if QRS_Ann_exit == 1
            plot(HR_loc_Ann); 
            legend('my detection','from R Annotation')
        end
        grid on
%         ylim([40 120])
        ylabel('HR (bpm)')
        xlabel('Index of valid RR intervals'); 
        title('HR computed from every 10 RRs')
end


%% return
output = {};
output.RRI = nRR_filt;% non-successive RRs with artifacts removed
output.HRVSet = HRV_featureSet;
output.HR_loc = HR_loc;
if QRS_Ann_exit == 1
    output.HR_loc_Ann = HR_loc_Ann;
    output.HRVSet_Ann = HRVSet_Ann;
end



