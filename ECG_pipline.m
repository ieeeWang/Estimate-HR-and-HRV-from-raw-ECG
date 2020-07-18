clc; clear; close all;
% change the following path in your PC.
MATLABdir = 'C:\Users\lwang\surfdrive\MATLAB\imec_caseECG_Lei_submit\';
addpath(genpath(MATLABdir));

%% read files
subjName1 = 'Chestpatch_subj_075';
subjName2 = 'Chestpatch_subj_246';
subjName3 = 'Chestpatch_subj_248';
subjName4 = 'Chestpatch_subj_323';

subj = subjName1; 

filename1 = [subj,'_a.xlsx']; 
filename2 = [subj,'_b.xlsx']; 
filename3 = [subj,'_c.xlsx']; 
filename4 = [subj,'_x.xlsx']; 

filepath = ['files\', filename1];
ECG = xlsread(filepath,'ecg','A:A');
QRS = xlsread(filepath,'qrs','A:A');
% Acc = xlsread(filepath,'accelerometer','A:C');
fs = 256;

%% preprocessing
MyUtil = MyUtilECG;% use defined function set: MyUtilECG
plt = 1; % plot details
[ecg_filter, ecg_smooth] = MyUtil.preprocessing(ECG, plt);

%% Detection of R (or S) peaks of ECG
% (option A) my implemented algrithm using MATLAB toolbox functions
peak2find = 'R'; % 'R' or 'S'
plt = 1; % plot details
locs = MyUtil.findQRS(ecg_smooth, peak2find, plt);

% (option B) use algrithms from the open-source toolbox HRV-master
% locs = MyUtil.HRVfindR(ecg_smooth, plt);

% mask the detected peaks where siganls are too noisy
locs2 = MyUtil.remove_artifact_peaks(locs, ecg_smooth, fs);

%% check the detected peaks on ECG
N=length(ecg_smooth);
t=[1:N]/fs; %time period(total sample/Fs )
ann = round(locs*fs); % index in raw ECG
ann2 = round(locs2*fs); % index in raw ECG

QRS_sec = QRS/fs; % sec
[recall, precition] = MyUtil.metric_peak_dtec(QRS_sec, locs);
[recall2, precition2] = MyUtil.metric_peak_dtec(QRS_sec, locs2);

figure
ax_1 = subplot(211);
    plot(t, ecg_smooth); hold on
    plot(t(ann),ecg_smooth(ann),'ro'); grid on
    plot(t(QRS),ecg_smooth(QRS),'k*')
    xlabel('time (s)')
%     ylim([-200 200])
    legend('ECG smooth','detected peaks', 'R peaks (expert)')
    title(['Without post-processing, recall: ', sprintf('%.2f%%',100*recall),...
        ' precision: ', sprintf('%.2f%%',100*precition)])    
ax_2 = subplot(212);
    plot(t, ecg_smooth); hold on
    plot(t(ann2),ecg_smooth(ann2),'ro'); grid on
    plot(t(QRS),ecg_smooth(QRS),'k*')
%     ylim([-200 200])
    xlabel('time (s)')
    legend('ECG smooth','detected peaks', 'R peaks (expert)')
    title(['After post-processing, recall: ', sprintf('%.2f%%',100*recall2),...
        ' precision: ', sprintf('%.2f%%',100*precition2)])   
linkaxes([ax_1 ax_2],'xy')

%% analyze RR intervals & extract HR and HRV features
% get RR intervals and filtering of artifacts
% RR = diff(QRS_sec); % Annotation
RR = diff(locs2); % my detection
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
RR_Ann = diff(QRS_sec); % Annotation
RR_filt_Ann = HRV.RRfilter(RR_Ann,20);
RR_filt_Ann(isnan(RR_filt_Ann))=[];
HR_loc_Ann = HRV.HR(RR_filt_Ann,nRR);
    
    
figure
    subplot(211)% visulize RR artifacts filtering
    plot(RR,'--'); hold on
    plot(RR_filt,'--.'); grid on
    xlabel('Index of RR intervals'); 
    ylabel('RR interval (s)')
    legend('RR raw','RR filtered')
    title('RR intervals before and after filter')
    subplot(212)
    plot(HR_loc,'--');hold on
    plot(HR_loc_Ann); grid on
    legend('my detection','from Annotation')
    ylabel('HR (bpm)')
    xlabel('Index of valid RR intervals'); 
    title('HR computed from every 10 RRs')

    

