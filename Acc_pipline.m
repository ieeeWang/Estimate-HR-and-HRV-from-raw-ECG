clc; clear; close all;
% correct the following path if needed.
MATLABdir = 'C:\Users\lwang\surfdrive\MATLAB\imec_case_Lei\';
addpath(genpath(MATLABdir));

%% choose input file
subj1 = '075';
subj2 = '246';
subj3 = '248';
subj4 = '323';

subj = subj1;

%% read files and processing
ECG_file={};
ECG_file{1} = ['Chestpatch_subj_',subj,'_a.xlsx']; 
ECG_file{2} = ['Chestpatch_subj_',subj,'_b.xlsx']; 
ECG_file{3} = ['Chestpatch_subj_',subj,'_c.xlsx']; 
ECG_file{4} = ['Chestpatch_subj_',subj,'_x.xlsx']; 

Acc_set={};
% ECG_set={};
for i = 1:length(ECG_file)
    filename = ECG_file{i};
    filepath = ['files\', filename];
    Acc_set{i} = xlsread(filepath,'accelerometer','A:C');
%     ECG_set{i} = xlsread(filepath,'ecg','A:A');
end

%%
fs=32;
ts = 100; % sec, starting time
squreAcc_set={};
% plot 4 sliding windows for each of a,b,c,x activities
figure
title_list = ['a','b','c','x'];
for i = 1: length(Acc_set)
    Acc=Acc_set{i};
    netAcc = sqrt(sum(Acc.^2,2));
    squreAcc_set{i}=netAcc;
    subplot(2,2,i)
        t=[1:size(Acc,1)]/fs;
        plot(t, Acc); grid on
        hold on
        plot(t, netAcc,'--')
        xlim([ts, ts+2.5])
        xlabel('t (s)')
        legend('x','y','z','net')
        title(title_list(i))
        hold off
end

figure
    for i=1:length(squreAcc_set)
        plot(squreAcc_set{i}); hold on
    end
    grid on
    xlim([fs*ts, fs*ts+2.5*fs])
    title('net Acc in a 2.5s window')
    legend('a','b','c','x')

%% compute mean and std of net Acc
squre_mean = []; squre_std = [];
for i=1:length(squreAcc_set)
    tmp= squreAcc_set{i};
    squre_mean(i) = mean(tmp);
    squre_std(i) = std(tmp);
end

figure
    subplot(2,1,1)
    bar(squre_mean) 
    ylabel('mean')
    grid on
    title(['subj# ',subj,': net Acc'])
    subplot(2,1,2)
    bar(squre_std)
    grid on
    ylabel('sd')
    