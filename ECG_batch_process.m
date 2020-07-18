clc; clear; close all;
% change the following path in your PC.
MATLABdir = 'C:\Users\lwang\surfdrive\MATLAB\imec_caseECG_Lei_submit\';
addpath(genpath(MATLABdir));

%% choose one subject to analyze
subj1 = '075';
subj2 = '246';
subj3 = '248';
subj4 = '323';

subj = subj4;
%% read files and processing
ECG_file={};
ECG_file{1} = ['Chestpatch_subj_',subj,'_a.xlsx']; 
ECG_file{2} = ['Chestpatch_subj_',subj,'_b.xlsx']; 
ECG_file{3} = ['Chestpatch_subj_',subj,'_c.xlsx']; 
ECG_file{4} = ['Chestpatch_subj_',subj,'_x.xlsx']; 

outputset={};
for ii = 1:length(ECG_file)
    filename = ECG_file{ii}
    filepath = ['files\', filename];
    plt = 1;
    QRSmethod = 'wavelet'; %'wavelet' or 'toolbox'
%     QRSmethod = 'toolbox'; 
    outputset{ii} = ECG_pipline_func(filepath, QRSmethod, plt);
end

%% save outputs for statistics
% savefolder = 'result\';
% savepath = [MATLABdir, savefolder];
% savefineName = ['outputset_',subj];
% save([savepath,savefineName], 'outputset')

%% plot my estimated HRs in a, b, c, x activites together
figure
    HR_global = [];
    HR_loc_mean = [];
    for i = 1:length(outputset)
        tmp = outputset{i};
        HR_global(i)= tmp.HRVSet.time(1);
        plot(tmp.HR_loc([10:end]),'--'); % starting at 10th RRs
        hold on
        HR_loc_mean(i)= mean(tmp.HR_loc);
    end
%     ylim([40, 140])
    ylabel('HR (bpm)')
    xlabel('Index of valid RR intervals (starting at 10th)'); 
    title(['Subj#',subj,': HR computed from every (previous) 10 valid RRs'])
    grid on
    legend(['a: gHR=',sprintf('%.1f',HR_global(1))],...
        ['b: gHR=',sprintf('%.1f',HR_global(2))],...
        ['c: gHR=',sprintf('%.1f',HR_global(3))],...
        ['x: gHR=',sprintf('%.1f',HR_global(4))]);

%% plot my estimated HRs (a, b, c, x) with annotation results
figure
    HR_global = [];
    for i = 1:length(outputset)
        tmp = outputset{i};
        HR_global(i)= tmp.HRVSet.time(1);
        plot(tmp.HR_loc([10:end]),'--'); % starting at 10th RRs
        hold on
    end
    HR_global_Ann = [];
    for i = 1:3
        tmp = outputset{i};
        HR_global_Ann(i)= tmp.HRVSet_Ann.time(1);
        plot(tmp.HR_loc_Ann([10:end])); % starting at 10th RRs
        hold on
    end
    grid on
%     ylim([40, 140])
    ylabel('HR (bpm)')
    xlabel('Index of valid RR intervals (starting at 10th)'); 
    title(['Subj#',subj,': HR computed from every (previous) 10 valid RRs'])    
    legend(['a (my): gHR=',sprintf('%.1f',HR_global(1))],...
           ['b (my): gHR=',sprintf('%.1f',HR_global(2))],...
           ['c (my): gHR=',sprintf('%.1f',HR_global(3))],...
           ['x (my): gHR=',sprintf('%.1f',HR_global(4))],...
           ['A (expert): gHR=',sprintf('%.1f',HR_global_Ann(1))],...
           ['B (expert): gHR=',sprintf('%.1f',HR_global_Ann(2))],...
           ['C (expert): gHR=',sprintf('%.1f',HR_global_Ann(3))]);










