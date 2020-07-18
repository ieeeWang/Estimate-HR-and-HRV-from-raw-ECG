clc; clear; close all;
% correct the following path if needed.
MATLABdir = 'C:\Users\lwang\surfdrive\MATLAB\imec_case_Lei\';
addpath(genpath(MATLABdir));

%% read files
subj={};
subj{1} = 'Chestpatch_subj_075_x.xlsx';
subj{2} = 'Chestpatch_subj_246_x.xlsx';
subj{3} = 'Chestpatch_subj_248_x.xlsx';
subj{4} = 'Chestpatch_subj_323_x.xlsx';

ECG_X = [];
for i = 1:4
    filename = subj{i}
    filepath = ['files\', filename];
    ECG_X(i,:) = xlsread(filepath,'ecg','A:A');
end

%%
figure
for i = 1:4
    subplot(4,1,i)
    plot(ECG_X(i,:))
end
