clc; clear; close all;
% change the following path in your PC.
MATLABdir = 'C:\Users\lwang\surfdrive\MATLAB\imec_caseECG_Lei_submit\';
addpath(genpath(MATLABdir));

%% load results
savefolder = 'result\';
resultpath = [MATLABdir, savefolder];

% search under path
contents = dir([resultpath '\*.mat']); % get all the files
N=size(contents,1);
result_set = {};
for j=1:N
    message_str = sprintf('Processing: %i of %i .mat data', j, N);
    disp(message_str);
    filename=contents(j).name;
    load ([resultpath, filename]); 
    result_set{j}=outputset;
end

%%
MyUtil = MyUtilECG;% for using defined function set: MyUtilECG
% get the max length of HR_loc for all subj
HR_loc_L=[];
k=0;
for j = 1:4
    subj1 = result_set{j};
    for i = 1:4
        k=k+1;
        acti = subj1{i};
        HR_loc_L(k) = length(acti.HR_loc);
    end
end
HR_loc_L_max = max(HR_loc_L);

% put all HR_loc into a 3D matrix, pading with nan for shorter ones
HR_loc_3D = [];
for j = 1:4 % 4 subj
    subj1 = result_set{j};
    for i = 1:4 % 4 acti
        acti = subj1{i};
        HR_loc_pading = MyUtil.pading_nan(acti.HR_loc, HR_loc_L_max);
        HR_loc_3D(j,i,:)=HR_loc_pading;
    end
end
% adjust dim to: samples*groups*items
Y = permute(HR_loc_3D, [3,1,2]);
X = [1 2 3 4];

%% use a nice toolbox to visulize 4*4 activities in one figure
figure
iosr.statistics.boxPlot(X,Y,...
    'symbolColor','k',...
    'medianColor','k',...
    'symbolMarker',{'o','o','o','o'},...
    'boxcolor','auto',...
     'groupLabels',{'a','b','c','x'},... % LABELS
     'showLegend',true,...
      'showScatter',true);
set(gca,'XTickLabel',{'Subj#075','Subj#246','Subj#248','Subj#323'}, 'FontSize', 10)
ylabel('10-RRs based HRs (bpm)') 
title('HRs (every ~10s) of 4 subj * 4 activities') 
box on
grid on