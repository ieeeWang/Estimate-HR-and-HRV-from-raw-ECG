clear all; close all;clc

%% Example 1: Basic grouped box plot with legend 
% @ ...\+iosr\+statistics\boxPlot.m  --> to change
% scatterMarker = 'x' % Marker used for scatter plots of underlying data.

  y = randn(50,4,3); % samples*groups*items(legend variable)
  x = [1 2 3 4];
  y(1:25) = NaN;

  figure
  h = iosr.statistics.boxPlot(x,y,...
      'symbolColor','k',...
      'medianColor','k',...
      'symbolMarker',{'+','o','d'},...
      'boxcolor',{[1 0 0]; [0 1 0]; [0 0 1]},...
      'groupLabels',{'y1','y2','y3'},...
      'showLegend',true);
  box on

  
figure
iosr.statistics.boxPlot(x,y,...
    'symbolColor','k',...
    'medianColor','k',...
    'symbolMarker',{'o','o','o'},...
    'boxcolor','auto',...
    'xSeparator',true,...
     'groupLabels',{'y1','y2','y3'},'showLegend',true,...
      'showScatter',true);
box on


% Grouped notched box plot with x separators and hierarchical labels
  figure;
  iosr.statistics.boxPlot(x,y,...
      'notch',true,...
      'medianColor','k',...
      'symbolMarker',{'+','o','d'},...
      'boxcolor','auto',...
      'style','hierarchy',...
      'xSeparator',true,...
      'groupLabels',{{'G1','G2','G3'}});
  box on
  
% return
%% Example 2: Box plot with legend labels from data

      % load data
      % (requires Statistics or Machine Learning Toolbox)
      load carbig

      % arrange data
      [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);

      % sort
      IX = [1 3 2]; % order
      g = g{1}(IX);
      y = y(:,:,IX);

      % plot
      figure
      h = iosr.statistics.boxPlot(x,y,...
          'symbolColor','k','medianColor','k','symbolMarker','+',...
          'boxcolor',{[1 1 1],[.75 .75 .75],[.5 .5 .5]},...
          'scalewidth',true,'xseparator',true,...
          'groupLabels',g,'showLegend',true);
      box on
      title('MPG by number of cylinders and period')
      xlabel('Number of cylinders')
      ylabel('MPG')

%    Box plot calculated from weighted quantiles

      % load data
      load carbig

      % random weights
      weights = rand(size(MPG));

      % arrange data
      [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
      weights_boxed = iosr.statistics.tab2box(Cylinders,weights,when);

      % plot
      figure
      h = iosr.statistics.boxPlot(x,y,'weights',weights_boxed);

%    Draw a violin plot
      y = randn(50,3,3);
      x = [1 2 3.5];
      y(1:25) = NaN;
      figure('color','w');
      h2 = iosr.statistics.boxPlot(x,y, 'showViolin', true, 'boxWidth', 0.025, 'showOutliers', false);
      box on
  