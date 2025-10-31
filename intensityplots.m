%% setup
clear all;
outputFolder = 'plots';
mkdir(outputFolder);
workingDir = pwd; %set wd to current

timeInterval = 10; % time interval between frames (seconds)

% load intensity files
CSVs = dir(fullfile(workingDir, '*.csv'));
eNum = length(CSVs);
intensities = cell(1, eNum);
for k = 1:length(CSVs)
    fileName = fullfile(workingDir, CSVs(k).name);
    intensities{k} = readtable(fileName);
end

% predetermined amounts to pad arrays by so as to align all data to a
% specific timepoint
pre = [0 11 24 16 27 13 6];
post = [8 0 7 15 3 9 5];

% total number of timepoints needed to represent the entire data set
% (after pre- and post-padding data for individual embryos).
numTimePoints = 134;

% index within the padded data arrays that corresponds to the alignment timepoint
alignPoint = 33;

% % number of frames to plot befoe the alignment timepoint
% preFrames = 0;
% 
% % number of frames to plot after the alignment timepoint
% postFrames = 134;

% clean up files
newHeaders = {'time', 'mean','int_den', 'rawint_den'};
for k = 1:eNum
    intensities{k}.Properties.VariableNames = newHeaders;
    intensities{k}.time = ((intensities{k}.time) - 1) * 10; %convert from frames to time (s)
end

%% normalize densities to maintenance phase, pad to align timepoints
% frame in which cytokinesis begins for each embryo (furrow visible)
cyto = [33 22 9 17 6 20 27]; %[e1 e2 e3 e4 e6]

% get average value during maintenance for each embryo
baseline = zeros(1,eNum);
for i = 1:eNum
    baseline(i) = mean(intensities{i}.rawint_den(1:cyto(i),:));
end

%normalize
for i = 1:eNum
    intensities{i}.normalized_int = (intensities{i}.rawint_den)/baseline(i);
end

%z-score normalization
for i = 1:eNum
    znorm = normalize(intensities{i}.rawint_den, "zscore");
    shifted_znorm = znorm - min(znorm);  % Shift so minimum is 0
    intensities{i}.znorm_int = shifted_znorm;
end

%padding to align to the beginning of cytokinesis
for i = 1:eNum
    n = pre(i);
    if n > 0
        prepad = array2table(nan(n, width(intensities{i})),'VariableNames', intensities{1}.Properties.VariableNames);
        intensities{i} = [prepad;intensities{i}];
    else
        %do nothing
    end
end

for i = 1:eNum
    x = post(i);
    if x > 0
        postpad = array2table(nan(x, width(intensities{i})),'VariableNames', intensities{1}.Properties.VariableNames);
        intensities{i} = [intensities{i};postpad];
    else
        %do nothing   
    end
end

%% plot 
f = figure;
hold on;

timeVec = 1:timeInterval:timeInterval*numTimePoints;
adjTimeVec = timeVec(alignPoint:end) - timeVec(alignPoint);

for i = 1:eNum
   cutInt = intensities{i}.znorm_int(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color',[0, 0.753, 1]); 
end

% formatting 
xlabel('time (s)');
ylabel('normalized intensity');
xlim([0 1000]);
ylim([-0.1 4]);
set(gca,'fontsize', 14) 
%darkmode(f); --- this was supposed to make a dark background but it's not
%compatible w this version of matlab
saveas(f,[pwd,'/wt_intensity_zscorenorm.svg']);
%% 
figure;
hold on;
for i = [3 4 ]
   cutInt = intensities{i}.int_den(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color',[0, 0.753, 1]);
end
