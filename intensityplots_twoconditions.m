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
pre = [20 18 0 21 16 17 6 0 3];
post = [32 19 25 17 19 0 26 29 32];

% total number of timepoints needed to represent the entire data set
% (after pre- and post-padding data for individual embryos).
numTimePoints = 152;

% index within the padded data arrays that corresponds to the alignment timepoint
alignPoint = 34;

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
cyto = [14 16 34 13 18 17 28 34 31]; 

% % get average value during maintenance for each embryo
% baseline = zeros(1,eNum);
% for i = 1:eNum
%     baseline(i) = mean(intensities{i}.rawint_den(1:cyto(i),:));
% end
% 
% %normalize
% for i = 1:eNum
%     intensities{i}.normalized_int = (intensities{i}.rawint_den)/baseline(i);
% end

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

%plot control embryos
for i = [1 2 9] 
   cutInt = intensities{i}.znorm_int(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color',[0, 0.753, 1]); 
end

%plot RNAi embryos
for i = [4 5 6 7]
   cutInt = intensities{i}.znorm_int(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color',' #f39c12'); 
end

xlabel('time (s)');
ylabel('normalized intensity');
ylim([-0.1 4.6]);

set(gca,'fontsize', 14) 
saveas(f,[pwd,'/nmy2RNAi_intensity_zscorenorm.svg']);

%% plot raw data
r = figure;
hold on;

%plot control embryos
for i = [1 2 9] 
   cutInt = intensities{i}.rawint_den(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color',[0, 0.753, 1]); 
end

%plot RNAi embryos
for i = [4 5 6 7]
   cutInt = intensities{i}.rawint_den(alignPoint:end);
   plot(adjTimeVec, cutInt, 'color','red'); 
end

xlabel('time (s)');
%ylabel('normalized intensity');
%xlim([0 1000]);
ylim([0.1 4.6]);

%saveas(f,[outputFolder '/nmy2RNAi_intensity.svg']);
