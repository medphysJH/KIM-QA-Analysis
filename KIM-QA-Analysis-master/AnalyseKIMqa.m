function AnalyseKIMqa(KIMmotion)
% AnalyseKIMqa(KIM);
%
% Purpose: Analyse KIM log files acquired during routine QA and compare
%   with the delivered motion files
% Use: Intended for use solely through the UI developed for KIM QA analysis
% Requirements: KIM variable produced by UI
%
% Authors: Jin, Chandrima Sengupta, Jonathan Hindmarsh
% Ver: Aug 2021
% Changes: combine dynamic, treatment interlock, Varian and Elekta versions
%   into one program with common modules and functions

%% setup internal variables

% Find KIM trajectory log files in the specified folder
listOfTrajFiles = ls([KIMmotion.KIMTrajFolder '\*GA*.txt']);
for n = size(listOfTrajFiles,1):-1:1
   if contains(listOfTrajFiles(n,:), 'ol', 'IgnoreCase', true)
       listOfTrajFiles(n,:) = [];
   end
end
noOfTrajFiles = size(listOfTrajFiles,1);

% Create output file name
prefix = datestr(now, 'yymmdd-HHMM');
[~, RobotFile, ~] = fileparts(KIMmotion.KIMRobotFile);
if length(RobotFile)<20
    middle = RobotFile;
else
    middle = RobotFile(1:20);
end
append = '_Dynamic.txt';
file_output = [prefix '_' middle append];
file_output = fullfile(KIMmotion.KIMOutputFolder, file_output);

% Original code included a latency value of either 0.2 (Dyn) or 0.35 (TxInt)
%   No documentation regarding source or reason for this value was included
%   and there was also no reason for a difference between the values used
latency = 0;

%% Read coordinate file
% First column x (RL), second column y (AP), third column z (IS)
% First 'n' rows are marker co-ordinates
% Last row is the isocentre

fid = fopen(KIMmotion.KIMcoordFile);
coordData = fscanf(fid, '%f %f %f');
fclose(fid);

nMar = (length(coordData)/3-1);

marker_x = sum(coordData(1:3:end-3))/nMar;
marker_y = sum(coordData(2:3:end-3))/nMar;
marker_z = sum(coordData(3:3:end-3))/nMar;

% Marker co-ordinates need to be transformed to machine space (Dicom to IEC)
% In IEC space y is inverted and y and z are switched
Avg_marker_x = 10*(marker_x - coordData(end-2));
Avg_marker_y = 10*(marker_z - coordData(end));
Avg_marker_z = -10*(marker_y - coordData(end-1));

%% Read parameter file
% this may be removed if it isn't needed
% Dynamic requires 3 parameters
if isfile(KIMmotion.KIMparamFile)
    fid = fopen(KIMmotion.KIMparamFile);
    paramData = fscanf(fid, '%f');
    fclose(fid);
    if length(paramData) ~= 3
        paramData = [-30 0.01 30];
    end
else
    paramData = [-30 0.01 30];
end

%% Read and extract motion data
% accepts both Robot 6DOF file and Hexamotion 3 DOF files

fid = fopen(KIMmotion.KIMRobotFile);
FirstLine = fgetl(fid);
if ~isnumeric(FirstLine) && FirstLine(1)=='t'
    % Hexamotion trajectory files start with 'trajectory'
    isrobot = 0;
    % Remainder of data is 3 columns of mm values specifying:
    %   LR|IS|PA
    %   where R, S, & A are positive
    rawMotionData = textscan(fid, '%f %f %f');
else
    % Robot trajectory files have no header and *should* start with '0'
    isrobot = 1;
    frewind(fid);
    % Robot data file has 7 columns of data:
    %   Time|x|y|z|rotx|roty|rotz
    %   Time is in seconds, position in mm, rotation in degrees
    %   Directions as per IEC 1217 definition
    rawMotionData = textscan(fid, '%f %f %f %f %f %f %f');
    % else
    %     print('Unrecognised motion input file type')
    %     msgbox('Unrecognised motion input file type','Motion File Type');
    %     return
end
fclose(fid);

if isrobot
    % Robot files are specified in IEC1217 format which is what the rest of
    %   the analysis expects so no adjustment necessary
    dataMotion.x = (1).*rawMotionData{2}(1:end);
    dataMotion.y = rawMotionData{3}(1:end);
    dataMotion.z = (1).*rawMotionData{4}(1:end);
    dataMotion.r = sqrt(dataMotion.x.^2 + dataMotion.y.^2 + dataMotion.z.^2);
    
    dataMotion.xOff = dataMotion.x - dataMotion.x(1);
    dataMotion.yOff = dataMotion.y - dataMotion.y(1);
    dataMotion.zOff = dataMotion.z - dataMotion.z(1);
    dataMotion.rOff = sqrt(dataMotion.xOff.^2 + dataMotion.yOff.^2 + dataMotion.zOff.^2);
    
    dataMotion.timestamps = rawMotionData{1};
else
    % Hexamotion co-ordinates need to be adapted to fit IEC definition
    %   specifically the x direction needs to be inverted so that L is pos+
    dataMotion.x = -1.*rawMotionData{1}(1:end);
    dataMotion.y = rawMotionData{2}(1:end);
    dataMotion.z = rawMotionData{3}(1:end);
    dataMotion.r = sqrt(dataMotion.x.^2 + dataMotion.y.^2 + dataMotion.z.^2);
    
    dataMotion.xOff = dataMotion.x - dataMotion.x(1);
    dataMotion.yOff = dataMotion.y - dataMotion.y(1);
    dataMotion.zOff = dataMotion.z - dataMotion.z(1);
    dataMotion.rOff = sqrt(dataMotion.xOff.^2 + dataMotion.yOff.^2 + dataMotion.zOff.^2);
    
    dataMotion.timestamps = [0:0.02:(length(dataMotion.x)-1)*0.02]';
end

%% Read couchshift file
if exist(fullfile(KIMmotion.KIMTrajFolder, 'couchShifts.txt'),'file') == 2
    fid=fopen(fullfile(KIMmotion.KIMTrajFolder, 'couchShifts.txt'));
    couch.Positions = textscan(fid, '%f,%f,%f\r', 'headerlines', 1);
    fclose(fid);
    
    couch.vrt = couch.Positions{1};
    couch.lng = couch.Positions{2};
    couch.lat = couch.Positions{3};
    couch.lat = couch.lat(couch.lat>950) - 1000;
    
    couch.NumShifts = length(couch.vrt)-1;
    couch.ShiftsAP = -diff(couch.vrt)*10;	% AP maps to couch -vert
    couch.ShiftsSI = diff(couch.lng)*10;    % SI maps to couch long
    couch.ShiftsLR = diff(couch.lat)*10;    % LR maps to couch lat
else
    couch.NumShifts = 0;
end

%% Read and extract KIM trajectory data
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableDescriptionsLine = 1; % first line contains varaiable descriptions
opts.DataLines = 2;  % data starts on the second line

if noOfTrajFiles > 1
    rawDataKIM = cell(noOfTrajFiles,1);
    for traj = 1:noOfTrajFiles
        logfilename = fullfile(KIMmotion.KIMTrajFolder, listOfTrajFiles(traj,:));
        rawDataKIM{traj} = readcell(logfilename, opts);
    end
    ShiftIndex_KIM = cellfun('size',rawDataKIM,1);
    ShiftIndex_KIM = cumsum(ShiftIndex_KIM);
    
    rawDataKIM = vertcat(rawDataKIM{:});
else
    logfilename = fullfile(KIMmotion.KIMTrajFolder, listOfTrajFiles);
    rawDataKIM = readcell(logfilename, opts);
end

dataKIM.timestamps = [rawDataKIM{:,2}]';
dataKIM.timestamps = dataKIM.timestamps - dataKIM.timestamps(1);
dataKIM.Gantry = [rawDataKIM{:,3}]';
dataKIM.index = [rawDataKIM{:,1}]';

% Calculate the number of arcs by looking at the change in gantry rotation
%   Make gantry angles in the file continuous
%   Calculate the change in gantry angle between points
%   Sum the number of times this changes sign (ie rotation direction)
%   Add one to give the number of arcs
KIMmotion.NumArcs = sum(abs(diff(diff(dataKIM.Gantry(dataKIM.Gantry<90)+360)>0)))+1;

% Determine the index for treatment start
d= diff(dataKIM.timestamps);
[~, d_index] = sort(d,'descend');
indexOfTreatStart = min(d_index(1:KIMmotion.NumArcs)) + 1;
dataKIM.indexOfTreatStart = indexOfTreatStart;

%% Trajectories for KIM data
% Index the markers by SI position where 1 is the most cranial and 3 the most caudal
% Note: for Varian kV panels x-y origin appears to be top-left of panel
%   (from source persepective); meaning lower y values are more sup

array = [rawDataKIM{1,6:3:3+3*nMar}];
[~, index] = sort(array, 'descend');

for n = 1:nMar
    dataKIM.x(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+2}]';   % LR maps to x
    dataKIM.y(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+3}]';   % SI maps to y
    dataKIM.z(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+1}]';   % AP maps to z
    
    % C# indexes from 0 to N-1 so a + 1 is added to each 2D trajectory for
    %	equivalent comparison to MATLAB
    dataKIM.x_pix(:,n) = [rawDataKIM{:,(3+3*nMar)+2*(index(n)-1)+1}]' + 1;
    dataKIM.y_pix(:,n) = [rawDataKIM{:,(3+3*nMar)+2*(index(n)-1)+2}]' + 1;
end

% Compute centroid for the 2D coordinates
dataKIM.xCent_pix = sum(dataKIM.x_pix,2)/nMar ;
dataKIM.yCent_pix = sum(dataKIM.y_pix,2)/nMar ;

% Compute centroid 3D trajectories for KIM data
dataKIM.r = sqrt(dataKIM.x.^2 + dataKIM.y.^2 + dataKIM.z.^2);

dataKIM.xCent = sum(dataKIM.x,2)/nMar - Avg_marker_x;
dataKIM.yCent = sum(dataKIM.y,2)/nMar - Avg_marker_y;
dataKIM.zCent = sum(dataKIM.z,2)/nMar - Avg_marker_z;
dataKIM.rCent = sqrt(dataKIM.xCent.^2 + dataKIM.yCent.^2 + dataKIM.zCent.^2);

dataKIM.xCentOff = dataKIM.xCent - dataKIM.xCent(1);
dataKIM.yCentOff = dataKIM.yCent - dataKIM.yCent(1);
dataKIM.zCentOff = dataKIM.zCent - dataKIM.zCent(1);
dataKIM.rCentOff = sqrt(dataKIM.xCentOff.^2 + dataKIM.yCentOff.^2 + dataKIM.zCentOff.^2);

%% Align KIM and motion traces
if couch.NumShifts >= 1
    dataKIM.unshifted = [dataKIM.xCent dataKIM.yCent dataKIM.zCent];
    dataMotion.unshifted = [dataMotion.x dataMotion.y dataMotion.z];
    for n = 1:couch.NumShifts
        dataKIM.yCent(ShiftIndex_KIM(n):end) = dataKIM.yCent(ShiftIndex_KIM(n):end) - couch.ShiftsSI(n);
        dataKIM.xCent(ShiftIndex_KIM(n):end) = dataKIM.xCent(ShiftIndex_KIM(n):end) - couch.ShiftsLR(n);
        dataKIM.zCent(ShiftIndex_KIM(n):end) = dataKIM.zCent(ShiftIndex_KIM(n):end) - couch.ShiftsAP(n);
    
        [~,ShiftIndex_Mot(n)] = min(abs(dataMotion.timestamps - dataKIM.timestamps(ShiftIndex_KIM(n))));
        
        dataMotion.x(ShiftIndex_Mot(n):end) = dataMotion.x(ShiftIndex_Mot(n):end) + couch.ShiftsLR(n);
        dataMotion.y(ShiftIndex_Mot(n):end) = dataMotion.y(ShiftIndex_Mot(n):end) + couch.ShiftsSI(n);
        dataMotion.x(ShiftIndex_Mot(n):end) = dataMotion.z(ShiftIndex_Mot(n):end) + couch.ShiftsAP(n);
    end
else
    dataKIM.unshifted = [dataKIM.xCent dataKIM.yCent dataKIM.zCent];
    dataMotion.unshifted = [dataMotion.x dataMotion.y dataMotion.z];
end

% TimeStart = -dataKIM.timestamps(1);
TimeStart = -dataKIM.timestamps(round(length(dataKIM.timestamps)/2));
TimeStep = mean(diff(dataMotion.timestamps))/4;
% TimeEnd = abs(dataKIM.timestamps(end) - dataMotion.timestamps(end));
TimeEnd = abs(dataKIM.timestamps(round(length(dataKIM.timestamps)/2)) - dataMotion.timestamps(end));
ShiftValues = TimeStart:TimeStep:TimeEnd;
% ShiftValues = paramData(1):paramData(2):paramData(3);

if dataKIM.timestamps(end) > dataMotion.timestamps(end)
    end_index = find(dataKIM.timestamps < dataMotion.timestamps(end), 1, 'last');
    KIM_time = dataKIM.timestamps(1:end_index);
else
    end_index = length(dataKIM.timestamps);
    KIM_time = dataKIM.timestamps;
end
MotionTime = dataMotion.timestamps;
Motion = smoothdata(sqrt(sum(dataMotion.unshifted.^2,2)),'gaussian',1/mean(diff(dataMotion.timestamps))); % Smooth the source motion using a 1 second window
rKIMmotion = sqrt(dataKIM.xCent.^2 + dataKIM.yCent.^2 + dataKIM.zCent.^2);
start_index = find(rKIMmotion<1,1); % Find the first datapoint < 1 to trim the start of the KIM signal incase there are artefacts when tracking starts
KIMmagnitude = smoothdata(rKIMmotion(start_index:end_index),'gaussian',5);  % smooth the KIM data using a 5 data point window
KIM_time = KIM_time(start_index:end);

Motion = Motion - min(Motion);
KIMmagnitude = KIMmagnitude - min(KIMmagnitude); 

rmseSI = nan(1,length(ShiftValues));
% figure(1)
% pause(0.1)
for a = 1:length(ShiftValues)
    mod_time = KIM_time + ShiftValues(a);
    interpMotion = interp1(MotionTime, Motion, mod_time(mod_time>0 & mod_time<MotionTime(end)));
    rmseSI(a) = sum((KIMmagnitude(mod_time>0 & mod_time<MotionTime(end))-interpMotion).^2)/end_index;
    
%     if rem(a,100)==0
%         plot(mod_time(mod_time>0 & mod_time<MotionTime(end)), interpMotion, 'k-', mod_time(mod_time>0 & mod_time<MotionTime(end)), KIMmotion(mod_time>0 & mod_time<MotionTime(end)), 'g.')
%         pause(0.05)
%     end
end
% close(1)

[~,I] = min(rmseSI);
TimeShift = ShiftValues(I);

%% Results and Output

dataKIM.CorrectedTime = dataKIM.timestamps + TimeShift + latency;

dataKIM.treat.time = dataKIM.CorrectedTime(dataKIM.indexOfTreatStart:end);
dataKIM.treat.x = dataKIM.xCent(dataKIM.indexOfTreatStart:end);
dataKIM.treat.y = dataKIM.yCent(dataKIM.indexOfTreatStart:end);
dataKIM.treat.z = dataKIM.zCent(dataKIM.indexOfTreatStart:end);

% Interpolate source motion to match KIM timepoints and ensure there are no
%   NaN values
% Timepoints corresponding to entire collected KIM data (pre-arc + treatment)
dataMotion.interp.x = fillmissing(interp1(dataMotion.timestamps, dataMotion.x, dataKIM.CorrectedTime),'nearest');
dataMotion.interp.y = fillmissing(interp1(dataMotion.timestamps, dataMotion.y, dataKIM.CorrectedTime),'nearest');
dataMotion.interp.z = fillmissing(interp1(dataMotion.timestamps, dataMotion.z, dataKIM.CorrectedTime),'nearest');
% Timepoints corresponding to KIM data during treatment
dataMotion.treat.x = fillmissing(interp1(dataMotion.timestamps, dataMotion.x, dataKIM.treat.time),'nearest');
dataMotion.treat.y = fillmissing(interp1(dataMotion.timestamps, dataMotion.y, dataKIM.treat.time),'nearest');
dataMotion.treat.z = fillmissing(interp1(dataMotion.timestamps, dataMotion.z, dataKIM.treat.time),'nearest');

% *For treatment only*
% Calculate the difference between KIM detected position and expected
%   position as specified by the motion source
% For ease of data processing rearrange positional data into column: x(LR), y(SI), z(AP)
dataKIM.analysis.TxMotionDiff(:,1) = dataKIM.treat.x - dataMotion.treat.x;
dataKIM.analysis.TxMotionDiff(:,2) = dataKIM.treat.y - dataMotion.treat.y;
dataKIM.analysis.TxMotionDiff(:,3) = dataKIM.treat.z - dataMotion.treat.z;

% Calculate mean, stdev and percentile for the positinal data
dataKIM.analysis.TxResults{1,:} = mean(dataKIM.analysis.TxMotionDiff,1);
dataKIM.analysis.TxResults{2,:} = std(dataKIM.analysis.TxMotionDiff,0,1);
dataKIM.analysis.TxResults{3,:} = tsprctile(dataKIM.analysis.TxMotionDiff,[5 95],1);

% *For all acquired KIM data (including prearc data)*
% Repeat above
% For ease of data processing rearrange positional data into column: x(LR), y(SI), z(AP)
dataKIM.analysis.AllMotionDiff(:,1) = dataKIM.xCent - dataMotion.interp.x;
dataKIM.analysis.AllMotionDiff(:,2) = dataKIM.yCent - dataMotion.interp.y;
dataKIM.analysis.AllMotionDiff(:,3) = dataKIM.zCent - dataMotion.interp.z;

% Calculate mean, stdev and percentile for the positinal data
dataKIM.analysis.AllResults{1,:} = mean(dataKIM.analysis.AllMotionDiff,1);
dataKIM.analysis.AllResults{2,:} = std(dataKIM.analysis.AllMotionDiff,0,1);
dataKIM.analysis.AllResults{3,:} = tsprctile(dataKIM.analysis.AllMotionDiff,[5 95],1);

failname = {' LR,', ' SI,', ' AP,'};
if any([dataKIM.analysis.AllResults{1,:}]>1)
    OutputText{1,1} = 'QA result: KIM FAILED in Dynamic test';
    OutputText{2,1} = ['Tested trajectory ', RobotFile, ': mean difference of',]; 
    OutputText{2,1} = [OutputText{2,1} failname{[dataKIM.analysis.AllResults{1,:}]>1} ' > or = 2 mm' newline];
elseif any([dataKIM.analysis.AllResults{2,:}]>2)
    OutputText{1,1} = 'QA result: KIM FAILED in Dynamic test';
    OutputText{2,1} = ['Tested trajectory ', RobotFile, ': standard deviation of difference of',]; 
    OutputText{2,1} = [OutputText{2,1} failname([dataKIM.analysis.AllResults{1,:}]>2) ' > or = 2 mm' newline];
else
    OutputText{1,1} = 'QA result: KIM PASSED in Dynamic test';
    OutputText{2,1} = ['Tested trajectory ' RobotFile newline]; 
end

OutputText{3,1} = sprintf('No. of couch shifts: %u\n', couch.NumShifts);
OutputText{4,1} = sprintf('Processing time per image (Online): %.3f\n', mean(diff(dataKIM.timestamps)));
OutputText{5,1} = sprintf('Mean\t\t\tStd\t\t\tPercentile(5,95)');
OutputText{6,1} = sprintf('LR\tSI\tAP\tLR\tSI\tAP\tLR\t\tSI\t\tAP');
OutputText{7,1} = sprintf('%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)', ...
                    [dataKIM.analysis.TxResults{1,:}], [dataKIM.analysis.TxResults{2,:}], [dataKIM.analysis.TxResults{3,:}]);
OutputText{8,1} = newline;
OutputText{9,1} = sprintf('Time shift to match motion traces is = %.4g seconds', TimeShift);

writecell(OutputText,file_output, 'Delimiter','space', 'QuoteStrings', false)

%% Plots
% Plot basic KIM x-y-z data
figure, plot(dataKIM.timestamps, dataKIM.xCent, 'bx', dataKIM.timestamps, dataKIM.yCent, 'gx', dataKIM.timestamps, dataKIM.zCent, 'rx')
ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);
title('KIM 3DoF motion', 'fontsize', 16);
legend('LR (KIM)', 'SI (KIM)', 'AP (KIM)', 'Location', 'best');
set(gca,'fontsize',16)

if couch.NumShifts > 1
    % Plot KIM SI with couch shifts
    figure, plot(dataKIM.unshifted(:,2), 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 1: KIM with couch shifts', 'fontsize', 16)
    set(gca,'fontsize',16)
    
    % Plot KIM SI corrected for couch shifts
    figure, plot(dataKIM.yCent, 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 2: KIM with couch shifts undone', 'fontsize', 16)
    set(gca,'fontsize',16)
    
    % Plot KIM with expected motion data (unmatched)
    figure, plot(dataMotion.timestamps, dataMotion.y, 'k-', dataKIM.timestamps, dataKIM.yCent, 'g.')
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 3: KIM before time shift', 'fontsize', 16)
    set(gca,'fontsize',16)
else
    figure, plot(dataKIM.yCent, 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('KIM SI trace', 'fontsize', 16)
    set(gca,'fontsize',16)
end

% Plot KIM synced with expected motion data
figure
hold on
plot(MotionTime, Motion, 'kx', KIM_time + TimeShift, KIMmagnitude, 'g.')
plot(MotionTime, sqrt(sum(dataMotion.unshifted.^2,2)), 'c-', KIM_time + TimeShift, rKIMmotion(start_index:end_index), 'r--')
xlabel('Index', 'fontsize',16)
ylabel('Position (mm)', 'fontsize',16)
title('Smoothed source motion vs KIM detected motion', 'fontsize', 16)
legend('Smoothed Source', 'Smoothed KIM', 'Source', 'KIM', 'Location', 'best');
set(gca,'fontsize',16)
hold off
ImageFilename = [prefix '_' middle '_KIM_Source_synced.jpg'];
% print('-djpeg','-r300',fullfile(KIMmotion.KIMOutputFolder,ImageFilename))

figure
hold on
plot(dataKIM.CorrectedTime, dataKIM.yCent,'gx', dataKIM.CorrectedTime, dataKIM.zCent,'rx', dataKIM.CorrectedTime, dataKIM.xCent,'bx', 'linewidth', 3)
plot(dataMotion.timestamps, dataMotion.y,'g-', dataMotion.timestamps, dataMotion.z,'r-', dataMotion.timestamps, dataMotion.x,'b-')
ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);
title('KIM vs Source motion', 'fontsize', 16);
legend('SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)','Location', 'best' );
set(gca,'fontsize',16)
hold off

ImageFilename = [prefix '_' middle '_3Dmotion.jpg'];
% print('-djpeg','-r300',fullfile(KIMmotion.KIMOutputFolder,ImageFilename))

end