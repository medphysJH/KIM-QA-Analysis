function AnalyseKIMqa(KIM)
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
listOfTrajFiles = ls([KIM.KIMTrajFolder '\*GA*.txt']);
for n = size(listOfTrajFiles,1):-1:1
   if contains(listOfTrajFiles(n,:), 'ol', 'IgnoreCase', true)
       listOfTrajFiles(n,:) = [];
   end
end
noOfTrajFiles = size(listOfTrajFiles,1);

% Create output file name
prefix = datestr(now, 'yymmdd-HHMM');
[~, RobotFile, ~] = fileparts(KIM.KIMRobotFile);
if length(RobotFile)<20
    middle = RobotFile;
else
    middle = RobotFile(1:20);
end
append = '_Dynamic.txt';
file_output = [prefix '_' middle append];
file_output = fullfile(KIM.KIMOutputFolder, file_output);

% Original code included a latency value of either 0.2 (Dyn) or 0.35 (TxInt)
%   No documentation regarding source or reason for this value was included
%   and there was also no reason for a difference between the values used
latency = 0;

%% Read coordinate file
% First column x (RL), second column y (AP), third column z (IS)
% First 'n' rows are marker co-ordinates
% Last row is the isocentre

fid = fopen(KIM.KIMcoordFile);
coordData = fscanf(fid, '%f %f %f');
fclose(fid);

nMar = (length(coordData)/3-1);

marker_x = sum(coordData(1:3:end-3))/nMar;
marker_y = sum(coordData(2:3:end-3))/nMar;
marker_z = sum(coordData(3:3:end-3))/nMar;

% Marker co-ordinates need to be transformed to machine space (Dicom to IEC)
% In IEC space y is inverted and y and z are switched
Avg_marker_x = 10*(marker_x - coordData(end-2));
Avg_marker_y = 10*(marker_z - coordData(end-1));
Avg_marker_z = -10*(marker_y - coordData(end));

%% Read parameter file
% this may be removed if it isn't needed
% Dynamic requires 3 parameters
if isfile(KIM.KIMparamFile)
    fid = fopen(KIM.KIMparamFile);
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

fid = fopen(KIM.KIMRobotFile);
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
if exist(fullfile(KIM.KIMTrajFolder, 'couchShifts.txt'),'file') == 2
    fid=fopen(fullfile(KIM.KIMTrajFolder, 'couchShifts.txt'));
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
        logfilename = fullfile(KIM.KIMTrajFolder, listOfTrajFiles(traj,:));
        rawDataKIM{traj} = readcell(logfilename, opts);
    end
    ShiftIndex_KIM = cellfun('size',rawDataKIM,1);
    ShiftIndex_KIM = cumsum(ShiftIndex_KIM);
    
    rawDataKIM = vertcat(rawDataKIM{:});
else
    logfilename = fullfile(KIM.KIMTrajFolder, listOfTrajFiles);
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
KIM.NumArcs = sum(abs(diff(diff(dataKIM.Gantry(dataKIM.Gantry<90)+360)>0)))+1;

% Determine the index for treatment start
d= diff(dataKIM.timestamps);
[~, d_index] = sort(d,'descend');
indexOfTreatStart = min(d_index(1:KIM.NumArcs)) + 1;
dataKIM.indexOfTreatStart = indexOfTreatStart;

%% Trajectories for KIM data
% Index the markers by SI position where 1 is the most cranial and 3 the most caudal
% Note: for Varian kV panels x-y origin appears to be top-left of panel
%   (from source persepective); meaning lower y values are more sup

array = [rawDataKIM{1,6:3:3+3*nMar}];
[~, index] = sort(array, 'descend');

for n = 1:nMar
    dataKIM.x_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+2}]';   % LR maps to x
    dataKIM.y_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+3}]';   % SI maps to y
    dataKIM.z_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+1}]';   % AP maps to z
    
    % C# indexes from 0 to N-1 so a + 1 is added to each 2D trajectory for
    %	equivalent comparison to MATLAB
    dataKIM.x_pix(:,n) = [rawDataKIM{:,(3+3*nMar)+2*(index(n)-1)+1}]' + 1;
    dataKIM.y_pix(:,n) = [rawDataKIM{:,(3+3*nMar)+2*(index(n)-1)+2}]' + 1;
end

% Compute centroid for the 2D coordinates
dataKIM.xCent_pix = sum(dataKIM.x_pix,2)/nMar ;
dataKIM.yCent_pix = sum(dataKIM.y_pix,2)/nMar ;

% Compute centroid 3D trajectories for KIM data
dataKIM.r_mm = sqrt(dataKIM.x_mm.^2 + dataKIM.y_mm.^2 + dataKIM.z_mm.^2);

dataKIM.xCent_mm = sum(dataKIM.x_mm,2)/nMar - Avg_marker_x;
dataKIM.yCent_mm = sum(dataKIM.y_mm,2)/nMar - Avg_marker_y;
dataKIM.zCent_mm = sum(dataKIM.z_mm,2)/nMar - Avg_marker_z;
dataKIM.rCent_mm = sqrt(dataKIM.xCent_mm.^2 + dataKIM.yCent_mm.^2 + dataKIM.zCent_mm.^2);

dataKIM.xCentOff = dataKIM.xCent_mm - dataKIM.xCent_mm(1);
dataKIM.yCentOff = dataKIM.yCent_mm - dataKIM.yCent_mm(1);
dataKIM.zCentOff = dataKIM.zCent_mm - dataKIM.zCent_mm(1);
dataKIM.rCentOff = sqrt(dataKIM.xCentOff.^2 + dataKIM.yCentOff.^2 + dataKIM.zCentOff.^2);

%% Align KIM and motion traces
if couch.NumShifts >= 1
    dataKIM.unshifted = [dataKIM.xCent_mm dataKIM.yCent_mm dataKIM.zCent_mm];
    dataMotion.unshifted = [dataMotion.x dataMotion.y dataMotion.z];
    for n = 1:couch.NumShifts
        dataKIM.yCent_mm(ShiftIndex_KIM(n):end) = dataKIM.yCent_mm(ShiftIndex_KIM(n):end) - couch.ShiftsSI(n);
        dataKIM.xCent_mm(ShiftIndex_KIM(n):end) = dataKIM.xCent_mm(ShiftIndex_KIM(n):end) - couch.ShiftsLR(n);
        dataKIM.zCent_mm(ShiftIndex_KIM(n):end) = dataKIM.zCent_mm(ShiftIndex_KIM(n):end) - couch.ShiftsAP(n);
    
        [~,ShiftIndex_Mot(n)] = min(abs(dataMotion.timestamps - dataKIM.timestamps(ShiftIndex_KIM(n))));
        
        dataMotion.x(ShiftIndex_Mot(n):end) = dataMotion.x(ShiftIndex_Mot(n):end) + couch.ShiftsLR(n);
        dataMotion.y(ShiftIndex_Mot(n):end) = dataMotion.y(ShiftIndex_Mot(n):end) + couch.ShiftsSI(n);
        dataMotion.x(ShiftIndex_Mot(n):end) = dataMotion.z(ShiftIndex_Mot(n):end) + couch.ShiftsAP(n);
    end
end

TimeStart = -dataKIM.timestamps(1);
TimeStep = mean(diff(dataMotion.timestamps))/4;
TimeEnd = abs(dataKIM.timestamps(end) - dataMotion.timestamps(end));
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
MotionY = dataMotion.y; %NEEDS TO BE UPDATED TO UNSHIFTED DATAMOTION
KIMy = dataKIM.yCent_mm(1:end_index);
rmseSI = nan(1,length(ShiftValues));
for a = 1:length(ShiftValues)
    interpMotionY = interp1(MotionTime, MotionY, KIM_time + ShiftValues(a));
    rmseSI(a) = sum((KIMy-interpMotionY).^2)/end_index;
end
[~,I] = min(rmseSI);
TimeShift = ShiftValues(I);



%% Results and Output

dataKIM.CorrectedTime = dataKIM.timestamps + TimeShift + latency;

dataKIM.treat.time = dataKIM.CorrectedTime(dataKIM.indexOfTreatStart:end);
dataKIM.treat.x = dataKIM.xCent_mm(dataKIM.indexOfTreatStart:end);
dataKIM.treat.y = dataKIM.yCent_mm(dataKIM.indexOfTreatStart:end);
dataKIM.treat.z = dataKIM.zCent_mm(dataKIM.indexOfTreatStart:end);

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
dataKIM.analysis.AllMotionDiff(:,1) = dataKIM.xCent_mm - dataMotion.interp.x;
dataKIM.analysis.AllMotionDiff(:,2) = dataKIM.yCent_mm - dataMotion.interp.y;
dataKIM.analysis.AllMotionDiff(:,3) = dataKIM.zCent_mm - dataMotion.interp.z;

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
figure, plot(dataKIM.timestamps, dataKIM.xCent_mm, 'bo', dataKIM.timestamps, dataKIM.yCent_mm, 'go', dataKIM.timestamps, dataKIM.zCent_mm, 'ro')
ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);
title('KIM 3DoF motion', 'fontsize', 16);
legend('LR (KIM)', 'SI (KIM)', 'AP (KIM)');
set(gca,'fontsize',16)

if couch.NumShifts > 1
    % Plot KIM SI with couch shifts
    figure, plot(dataKIM.unshifted(:,2), 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 1: KIM with couch shifts', 'fontsize', 16)
    set(gca,'fontsize',16)
    
    % Plot KIM SI corrected for couch shifts
    figure, plot(dataKIM.yCent_mm, 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 2: KIM with couch shifts undone', 'fontsize', 16)
    set(gca,'fontsize',16)
    
    % Plot KIM with expected motion data (unmatched)
    figure, plot(dataMotion.timestamps, dataMotion.y, 'k-', dataKIM.timestamps, dataKIM.yCent_mm, 'g.')
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 3: KIM before time shift', 'fontsize', 16)
    set(gca,'fontsize',16)
    
    % Plot KIM synced with expected motion data
    figure, plot(dataMotion.timestamps, dataMotion.y, 'k-', dataKIM.CorrectedTime, dataKIM.yCent_mm, 'g.')
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('Step 4: KIM with couch shifts undone and with time shift', 'fontsize', 16)
    set(gca,'fontsize',16)
else
    figure, plot(dataKIM.yCent_mm, 'g.', 'linewidth', 3)
    xlabel('Index', 'fontsize',16)
    ylabel('SI position (mm)', 'fontsize',16)
    title('KIM SI trace', 'fontsize', 16)
    set(gca,'fontsize',16)
end

figure
hold on
plot(dataKIM.timestamps, dataKIM.yCent_mm,'gx', dataKIM.timestamps, dataKIM.zCent_mm,'rx', dataKIM.timestamps, dataKIM.xCent_mm,'bx', 'linewidth', 3)
plot(dataMotion.timestamps, dataMotion.y,'g-', dataMotion.timestamps, dataMotion.z,'r-', dataMotion.timestamps, dataMotion.x,'b-')
ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);
title('KIM vs Source motion', 'fontsize', 16);
legend( 'SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)','Location','NorthEastOutside' );
set(gca,'fontsize',16)
hold off

end