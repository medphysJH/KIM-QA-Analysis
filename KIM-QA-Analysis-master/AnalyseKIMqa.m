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

%% Read coordinate file
% First column x (RL), second column y (AP), third column z (IS)
% First 'n' rows are marker co-ordinates
% Last row is the isocentre

fid = fopen(KIM.KIMcoordFile);
coordData = fscanf(fid, '%f %f %f');
fclose(fid);

num_markers = (length(coordData)/3-1);

marker_x = sum(coordData(1:3:end-3))/num_markers;
marker_y = sum(coordData(2:3:end-3))/num_markers;
marker_z = sum(coordData(3:3:end-3))/num_markers;

% Marker co-ordinates need to be transformed to machine space (Dicom to IEC)
% In IEC space y is inverted and y and z are switched
Avg_marker_x = 10*(marker_x - coordData(end-2));
Avg_marker_y = 10*(marker_z - coordData(end-1));
Avg_marker_z = -10*(marker_y - coordData(end));

%% Read parameter file
% this may be removed if it isn't needed
% Dynamic requires 3 parameters, tx interrupt requires 4
fid = fopen(KIM.KIMparamFile);
paramData = fscanf(fid, '%f');
fclose(fid);

%% Read and extract motion data
% accepts both Robot 6DOF file and Hexamotion 3 DOF files

fid=fopen(KIM.KIMRobotFile);
FirstLine = fgetl(fid);
if ~isnumeric(FirstLine) && FirstLine(1)=='t'
    % Hexamotion trajectory files start with 'trajectory'
    isrobot = 0;
    rawMotionData = textscan(fid, '%f %f %f');
else
    % Robot trajectory files *should* start with '0'
    isrobot = 1;
    frewind(fid);
    rawMotionData = textscan(fid, '%f %f %f %f %f %f %f');
    % else
    %     print('Unrecognised motion input file type')
    %     msgbox('Unrecognised motion input file type','Motion File Type');
    %     return
end
fclose(fid);

if isrobot
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

%% Read and extract KIM trajectory data
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableDescriptionsLine = 1; % first line contains varaiable descriptions
opts.DataLines = 2;  % data starts on the second line

if noOfTrajFiles > 1
    for traj = 1:noOfTrajFiles
        logfilename = fullfile(KIM.KIMTrajFolder, listOfTrajFiles(traj));
        %     fopen(logfilename);
        %     rawDataKIM = textscan(fid,    '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s', 'headerLines', 1);
        rawDataKIM{traj} = readmatrix(logfilename, opts);
        %     fclose(fid);
    end
else
    logfilename = fullfile(KIM.KIMTrajFolder, listOfTrajFiles);
    rawDataKIM = readmatrix(logfilename, opts);
end

end