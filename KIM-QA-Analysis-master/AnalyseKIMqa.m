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
% Dynamic requires 3 parameters, tx interrupt requires 4
fid = fopen(KIM.KIMparamFile);
paramData = fscanf(fid, '%f');
fclose(fid);

%% Read and extract motion data
% accepts both Robot 6DOF file and Hexamotion 3 DOF files

fid = fopen(KIM.KIMRobotFile);
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

%% Read couchshift file
if exist(fullfile(KIM.KIMTrajFolder, 'couchShifts.txt'),'file') == 2
    fid=fopen(fullfile(KIM.KIMTrajFolder, 'couchShifts.txt'));
    couch.Positions = textscan(fid, '%f,%f,%f\r', 'headerlines', 1);
    fclose(fid);
    
    couch.vrt = couch.Positions{1};
    couch.lng = couch.Positions{2};
    couch.lat = couch.Positions{3};
    couch.lat = couch.lat(couch.lat>950) - 1000;
    
    couch.NumShifts = length(couch.vrt);
    couch.ShiftsAP = -diff(couch.vrt)*10;	% AP maps to couch -vert
    couch.ShiftsSI = diff(couch.lng)*10;    % SI maps to couch long
    couch.ShiftsLR = diff(couch.lat)*10;    % LR maps to couch lat
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
    ShiftIndex = cellfun('size',rawDataKIM,1);
    ShiftIndex = cumsum(ShiftIndex);
    
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
array = [rawDataKIM{1,6:3:3+3*nMar}];
[~, index] = sort(array, 'descend');

for n = 1:nMar
    dataKIM.x_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+2}]';   % LR maps to x
    dataKIM.y_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+3}]';   % SI maps to y
    dataKIM.z_mm(:,n) = [rawDataKIM{:,3+3*(index(n)-1)+1}]';   % AP maps to z
    
    %   C# indexes from 0 to N-1 so a + 1 is added to each 2D trajectory for
    %   equivalent comparison to MATLAB
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
dataKIM.rCent_mm = sqrt(dataKIM.xCent.^2 + dataKIM.yCent.^2 + dataKIM.zCent.^2);

dataKIM.xCentOff = dataKIM.xCent_mm - dataKIM.xCent_mm(1);
dataKIM.yCentOff = dataKIM.yCent_mm - dataKIM.yCent_mm(1);
dataKIM.zCentOff = dataKIM.zCent_mm - dataKIM.zCent_mm(1);
dataKIM.rCentOff = sqrt(dataKIM.xCentOff.^2 + dataKIM.yCentOff.^2 + dataKIM.zCentOff.^2);
disp(KIM)

end