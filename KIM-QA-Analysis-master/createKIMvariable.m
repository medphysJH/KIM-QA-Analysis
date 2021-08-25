function KIM = createKIMvariable
% KIM = createKIMvariable

%% Create a KIM variable for testing AnalyseKIMqa

% KIM.KIMTrajFolder = 'E:\KIM\LARK\Westmead QA\Dyn_Large_SI';
KIM.KIMTrajFolder = 'E:\KIM\LARK\Westmead QA\Tx interrupt part 3\Large_SI_AP_Breathhold';

KIM.KIMRobotFile = 'E:\GitHub\KIM-QA-Analysis\KIM-QA-Analysis-master\Robot traces\Stitched traces\LiverTraj_LargeSI70s_robot_20min.txt';
% KIM.KIMRobotFile = 'E:\GitHub\KIM-QA-Analysis\KIM-QA-Analysis-master\Robot traces\Stitched traces\LiverTraj_LargeSIandAPWithBreathHold_robot_20min.txt';
% KIM.KIMRobotFile = 'E:\GitHub\KIM-QA-Analysis\KIM-QA-Analysis-master\Robot traces\Hexa_TypicalLung_Tx.txt';

KIM.KIMcoordFile = 'E:\KIM\LARK\Westmead QA\co-ords.txt';

KIM.KIMparamFile = 'E:\KIM\LARK\Westmead QA\param.txt';

KIM.KIMOutputFolder = 'E:\KIM\LARK\';