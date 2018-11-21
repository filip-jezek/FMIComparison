function [time, outs, t_elapsed] = RunOMSim(filename, outputs, stopTime, stepSize)
% Runs a simulation of exe simulator from openmodelica

% command = [ ...
%     'C:/Users/User/AppData/Local/Temp/OpenModelica/OMEdit/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow.exe ', ... 
%     '-port=55507 -logFormat=xmltcp ', ...
%     '-override=startTime=0,stopTime=300,stepSize=0.06,tolerance=1e-06,', ...
%     'solver=dassl,outputFormat=mat,variableFilter=.* ', ...
%     '-r=C:/Users/User/AppData/Local/Temp/OpenModelica/OMEdit/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow_res.mat ', ...
%     '-w -lv=LOG_STATS ', ...
%     '-inputPath=C:/Users/User/AppData/Local/Temp/OpenModelica/OMEdit/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow ', ...
%     '-outputPath=C:/Users/User/AppData/Local/Temp/OpenModelica/OMEdit/Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow ' ...
%           ]
% 
    
%% SIMULATE
command = [ ...
    filename, ' ', ... 
    '-override startTime=0,stopTime=', num2str(stopTime), ',stepSize=', num2str(stepSize), ',tolerance=1e-07,', ...
    'solver=dassl,outputFormat=mat,variableFilter=.* ', ...
          ];


tic
[s, a] = system(command);
disp('OpenModelica exe simulation time:')
toc

%% LOAD RESULTS
m = load('FMITest.Smith_patSpec_res.mat');

% unpack the modelica mat structure
names = cellstr(m.name');
% names(find(contains(names, 'pressure')))
plot_var_i = find(contains(names, outputs));
% get the position in data_2
data2_i = m.dataInfo(2, plot_var_i);

% get the time
time = m.data_2(1, :);
% get the var
outs = m.data_2(data2_i, :)';

disp('OpenMOdelica loading results:')
toc
t_elapsed = toc;