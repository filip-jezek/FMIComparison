

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

      
clear;

%% PARAMETRIZATION
% all params checked manually in the model
% Following params are not inlcuded in the Modelica model:

% Left ventricle free wall parameters
V_d_lvf = 0;                                % LV ES zero P volume (mL)        
V_0_lvf = 0;                                % LV ED pressure param (mL)    
% Right ventricle free wall parameters
V_d_rvf = 0;                                % RV ES zero P volume (mL)
V_0_rvf = 0;                                % RV ED pressure param (mL)
%% SIMULATE

command = [ ...
    'Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow.exe ', ... 
    '-override=startTime=0,stopTime=30,stepSize=0.01,tolerance=1e-06,', ...
    'solver=dassl,outputFormat=mat,variableFilter=.* ', ...
          ]


tic
[s, a] = system(command)
toc

%% LOAD RESULTS
m = load('Cardiovascular.Model.Smith2004.HemodynamicsSmith_shallow_res.mat');
toc
% unpack the modelica mat structure

names = cellstr(m.name')
% names(find(contains(names, 'pressure')))
plot_var_i = find(contains(names, 'Rsys.q_in.pressure'));
% get the position in data_2
data2_i = m.dataInfo(2, plot_var_i);

% get the time
time = m.data_2(1, :);
% get the var
pressure_ao = m.data_2(data2_i, :);
%% VIEW RESULTS
plot(time, pressure_ao);