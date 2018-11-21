%% init
clear
clf
clc
current_folder = fileparts(mfilename('fullpath'));


% timesn(1, :) = cell2mat(times)

t = 50;

rep = 10;
% timeseries = [1*ones(rep, 1); 5*ones(rep, 1); 10*ones(rep, 1)];
timeseries = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89];
for i = 1:length(timeseries)
    for j = 1:rep
        k = i*j;
        t = timeseries(i);
        s = 0.01;
        clear names times;
        names = {};
        times = {};
        clf
        clc
        disp(['COMMENCING RUN n.', num2str(k), ', simulating for ', num2str(t), ' seconds for ', num2str(j), 'x time'])

        %% DYMOLA ME
        fmu_path=fullfile(current_folder,'SmithDymolaME.fmu');
        tic
        fmu=FMUModelME2(fmu_path);
        disp('Loading FMU - Dymola ME:');
        toc
        output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        [time, yout, yname] = simulate(fmu, [0,t], 'Output', output);
        disp('simulating FMU Dymola ME:')
        toc
        times{end+1} = toc;
        names{end+1} = "DYMOLA ME";

        figure(1);subplot(211);hold off;plot(time, yout(:, 1)/133.32);
        figure(1);subplot(212);hold off;plot(time, yout(:, 2)*1000*1000/60);
        titlies(names);

        %% DYMOLA CS
        fmu_path=fullfile(current_folder,'SmithDymolaCS.fmu');
        tic
        fmu=FMUModelCS2(fmu_path);
        disp('FMU - Dymola CS loading:');
        toc
        output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        [time, yout, yname] = simulate(fmu, 0:s:t, 'Output', output);
        disp('Dymola CS:')
        toc
        times{end+1} = toc;
        names{end+1} = "DYMOLA CS";

        figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);
        titlies(names);

        %% DYMOLA CS CVODE
        fmu_path=fullfile(current_folder,'SmithDymolaCvodeCS.fmu');
        tic
        fmu=FMUModelCS2(fmu_path);
        disp('FMU - Dymola Cvode CS:');
        toc
        output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        [time, yout, yname] = simulate(fmu, 0:s:t, 'Output', output);
        disp('FMU Dymola Cvode CS:')
        toc
        times{end+1} = toc;
        names{end+1} = "DYMOLA CS CVODE";

        figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);
        titlies(names);

        %% SIMULINK FMIToolbox CS - output
        tic
        simout = sim('SmithCS', 'StopTime', num2str(t));
        disp('Simulink FMIToolbox CS')
        toc
        times{end+1} = toc;
        names{end+1} = "Simulink FMIToolbox CS";

        figure(1);subplot(211);hold on;plot(simout.simout.Time, simout.simout.Data/133.32);
        figure(1);subplot(212);hold on;plot(simout.simout.Time, zeros(size(simout.simout.Time)));
        titlies(names);
        %% OpenModelica simulator

        [time, yout, elapsed] = RunOMSim('FMITest.Smith_patSpec.exe', {'aorta.q_in.pressure','aorticValve.q_in.q'}, t, s);
        times{end+1} = elapsed;
        names{end+1} = "OpenModelica EXE";

        figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);
        titlies(names);

        %% REFERENCE MATLAB IMLPEMENTATION

        [time, yout, elapsed] = Smith_Script(t);
        times{end+1} = elapsed;
        names{end+1} = "Matlab ODE15";

        figure(1);subplot(211);hold on;plot(time, yout(:, 1));
        figure(1);subplot(212);hold on;plot(time, zeros(size(time)));
        titlies(names);
        %% commented out - crashes or not working
        %% SIMULINK FMIToolbox ME - output is rubbish and error occures
        % An error occurred while running the simulation and the simulation was terminated
        % Caused by:
        % Solver encountered an error while simulating model 'SmithME' at time 9.1569641639858244 and cannot continue. Please check the model for errors.
        % Nonlinear iteration is not converging with step size reduced to hmin (0.0000000000000000E+000) at time 9.1569641639857900E+000.  Try reducing the minimum step size and/or relax the relative error tolerance.

        %% SIMULINK 2018a FMI block - not working
        % % Value type mismatch for parameter 'aorta' in 'SmithCS/FMU'.
        % % Caused by:
        % % Value type does not match the structure of variable 'aorta' defined in the modelDescription.xml file.
        % % Component:Simulink | Category:Block error
        % 
        %% OPENMODELICA ME1 - broken output
        % 
        % fmu_path=fullfile(current_folder,'SmithOMME1.fmu');
        % tic
        % fmu=FMUModelME1(fmu_path);
        % disp('Loading FMU:');
        % toc
        % output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        % [time, yout, yname] = simulate(fmu, 0:s:t, 'Output', output);
        % disp('simulating FMU Dymola Cvode CS:')
        % toc
        % figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        % figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);
        % 
        %% OPENMODELICA ME2 - not working
        % 
        % fmu_path=fullfile(current_folder,'SmithOMME2.fmu');
        % tic
        % fmu=FMUModelME2(fmu_path);
        % disp('Loading FMU:');
        % toc
        % output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        % [time, yout, yname] = simulate(fmu, [0 t], 'Output', output);
        % disp('simulating FMU Dymola Cvode CS:')
        % toc
        % figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        % figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);
        % 
        % 
        %% OPENMODELICA CS2 - crashes the matlab (2018a)
        % 
        % fmu_path=fullfile(current_folder,'SmithOMCS2.fmu');
        % tic
        % fmu=FMUModelME2(fmu_path);
        % disp('Loading FMU:');
        % toc
        % output.name={'aorta.q_in.pressure','aorticValve.q_in.q'};
        % [time, yout, yname] = simulate(fmu, 0:s:t, 'Output', output);
        % disp('simulating FMU Dymola Cvode CS:')
        % toc
        % figure(1);subplot(211);hold on;plot(time, yout(:, 1)/133.32);
        % figure(1);subplot(212);hold on;plot(time, yout(:, 2)*1000*1000/60);

        %% output legends

        %% populate timeseries
        timesn(i, j, :) = cell2mat(times);
    end
end
%% Plot the time results
colors = get(gca,'ColorOrder');
figure(2);clf;hold on;
markers = {'.', 'o', 'x', '+', 'v', '*', 's', 'd', '^', '<', '>', 'p', 'h'};
handles = [];
for method = 1:length(names)
    for j = 1:size(timesn, 2)
        plot(timeseries, squeeze(timesn(:, j, method)), markers{method}, 'Color', colors(method, :));  
    end
    handles(method) = plot(timeseries, mean(squeeze(timesn(:, :, method)), 2), 'LineWidth', 2, 'Color', colors(method, :));
end

title('Simulation times');
xlabel('Simulated time [s]');
ylabel('Elapsed time [s]')
legend(handles, names);
