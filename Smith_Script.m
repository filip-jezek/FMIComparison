function [T_Sim, P_AOSim, elapsed] = Smith_Script(stopTime)
% ***********************************************************************************
%          S M I T H   C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
%
%   This script uses a patient specific cardiovascular systems model which is
%   based on the Smith el al. model (Med Eng Phys 26:131, 2004) to define 
%   the cardiovascular state of an individual. Currently (30 November 2017) ten 
%   adjustable parameter are manipulated to match eight (right heart cath,
%   RHC, data only) or ten (RHC plus echocardiagraphy data) clinical measures.
%   In this current version the parameters are hand tuned however other versions
%   of this model have been created to methodically and automatically optimize
%   the parameters to have the smallest difference between simulation output and 
%   the clinical measures. 
%  
%   The set of equations in the Smith model is a set of differential - algebraic 
%   equations (DAEs) and the dXdT called contains a single expression where the
%   right hand side is equal to 0 not the dX/dt of a state variable. This means 
%   when ode15s is called it loads a singular mass matrix M with ones along the 
%   diagonal except for a zero in the position of the implicit expression. 
%
%   Model originally created on     17 January 2016
%   Model last modfied on           25 January 2018

%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%   Modified by Filip Ježek, Charles University in Prague, 2018
%
% ***********************************************************************************
%  Start of         P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************

%% **********************************************************************************
%  Load Data for    P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************

  
       
    ID_Num = 'Normal';
    Ave_HR = 80;                            % Average heart rate (beats/min)
    BW = 72.3;                              % Body weight (kg, ~160 lbs)
    Hgt = 178;                              % Height (cm, ~5'10")
    Gender = 'M';                           % Assuming male since TotBV ~ 5L
    
%% **********************************************************************************
%  Fix Params of    P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************

    % Calculate total blood volume based on height, weight and gender. 
    %  This expression is from Nadler et al. Surgery 51:224,1962.
    if (Gender == 'M')
        TotBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000;
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000;
    end
    
    % Elastance function driver parameters
    period = 60/Ave_HR;                             % Period of heart beat (s)
    A = 1;                                          % Elastance function param (uls)
    B = Ave_HR;                                     % Elastance fctn param (1/s^2)
    C = period/2;                                   % Elastance fctn param (s)

    % SMITH ET AL. MODEL PARAMETERS
    % Pericardium calculation parameters
    P_0_pcd = 0.5003;                           % Unstress pericard press (mmHg)
    lambda_pcd = 0.03;                          % Pericard press exp term (1/mL)
    V_0_pcd = 200;                              % Pericard zero P volume (mL)
    P_th = -4;                                  % Thoracic pressure (mmHg)
    % Left ventricle free wall parameters
    E_es_lvf = 2.8798;                          % LV free wall elast (mmHg/mL) 
    V_d_lvf = 0;                                % LV ES zero P volume (mL)        
    P_0_lvf = 0.1203;                           % LV ED pressure param (mmHg)
    lambda_lvf = 0.033;                         % LV ED pressure param (1/mL)
    V_0_lvf = 0;                                % LV ED pressure param (mL)    
    % Right ventricle free wall parameters
    E_es_rvf = 0.585;                           % RV free wall elast (mmHg/mL) 
    V_d_rvf = 0;                                % RV ES zero P volume (mL)
    P_0_rvf = 0.2157;                           % RV ED pressure param (mmHg)
    lambda_rvf = 0.023;                         % RV ED pressure param (1/mL)
    V_0_rvf = 0;                                % RV ED pressure param (mL)
    % Pulmonary artery and vein parameters
    E_es_pa = 0.369;                            % Pulm artery elastance (mmHg/mL)
    V_d_pa = 0;                                 % Pulm artery zero P volume (mL)
    E_es_pu = 0.0073;                           % Pulm vein elastance (mmHg/mL)
    V_d_pu = 0;                                 % Pulm vein zero P volume (mL)
    R_pul = 0.1552;                             % Pulm vasc resist (mmHg*s/mL)
    % Aortic and vena cava parameters
    E_es_ao = 0.6913;                           % Aorta elastance (mmHg/mL)
    V_d_ao = 0;                                 % Aorta zero P volume (mL)
    E_es_vc = 0.0059;                           % Vena cava elastance (mmHg/mL)
    V_d_vc = 0;                                 % Vena cava zero P volume (mL)
    R_sys = 1.0889;                             % Syst art resistance (mmHg*s/mL)
    % Heart valve paramenters
    R_mt = 0.0158;                              % Mitral valve resist (mmHg*s/mL)
    L_mt = 7.6968e-5;                           % Mitrl vlve inert (mmHg*s^2/mL)
    R_av = 0.018;                               % Aortic valve resist (mmHg*s/mL)
    L_av = 1.2189e-4;                           % Aortic vlve inert (mmHg*s^2/mL)
    R_tc = 0.0237;                              % Tricspd vlv resist (mmHg*s/mL)
    L_tc = 8.0093e-5;                           % Tricspd vlv inert (mmHg*s^2/mL)
    R_pv = 0.0055;                              % Pulmon vlv resist (mmHg*s/mL)
    L_pv = 1.4868e-4;                           % Pulmon vlv inert (mmHg*s^2/mL)
    % Septum free wall parameters
    E_es_spt = 48.754;                          % Septum FW elstnce (mmHg/mL)
    V_d_spt = 2;                                % Septum zero P volume (mL)
    P_0_spt = 1.1101;                           % Septum ED pressure param (mmHg)
    lambda_spt = 0.435;                         % Septum ED pressure param (1/mL)
    V_0_spt = 2;                                % Septum ED pressure param (mL)

   
    % Save all fixed parameters in a structure 
    CVParam_Values = {period A B C P_0_pcd lambda_pcd V_0_pcd P_th E_es_lvf ...
        V_d_lvf P_0_lvf lambda_lvf V_0_lvf E_es_rvf V_d_rvf P_0_rvf ...
        lambda_rvf V_0_rvf E_es_pa V_d_pa E_es_pu V_d_pu R_pul E_es_ao ...
        V_d_ao E_es_vc V_d_vc R_sys R_mt L_mt R_av L_av R_tc L_tc R_pv L_pv ...
        E_es_spt V_d_spt P_0_spt lambda_spt V_0_spt};
    CVParam_Fields = {'period' 'A' 'B' 'C' 'P_0_pcd' 'lambda_pcd' 'V_0_pcd' ...
        'P_th' 'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' 'V_0_lvf' ...
        'E_es_rvf' 'V_d_rvf' 'P_0_rvf' 'lambda_rvf' 'V_0_rvf' 'E_es_pa' ...
        'V_d_pa' 'E_es_pu' 'V_d_pu' 'R_pul' 'E_es_ao' 'V_d_ao' 'E_es_vc' ...
        'V_d_vc' 'R_sys' 'R_mt' 'L_mt' 'R_av' 'L_av' 'R_tc' 'L_tc' 'R_pv' ...
        'L_pv' 'E_es_spt' 'V_d_spt' 'P_0_spt' 'lambda_spt' 'V_0_spt'};
    CVParam_Struct = cell2struct(CVParam_Values, ...
            CVParam_Fields,2);
        
     % The original Smith model only circulated a portion of the blood so aortic
    %  pressure dynamics are not lumped into a general arterila systemic 
    %  compartment. Assuming they were simulating a typical 5000 mL total blood 
    %  volume they included only 1500 mL (or 30%) in the circulating volume
    %  therefore we will multiply our calculated TotBV value by 0.3 to yield 
    %  circulating blood volume. The additional factor is used to adjust 
    %  circulating volume in heart failure
    CircBV = 0.30 * TotBV; 
       
    
%% **********************************************************************************
%  Sim Params of    P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************

    TSpan = [0 stopTime];
    % Setting state variable initial conditions 
    %  Note that initial volume division is scaled as a fraction
    %  of circulating blood volume calculated earlier
    V_lv0 = (94.6812/1500) * CircBV;
    V_rv0 = (90.7302/1500) * CircBV;
    V_pa0 = (43.0123/1500) * CircBV;
    V_pu0 = (808.458/1500) * CircBV;
    V_ao0 = (133.338/1500) * CircBV;
    V_vc0 = (329.780/1500) * CircBV;
    Q_mt0 = 245.581;
    Q_av0 = 0;
    Q_tc0 = 190.066;
    Q_pv0 = 0;

    % Calculating the initial condition on the volume due to the septum 
    %  deflection which is an implicit function requiring the use of fsolve
    time = 0;
    SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0,V_rv0,time,CVParam_Struct);
    [V_spt0,SeptZF_Val,ExitFlag,OutData] = fsolve(SeptZF_Hndl,-15);
    % Put into vector to pass to ode15s
    X0(1) = V_lv0;
    X0(2) = V_rv0;
    X0(3) = V_pa0;
    X0(4) = V_pu0;
    X0(5) = V_ao0;
    X0(6) = V_vc0;
    X0(7) = Q_mt0;
    X0(8) = Q_av0;
    X0(9) = Q_tc0;
    X0(10) = Q_pv0;
    X0(11) = V_spt0;
%     % Building the simulation parameter structure
%     SimParam_Values = {TSpan_SS TSpan_Sim X0};
%     SimParam_Fields = {'TSpan_SS' 'TSpan_Sim' 'X0'};
%     SimParam_Struct = cell2struct(SimParam_Values, ...
%             SimParam_Fields,2);
%     
    
%% **********************************************************************************
%  Integrate        P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************
tic
    % Build mass matrix for DAE solver
    M = eye(11);                                    % Put identity on diagonal
    M(11,11) = 0;                                   % Set last expression as a ZFun
    % Set ODE/DAE options
%     ODE_Opts = odeset('Mass',M);
    ODE_Opts = odeset('AbsTol',1e-7,'RelTol',1e-4,'Mass',M); 
    % Solve over the steady state time span with ode15s
%     [T_Sim_SS,X_Sim_SS] = ...
%         ode15s(@dXdT_Smith,TSpan_SS,X0,ODE_Opts,CVParam_Struct);
    % Now solve over the simulation time span
    [T_Sim,X_Sim] = ...
        ode15s(@dXdT_Smith,TSpan, X0,ODE_Opts,CVParam_Struct);
disp('ODE15 simulation:')
    toc    
    
%% **********************************************************************************
%  Plot Fig for     P A T I E N T   S P E C I F I C   S M I T H   C V S   M O D E L
% ***********************************************************************************

    % Run model to get intermediate pressures 
    Num_TSim = size(T_Sim,1);                   % Number of timepoints
    P_LVSim = zeros(Num_TSim,1);                % Preallocate matrices
    P_RVSim = zeros(Num_TSim,1);
    P_AOSim = zeros(Num_TSim,1);
    P_VCSim = zeros(Num_TSim,1);
    P_PASim = zeros(Num_TSim,1);
    P_PUSim = zeros(Num_TSim,1);
    P_RVsyst_DPlot = zeros(Num_TSim,1);
    P_RVdiast_DPlot = zeros(Num_TSim,1);
    P_AOsyst_DPlot = zeros(Num_TSim,1);
    P_AOdiast_DPlot = zeros(Num_TSim,1);
    P_PAsyst_DPlot = zeros(Num_TSim,1);
    P_PAdiast_DPlot = zeros(Num_TSim,1);
    P_PCWave_DPlot = zeros(Num_TSim,1);
    V_LVsystSim = 300;
    V_LVdiastSim = 0;
    % Calculate and store intermediate pressures and volumes
    for i = 1:Num_TSim
        VarOut = dXdT_Smith(T_Sim(i), ...
            X_Sim(i,:),CVParam_Struct,1);
        P_LVSim(i) = VarOut(1);
        P_RVSim(i) = VarOut(2);
        P_AOSim(i) = VarOut(3);
        P_VCSim(i) = VarOut(4);
        P_PASim(i) = VarOut(5);
        P_PUSim(i) = VarOut(6);
        V_LVsystSim = min(V_LVsystSim,X_Sim(i,1));
        V_LVdiastSim = max(V_LVdiastSim,X_Sim(i,1));
    end
    % Calculate simulated cardiac output
    CO_Sim = ((V_LVdiastSim - V_LVsystSim) * Ave_HR) / 1000;
   
    disp('ODE15: getting results')
    toc
    elapsed = toc;
  