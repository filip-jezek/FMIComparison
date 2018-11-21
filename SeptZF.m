% ***********************************************************************************
%              S E P T U M   V O L U M E   Z E R O   F U N C T I O N   for
%         S M I T H   C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
%
%   This function calculates the zero function that solves for the left ventricular
%   volume attributed to the septal wall (bowed out --> postive, bowed in -->
%   negative) knowing the left ventricular volume attributed to the LV free wall,
%   V_lvf and the right ventricular volume attributed to the RV free wall, V_rvf.
%   This solution is implicit so here we are solving it for the initial conditions
%   while this same equation appears in the dXdT and is solved within the ode solver
%   for each time step.
%
%   Model originally created on     21 October 2016
%   Model last modfied on           21 October 2016
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF  	     SeptZF  for   S M I T H   C A R D I O V A S C U L A R   
%                                         S Y S T E M S   M O D E L
% ***********************************************************************************

    function SeptZF_Val = SeptZF(V_spt,V_lv,V_rv,time,CVParam_Struct)
    
    % Unpack all fixed parameters needed
    period = CVParam_Struct.period;
    A = CVParam_Struct.A;
    B = CVParam_Struct.B;
    C = CVParam_Struct.C;
    E_es_spt = CVParam_Struct.E_es_spt;
    V_d_spt = CVParam_Struct.V_d_spt;
    P_0_spt = CVParam_Struct.P_0_spt;
    lambda_spt = CVParam_Struct.lambda_spt;
    V_0_spt = CVParam_Struct.V_0_spt;
    E_es_lvf = CVParam_Struct.E_es_lvf;
    V_d_lvf = CVParam_Struct.V_d_lvf;
    P_0_lvf = CVParam_Struct.P_0_lvf;
    lambda_lvf = CVParam_Struct.lambda_lvf;
    V_0_lvf = CVParam_Struct.V_0_lvf;
    E_es_rvf = CVParam_Struct.E_es_rvf;
    V_d_rvf = CVParam_Struct.V_d_rvf;
    P_0_rvf = CVParam_Struct.P_0_rvf;
    lambda_rvf = CVParam_Struct.lambda_rvf;
    V_0_rvf = CVParam_Struct.V_0_rvf;
    
    % Calculate the driver function at the given time
    tau = time - (floor(time/period) * period);
	e_t = A * exp((-1) * B * (tau-C)^2);
    
    % Calculate the zero function which is derived from the fact that
    %  the pressure attributed to the septal wall is the difference between
    %  the pressure attributed to the left ventricular free wall minus the 
    %  pressure attributed to the left ventricular free wall --> 
    %                    P_spt = P_lvf - P_rvf
    SeptZF_Val = (e_t * E_es_spt * (V_spt-V_d_spt)) + ...
        ((1-e_t) * P_0_spt * (exp(lambda_spt * (V_spt-V_0_spt)) - 1)) - ...
        (e_t * E_es_lvf * (V_lv-V_spt-V_d_lvf)) - ...
        ((1-e_t) * P_0_lvf * (exp(lambda_lvf * (V_lv-V_spt-V_0_lvf)) - 1)) + ...
        (e_t * E_es_rvf * (V_rv+V_spt-V_d_rvf)) + ...
        ((1-e_t) * P_0_rvf * (exp(lambda_rvf * (V_rv+V_spt-V_0_rvf)) - 1));

    end

