%----------------------------------------------
% AER 715 Avionics and Systems
% Lab 3 – “Flight Control - Estimation of Model Parameters and Simulation”
% Bryan Kikuta - 501039179
% 
% 
%----------------------------------------------
%
%% Introduction
% Type your introduction in this section
%
%% Post Lab Exercises –
% Put your exercises in this section
%
%% Conclusion
% Write your lab conclusion for the WHOLE lab in this section.
%% Question 1
inch_to_m = 0.0254;
M_h = 1.422;  % [kg]
M_c = 1.916;  % [kg]
L_a = 25.75 * inch_to_m;  % [m]
L_b = 18.125 * inch_to_m;  % [m]

% Compute moment of inertia of elevation axis
J_e = M_h * L_a^2 + M_c * L_b^2;
J_e = J_e * 1.05;  % Add 5% for the lever arm
fprintf("The polar moment of intertia for the elevation axis is: %d [kg*m^2]\n", J_e);

%% Question 2
% Assume the net moment is zero at steady level flight
% Net moment equation: 0 = W_h*L_a - F_t*L_a - W_c*L_b
% Solving for force yields: F_t = (W_h*L_a - W_c*L_b) / L_a
g = 9.81;  % [m/s^2]
W_h = M_h * g;  % [N]
W_c = M_c * g;  % [N]
F_t = (W_h*L_a - W_c*L_b) / L_a;
fprintf("The force required for steady level flight is: %d [N]\n", F_t);

%% Question 3
K_f = 0.140;
V_sum = F_t / K_f;
fprintf("The voltage required to keep the helicopter at steady level is: %d [V]", V_sum)

%% Question 4

%% Question 5

%% Question 6

%% Question 7

%% Question 8