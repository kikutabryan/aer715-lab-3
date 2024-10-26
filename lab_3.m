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
fprintf("The voltage required to keep the helicopter at steady level is: %d [V]\n", V_sum)

%% Question 4
% The transfer function from equation 2 is found by subsituting F_t =
% K_f*V_sum into it, then taking the Laplace transforms of both sides.
% Afterwards, T_g is ignored, and E(s)/V_sum(s) is solved for yielding:
% G_elev(S) = E(s)/V_sum(s) = (L_a * K_f) / (J_e * s^2)
numerator = L_a * K_f;
denominator = [J_e, 0, 0];  % Represents J_e * s^2
G2_elev1 = tf(numerator, denominator);
disp("The transfer function G2_elev1 is:");
G2_elev1

%% Question 5
% Load first set of elevation data
elevation_data_1 = load("elevationData1.mat").elev1;
time = elevation_data_1(1, 121:2500);  % Select to t=24.99 [s]
volts1 = elevation_data_1(2, 121:2500);
elev1 = elevation_data_1(3, 121:2500);

% Load second set of elevation data
elevation_data_2 = load("elevationData2.mat").elev2;
volts2 = elevation_data_2(2, 121:2500);
elev2 = elevation_data_2(3, 121:2500);

% Load third set of elevation data
elevation_data_3 = load("elevationData3.mat").elev3;
volts3 = elevation_data_3(2, 121:2500);
elev3 = elevation_data_3(3, 121:2500);

% Create plot
plot(time, elev1, "-", time, elev2, "-", time, elev3, "-");
xlabel("Time (s)");
ylabel("Elevation");
title("Elevation vs Time");
legend("elev1", "elev2", "elev3");
grid on;

% Create iddata objects for each elevation dataset
data_1 = iddata(elev1', volts1', 0.01);
data_2 = iddata(elev2', volts2', 0.01);
data_3 = iddata(elev3', volts3', 0.01);

% Combine the iddata objects into one dataset
combined_data = merge(data_1, data_2, data_3);

% Use tfest to estimate a transfer function with 2 poles
G_elev2 = tfest(combined_data, 2, 1); % 1 zero
disp("The transfer function G2_elev2 is:");
G_elev2

% Use tfest to estimate a transfer function with 3 poles
G_elev3 = tfest(combined_data, 3, 1); % 1 zero
disp("The transfer function G2_elev3 is:");
G_elev3

% Extract numerator and denominator polynomials for G_elev2
[num2, den2] = tfdata(G_elev2, 'v'); % 'v' for vector format
disp('G_elev2 Numerator:');
disp(num2);
disp('G_elev2 Denominator:');
disp(den2);

% Extract numerator and denominator polynomials for G_elev3
[num3, den3] = tfdata(G_elev3, 'v'); % 'v' for vector format
disp('G_elev3 Numerator:');
disp(num3);
disp('G_elev3 Denominator:');
disp(den3);

%% Question 6

%% Question 7

%% Question 8