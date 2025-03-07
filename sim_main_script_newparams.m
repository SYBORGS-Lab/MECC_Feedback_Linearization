% Simulation main script

%% Multiple Initial Conditions - no control input 
clear all; close all; clc; 

tic 
% virtual input gain
gain.k1 = 0; 
gain.k2 = 0; 
gain.k3 = 0; 

% set u = 0 and v = 0 using Simulink file 

% Define a set of initial conditions (each row represents [a, b, c])
initial_conditions = [0, 0,  0
                      1, 3,  0.05, 
                      4, 12, 0.2, 
                      6, 32, 0.5];

% Load the Simulink model
model_name = "Feedback_Lin_MECC2025";
load_system(model_name);

% Simulation stop time
stop_time = "60";
% Initialize figure
figure();

% Loop through each initial condition
for i = 1:size(initial_conditions, 1)
    % Set the initial condition
    x0 = initial_conditions(i, :);

    % Assign the initial condition to the MATLAB base workspace
    assignin('base', 'x0', x0);

    % Run the simulation
    simout = sim(model_name, 'StopTime', stop_time);

    % Extract time and state data
    t = simout.tout;         % Time vector
    x = simout.x.Data;       % State matrix
    x1 = x(:,1);             % First state
    x2 = x(:,2);             % Second state
    x3 = x(:,3);             % Third state

    % Format for short legend
    legend_str = sprintf('%.3f, %.3f, %.3f', x0(1), x0(2), x0(3));

    % Plot x1 over time
    subplot(3,1,1); hold on;
    plot(t, x1, 'LineWidth', 3, 'DisplayName', legend_str);
    yline(5.836, 'k--', 'HandleVisibility','off', 'LineWidth', 3); 
    ylabel('x_1 [nM]','Fontsize', 14);
    title('State Response x_1 vs Time','Fontsize', 14);

    % Plot x2 over time
    subplot(3,1,2); hold on;
    plot(t, x2, 'LineWidth', 3, 'DisplayName', legend_str);
    yline(31.248, 'k--', 'HandleVisibility','off', 'LineWidth', 3); 
    ylabel('x_2 [nM]','Fontsize', 14);
    title('State Response x_2 vs Time','Fontsize', 14);

    % Plot x3 over time
    subplot(3,1,3); hold on;
    plot(t, x3, 'LineWidth', 3, 'DisplayName', legend_str);
    yline(0.145, 'k--', 'HandleVisibility','off', 'LineWidth', 3);
    xlabel('Time (min)','Fontsize', 14);
    ylabel('x_3 [nM]','Fontsize', 14);
    title('State Response x_3 vs Time','Fontsize', 14);

    % Store results
    results{i} = struct('time', t, 'x1', x1, 'x2', x2, 'x3', x3);

    % Display progress
    fprintf('Simulation %d completed with x0 = [%.3f, %.3f, %.3f]\n', i, x0(1), x0(2), x0(3));
end

% Add legends to plots
subplot(3,1,1); legend show;
subplot(3,1,2); legend show;
subplot(3,1,3); legend show;
toc

%% Multiple Gains (Tracking Controller) 

% set u = v, and let v = the linear error tracking in Simulink file 

clear all; close all; clc; 

tic 
% Initial Condition 
x0 = 4.*[1, 3, 0.05]; 

% Load the Simulink model
model_name = "Feedback_Lin_MECC2025";
load_system(model_name);

% Simulation stop time
stop_time = "60";

% Define the range of gains to test

gain_values = [1, 3, 2; 
               4, 3, 2; 
               4, 3, 4; 
               2, 10, 1/2; 
               2, 10, 2; 
               4, 10, 1/2]; 

% Initialize figure
figure();

% Loop through each set of gain values
for i = 1:size(gain_values, 1)
    % Set gain values
    gain.k1 = gain_values(i, 1);
    gain.k2 = gain_values(i, 2);
    gain.k3 = gain_values(i, 3);

    % Assign gains and initial condition to MATLAB base workspace
    assignin('base', 'x0', x0);
    assignin('base', 'gain', gain);

    % Run the simulation
    simout = sim(model_name, 'StopTime', stop_time);

    % Extract time and state data
    t = simout.tout;         % Time vector
    x = simout.x.Data;       % State matrix
    x1 = x(:,1);             % First state
    x2 = x(:,2);             % Second state
    x3 = x(:,3);             % Third state

    e = simout.e.Data;   % Error 
    yd = simout.yd.Data; % Desired Fibrin 
    u = simout.u.Data;   % Linearizing Control 
    v = simout.v.Data;   % Tracking Control 
    
    LgLf2h = simout.LgLf2h.Data; % denominator of linearizing controller 
    Lf3h = simout.Lf3h.Data; % numerator of linearizing controller (negative part before v) 
    % Format for short legend (one decimal place)
    legend_str = sprintf('k1=%.1f, k2=%.1f, k3=%.1f', gain.k1, gain.k2, gain.k3);

    % Plot x1 over time
    subplot(4,1,1); hold on;
    plot(t, x1, 'LineWidth', 3, 'DisplayName', legend_str);
    ylabel('x_1 [nM]','Fontsize', 14);
    title('State Response x_1 vs Time','Fontsize', 14);

    % Plot x2 over time
    subplot(4,1,2); hold on;
    plot(t, x2, 'LineWidth', 3, 'DisplayName', legend_str);
    ylabel('x_2 [nM]','Fontsize', 14);
    title('State Response x_2 vs Time','Fontsize', 14);

    % Plot x3 over time
    subplot(4,1,3); hold on;
    plot(t, x3, 'LineWidth', 3, 'DisplayName', legend_str);
    plot(t, yd, 'k--', 'LineWidth', 3, 'HandleVisibility','off'); % Desired output
    ylabel('x_3 [nM]','Fontsize', 14);
    title('State Response x_3 vs Time','Fontsize', 14);

     % Plot e over time
    subplot(4,1,4); hold on;
    plot(t, e, 'LineWidth', 3, 'DisplayName', legend_str);
%     axis([0 60 0 0.03]); 
    ylabel('e','Fontsize', 14);
    xlabel('Time (min)','Fontsize', 14);
    title('Error Signal e vs Time','Fontsize', 14);
    
    % Store results
    results{i} = struct('time', t, 'x1', x1, 'x2', x2, 'x3', x3, 'yd', yd, 'e', e, 'u', u, 'v', v, 'LgLf2h', LgLf2h, 'Lf3h', Lf3h);

    % Display progress
    fprintf('Simulation %d completed with k1=%.1f, k2=%.1f, k3=%.1f\n', i, gain.k1, gain.k2, gain.k3);
end

% Add legends to plots
subplot(4,1,1); legend show;
% subplot(4,1,2); legend show;
% subplot(4,1,3); legend show;
% subplot(4,1,4); legend show;

% New figure for e, u, v
figure();

% Loop again to plot e, u, and v
for i = 1:size(gain_values, 1)
    t = results{i}.time;
    e = results{i}.e; % error is: e = x3 - yd
    u = results{i}.u;
    v = results{i}.v;
    legend_str = sprintf('k1=%d, k2=%d, k3=%d', gain_values(i, 1), gain_values(i, 2), gain_values(i, 3));

    % Plot e over time
    subplot(3,1,1); hold on;
    plot(t, e, 'LineWidth', 3, 'DisplayName', legend_str);
%     axis([0 60 0 0.03]); 
    ylabel('e');
    title('Error Signal e vs Time','Fontsize', 14);

    % Plot u over time
    subplot(3,1,2); hold on;
    plot(t, u, 'LineWidth', 3, 'DisplayName', legend_str);
    ylabel('u','Fontsize', 14);
    title('Linearizing Control Input u vs Time','Fontsize', 14);

    % Plot v over time
    subplot(3,1,3); hold on;
    plot(t, v, 'LineWidth', 3, 'DisplayName', legend_str);
    xlabel('Time (min)','Fontsize', 14);
    ylabel('v','Fontsize', 14);
    title('Tracking Control Input v vs Time','Fontsize', 14);
end

% Add legends to e, u, v plots
subplot(3,1,1); legend show;
subplot(3,1,2); legend show;
subplot(3,1,3); legend show;


% New figure for Linearizing controller denominator 
figure();

% Loop again to plot e, u, and v
for i = 1:size(gain_values, 1)
    t = results{i}.time;
    LgLf2h = results{i}.LgLf2h;
    Lf3h = results{i}.Lf3h;
    
    combined = Lf3h./LgLf2h; 
    
    legend_str = sprintf('k1=%d, k2=%d, k3=%d', gain_values(i, 1), gain_values(i, 2), gain_values(i, 3));

    % Plot denomintor over time (check for if it crosses zero) 
    subplot(3,1,1); hold on; 
    plot(t, LgLf2h, 'LineWidth', 3, 'DisplayName', legend_str);
    ylabel('LgLf2h','Fontsize', 14);
    title('Linearizing Controller Denominator vs Time','Fontsize', 14);
    
    % plot numerator over time (compare to value of v) 
    subplot(3,1,2); hold on; 
    plot(t, Lf3h,'Linewidth', 3, 'DisplayName', legend_str); 
    xlabel('Time (min)','Fontsize', 14);
    ylabel('Lf3h','Fontsize', 14); 
    title('Linearizing Controller Numerator vs Time','Fontsize', 14); 
    
     % plot multiplied values 
    subplot(3,1,3); hold on; 
    plot(t, combined,'Linewidth', 3, 'DisplayName', legend_str); 
    xlabel('Time (min)','Fontsize', 14);
    ylabel('Lf3h/LgLf2h','Fontsize', 14); 
    title('Numerator/Denomintor vs Time','Fontsize', 14); 
end 
subplot(3,1,1); legend show;
subplot(3,1,2); legend show;
subplot(3,1,3); legend show;

toc

%% Final graphs 
clear all; close all; clc; 

% Note: all concentrations are in nM and time is in minutes 

% set u = FBL code, make sure v = error tracking 

% Initial Condition 
x0 = 4.*[1, 3, 0.05]; 

% virtual input gain
gain.k1 = 2; 
gain.k2 = 10; 
gain.k3 = 0.5; 

simout = sim("Feedback_Lin_MECC2025.slx", StopTime="60");

% plots
x = simout.x.Data;
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

e = simout.e.Data; 
t = simout.tout; 

yd = simout.yd.Data; 

figure(); 
subplot(4, 1, 1) 
plot(t, x1, 'LineWidth', 3, 'DisplayName', 'x_1')
grid on
ylabel("x_1 [nM]", 'Fontsize', 14)
title('Prothrombinase Response', 'Fontsize', 14) 

subplot(4,1,2) 
plot(t, x2, 'LineWidth', 3, 'DisplayName', 'x_2')
grid on
ylabel("x_2 [nM]", 'Fontsize', 14)
title('tPA Response', 'Fontsize', 14) 

subplot(4,1,3) 
hold on 
plot(t, x3, 'LineWidth', 3, 'DisplayName', 'x_3')
plot(t, yd, '--', 'LineWidth', 3, 'DisplayName', 'y_d')
grid on
ylabel("x_3 [nM]", 'Fontsize', 14)
legend('x_3','y_d','Location', 'Southeast')
title('Fibrin Response', 'Fontsize', 14) 

subplot(4,1,4) 
plot(t, e, 'LineWidth', 3, 'DisplayName', 'e')
grid on
ylabel("e", 'Fontsize', 14)
title('Tracking Error', 'Fontsize', 14) 
xlabel('Time (min)', 'Fontsize', 14) 