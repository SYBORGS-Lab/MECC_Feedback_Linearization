%% Controllability and Observability Check 
clear all; close all; clc; 

% Define the symbolic variables 
% z1 = sym('z1','real'); 
% g1 = sym('g1','real');
% a1 = sym('a1','real');
% d1 = sym('d1','real');
% K2 = sym('K2','real');
% a2 = sym('a2','real');
% n2 = sym('n2','real');
% z2 = sym('z2','real');
% d2 = sym('d2','real');
% g3 = sym('g3','real');
% z3 = sym('z3','real');
% d3 = sym('d3','real');
% D1 = sym('D1','real'); % IL6 disturbance 
% U1 = sym('U1','real'); % APC input 

z1 = 0.8;  % 1/ nM min      -- activation/upregulation of Xa_Va by Fbn 
g1 = 10;   % 1/min      -- upregulation of Xa_Va by IL6
a1 = 1;    % 1/nM min   -- downregulation by APC input and Xa_Va 
d1 = 0.12; % 1/min      -- natural degradation of Xa_Va -- based on HM for prothrombinase 
K2 = 0.034 * log(2)/4; % equals 0.0059;  % nM/min     -- natural production of tPA
a2 = 15;   % 1/nM min   -- downregulation by IL6 and tPA
n2 = 1.5;  % 1/nM min   -- upregulation by APC input and tPA
z2 = 0.5;  % 1/nM min 
d2 = log(2)/4; % equals 0.1733;  % 1/min      -- natural degradation of tPA 
g3 = 20;   % 1/min      -- upregulation by Xa_Va
z3 = 0.01; % 1/nM min 
d3 = 0.0046; % 1/min      -- natural degradation of fibrin  

D1 = 0.0023; % nM, disease case 
% D1 = 0; % no disturbance in healthy case
U1 = 0; % no APC input 

% x1 = sym('x1','real'); % Xa_Va
% x2 = sym('x2','real'); % tPA 
% x3 = sym('x3','real'); % Fbn 
syms x1 x2 x3 U1 positive

% system: 
sys =  [(z1*x1*x3 + g1*D1 - d1*x1 - a1*x1*U1), 
        (K2 - a2*D1*x2 + n2*x1 - z2*x2*x3 - d2*x2),
        (g3*D1 - z3*x2*x3 - d3*x3)]; 

% Symbolic variables 
% x1_eqm = sym('x1_eqm','positive');
% x2_eqm = sym('x2_eqm','positive');
% x3_eqm = sym('x3_eqm','positive');    
% U1_eqm = sym('U1_eqm','positive');  

% Disease state equilibrium values (nM)
x1_eqm = 5.836; % 0.21;
x2_eqm = 31.248; % 1.50;
x3_eqm = 0.145; % 0.014;    
U1_eqm = 0; 

% Healthy state equilibrium values (nM)
% x1_eqm = 0;
% x2_eqm = 0.034;
% x3_eqm = 0;    
% U1_eqm = 0; 

a = [x1_eqm, x2_eqm, x3_eqm, U1_eqm]; 
x = [x1, x2, x3, U1]; 

% first need to linearize the system using taylor series expansion 
taylor = taylor(sys, x, a, 'Order',2)
  
% x1_dotL = D1*g1 - (x1 - x1_eqm)*(d1 + U1*a1 - x3_eqm*z1) - d1*x1_eqm + x1_eqm*z1*(x3 - x3_eqm) - U1*a1*x1_eqm + x1_eqm*x3_eqm*z1
% x2_dotL = K2 - (x2 - x2_eqm)*(d2 + D1*a2 + x3_eqm*z2) - d2*x2_eqm + n2*x1_eqm + n2*(x1 - x1_eqm) - x2_eqm*z2*(x3 - x3_eqm) - D1*a2*x2_eqm - x2_eqm*x3_eqm*z2
% x3_dotL = D1*g3 - (d3 + x2_eqm*z3)*(x3 - x3_eqm) - d3*x3_eqm - x3_eqm*z3*(x2 - x2_eqm) - x2_eqm*x3_eqm*z3
% x = [x1; x2; x3]; 
  

% put in x_dot = A*x + B*u + D, y = C*x form: 
A = [(-d1 - U1_eqm * a1 + x3_eqm * z1),   0,                            x1_eqm * z1; 
       n2,                              (-d2 - D1 * a2 - x3_eqm * z2), (-x2_eqm * z2); 
       0,                               (-x3_eqm * z3),               (-d3 - x2_eqm * z3)]; 

B = [-a1*x1_eqm; 0; 0]; 
    

D = [(D1*g1 + U1_eqm * a1 * x1_eqm - x1_eqm * x3_eqm * z1); 
    (K2 + x2_eqm * x3_eqm * z2); 
    (D1 * g3 + x2_eqm * x3_eqm * z3)]; 

C = [0 0 1]; % only output is fibrin, x3 -- C is a 1x3 matrix 

Cont = [B, A*B, (A*A)*B]
Cont_rank = rank(Cont) 

Obs = [C; C*A; C*(A*A)]
Obs_rank = rank(Obs)    
    
% check stability 
eigA = eig(A) % all eigen values are negative therefore the system is stable at the diseased equilibrium point 