%% Feedback Linearization Check 
clear all; close all; clc; 
% ------------------------------------------------------------------------
% Feedback Linearization Checking Steps: 
% 1. Define the System 
% 2. Put system in the form x_dot(t) = f(x) + g(x)*u 
% 3. Check if rank of [g(x0), adfg(x0), adf2g(x0), ..., adf(n-1)g(x0)] = n 
%    You want a linearly independent set of vectors 
% 4. Check involutivity 
%    You want the gradient of x to be in the span of g(x) and adf(n-2)g(x) 
%       (a) Find the Lie Bracket: [g(x), adf(n-2)g]
%       (b) add to span: {g(x), adf(n-2)g(x), [g(x), adf(n-2)g]}
%       (c) check the rank: want rank[g(x), adf(n-2)g] = rank[g(x),
%       adf(n-2)g(x), [g(x), adf(n-2)g]]
% ------------------------------------------------------------------------

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
% U1 = 0; % sym('U1','real'); % APC input 

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
syms x1 x2 x3 real

ad_fng = sym('ad_fng','real'); 

% Define the  system -----------------------------------------------------
f = [z1*x1*x3 + g1*D1 - d1*x1; 
     K2 - a2*D1*x2 + n2*x1 - z2*x2*x3 - d2*x2; 
     g3*D1 - z3*x2*x3 - d3*x3]; 
 
g = [-a1*x1; 
      0; 
      0]; 

h = [0, 0, 1]; 

% Find the equilibrium points --------------------------------------------
sys =  [(z1*x1*x3 + g1*D1 - d1*x1 - a1*x1*U1) == 0, 
        (K2 - a2*D1*x2 + n2*x1 - z2*x2*x3 - d2*x2) == 0,
        (g3*D1 - z3*x2*x3 - d3*x3) == 0]; 

x_eqm = solve(sys,[x1 x2 x3],'Real',true); 
x1_eqm = vpa(x_eqm.x1)
x2_eqm = vpa(x_eqm.x2)
x3_eqm = vpa(x_eqm.x3)

% diseased equilibrium = [0.21, 1.50, 0.014]; (only fully positive answer) 

U1_eqm = 0; % no input 

% relative degree of the system ------------------------------------------
r = 3; % calculated by hand 

%% Controllability/Observability Check (need Controllability for FBL)
% Check in Cont_Obs_3state.m

%%  Step 3: Check rank[g(x0), adfg(x0), adf2g(x0), ..., adf(n-1)g(x0)] = n 
%    You want a linearly independent set of vectors 
X = [x1, x2, x3]; 

% Call the Lie Bracket function ------------------------------------------
ad_fng = liebracket(f,g,X,r) ; % output is: 
% (g, [f,g], [f,[f,g]], [f,[[f,[f,g]]]) == (g, adf1g, adf2g, adf3g) 

% Now sub in equilibrium points: 
ad_fng_eqm = subs(ad_fng,[x1 x2 x3],[x1_eqm(1), x2_eqm(1), x3_eqm(1)]) 

% Check the rank of the adjoint: 
n_rank = rank(ad_fng_eqm) % rank = 3, therefor the vectors are lin. ind.

%% Step 4: Check involutivity

% Now check for involutivity ---------------------------------------------
adf0g = ad_fng(:,1) % this should be g 
adf1g = ad_fng(:,2) % this is [f,g]
adf2g = ad_fng(:,3) % this is [f,[f,g]]
adf3g = ad_fng(:,4) % this is [f,[f2g]]

% Find the Lie bracket of g(x) with adf(n-2)g ----------------------------
adg_adf1g = liebracket(g,adf1g,X,1);  % output is: (adf1g, [g,adf1g])

% add this to span: {g(x), adf(n-2)g(x), [g(x), adf(n-2)g]}
M_invol = [g, adf1g, adg_adf1g(:,2)]
M_invol_eqm = subs(M_invol,[x1 x2 x3],[x1_eqm(1), x2_eqm(1), x3_eqm(1)]) % sub in eqm points 

% sub in equilibrium points 
rank_invol_eqm = rank(M_invol_eqm)
rank_check_eqm = rank(subs([g, adf1g],[x1 x2 x3],[x1_eqm(1), x2_eqm(1), x3_eqm(1)]))
% both ranks match each other so the system is involutive! 


%% Construct the linearizing controller, u 
% u = 1/(LgLf^(n-1)h(x)) * (-Lf^nh(x) + v) % v is the tracking controller
% (Designed separately) 





