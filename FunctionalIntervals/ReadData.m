addpath ../COVID19/Results/Estimations/
addpath ../Dengue/Results/Estimations/


%% Dengue estimations
load('../Dengue/Results/Estimations/Estimations_Itagui.mat');

Rows = ["M_e(0)"; 'H_r(0)'; 'H_s(0)'; 'H_e(0)'; 'M_s(0)'; 'M_i(0)';...
         '\beta_h';'\beta_m'; '\gamma_h';'\Lambda'; '\mu_m'; '\theta_h'; '\theta_m'];

csvwrite('DengueEstimations.csv', EstimatedParams(1:1100,:))

%% COVID-19 estimationsclear

load("COVID19/Results/Estimations/Estimations_Amazonas.mat")
csvwrite('COVIDEstimations.csv', EstimatedParams)
