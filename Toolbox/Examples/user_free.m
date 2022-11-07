function odes = user_free(params)
%R_0
%    ODES = R_0(BETA_M,BETA_H,THETA_M,THETA_H,ALPHA,MU_M,GAMMA_H,MU_H)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    04-Aug-2019 17:48:42
beta_m=params(1,:);
beta_h=params(2,:);
theta_m=params(3,:);
theta_h=params(3,:);
alpha=params(3,:);
mu_m=params(3,:);
gamma_h=params(3,:);
mu_h=params(3,:);

t2 = alpha.*mu_m;
odes = (beta_h.*beta_m.*theta_h.*theta_m)./(t2.*(gamma_h+mu_h).*(mu_h+theta_h).*(t2+theta_m));
