function expr = pend(t,in2,in3)
%PEND
%    EXPR = PEND(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    25-Nov-2022 19:27:28

param16 = in3(:,1);
param17 = in3(:,2);
param18 = in3(:,3);
param19 = in3(:,4);
x_1 = in2(1,:);
x_2 = in2(2,:);
expr = [x_2;-(param17.*param19.*sin(x_1)+param16.*param18.*x_2-1.0)./(param18.*param19)];