function expr = Model(t,in2,in3)
%Model
%    EXPR = Model(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    01-Dec-2022 15:08:47

A = in2(5,:);
C = in2(6,:);
D = in2(2,:);
E = in2(7,:);
I = in2(3,:);
Q = in2(8,:);
R = in2(1,:);
S = in2(4,:);
param8 = in3(:,1);
param9 = in3(:,2);
param10 = in3(:,3);
param11 = in3(:,4);
param12 = in3(:,5);
param13 = in3(:,6);
param14 = in3(:,7);
param15 = in3(:,8);
param16 = in3(:,9);
param17 = in3(:,10);
param18 = in3(:,11);
t2 = A.*param18;
t3 = I.*param10;
t4 = Q.*param14;
t5 = Q.*param15;
t6 = E.^2;
t7 = S.^2;
t8 = A.*S.*param12;
t9 = A.*S.*param17;
t10 = C.*S.*param12;
t11 = D.*S.*param12;
t12 = E.*S.*param12;
t13 = I.*S.*param9;
t14 = I.*S.*param12;
t15 = Q.*S.*param12;
t16 = R.*S.*param12;
t18 = A+C+D+E+I+Q+R+S;
t17 = param12.*t7;
t19 = 1.0./t18;
mt1 = [t5;t4;t2-t3+E.*param13;-t19.*(t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+param8.*t7+A.*S.*param8+C.*S.*param8+D.*S.*param8+E.*S.*param8+I.*S.*param8+Q.*S.*param8+R.*S.*param8);-t2+E.*param11;-C.*param16+S.*param8];
mt2 = [t19.*(t8+t9+t10+t11+t12+t13+t14+t15+t16+t17-param11.*t6-param13.*t6+C.^2.*param16+A.*C.*param16-A.*E.*param11-A.*E.*param13+C.*D.*param16-C.*E.*param11-C.*E.*param13+C.*E.*param16-D.*E.*param11-D.*E.*param13+C.*I.*param16-E.*I.*param11-E.*I.*param13+C.*Q.*param16+C.*R.*param16+C.*S.*param16-E.*Q.*param11-E.*Q.*param13-E.*R.*param11-E.*R.*param13-E.*S.*param11-E.*S.*param13);t3-t4-t5];
expr = [mt1;mt2];
