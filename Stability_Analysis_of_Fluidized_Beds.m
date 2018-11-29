clc;
clear all:
format long;
syms T Tp P Pp
%% finding Steady State Tp
figure(1)
fplot(@(Tp) 0.48*exp(20.7-15000/Tp)/(0.1926*exp(20.7-15000/Tp)+1),[500 1000])
hold on
fplot(@(Tp) (2.6*Tp-1752)/269.267,[500 1000])
hold off
xlabel('Tp')
ylabel('Q1 and Q11')
title('Q1 and Q11 vs Tp')
legend('Q11','Q1')
z = @(Tp) (0.48*exp(20.7-15000/Tp)/(0.1926*exp(20.7-15000/Tp)+1))-(2.6*Tp-1752)/269.267;
Tp1 = fzero(z,[650,725]);
Tp2 = fzero(z,[725,800]);
Tp3 = fzero(z,[900,1000]);
%% finding First Steady State Values 
eq1 = 320*Pp-321*P +0.1;
eq2 = 1752 + 266.667*Tp1-269.267*T;
eq3 = 1866.76*P-(1.12*exp(20.7-15000/Tp1)+1866.76)*Pp;
sol = solve(eq1,eq2,eq3);  
T1=double(sol.T);
P1=double(sol.P);
Pp1=double(sol.Pp);
%% finding Second Steady State Values 
eq1 = 320*Pp-321*P +0.1;
eq2 = 1752 + 266.667*Tp2-269.267*T;
eq3 = 1866.76*P-(1.12*exp(20.7-15000/Tp2)+1866.76)*Pp;
sol = solve(eq1,eq2,eq3);  
T2=double(sol.T);
P2=double(sol.P);
Pp2=double(sol.Pp);
%% finding Third Steady State Values 
eq1 = 320*Pp-321*P +0.1;
eq2 = 1752 + 266.667*Tp3-269.267*T;
eq3 = 1866.76*P-(1.12*exp(20.7-15000/Tp3)+1866.76)*Pp;
sol = solve(eq1,eq2,eq3);  
T3=double(sol.T);
P3=double(sol.P);
Pp3=double(sol.Pp);
%% Steady State Values
SS = [P1 P2 P3;T1 T2 T3;Pp1 Pp2 Pp3;Tp1 Tp2 Tp3];
%% Jacobian Matrix
E1 = 320*Pp-321*P +0.1;
E2 = 1752 + 266.667*Tp-269.267*T;
E3 = 1866.76*P-(1.12*exp(20.7-15000/Tp)+1866.76)*Pp;
E4 = (6.22*exp(20.7-15000/Tp))*Pp+1.296*T-1.296*Tp;
J = jacobian([E1; E2; E3; E4], [P T Pp Tp]);
J1 = subs(J,{P,T,Pp,Tp},{P1,T1,Pp1,Tp1});
Lambda1 = double(eig(J1));
J2 = subs(J,{P,T,Pp,Tp},{P2,T2,Pp2,Tp2});
Lambda2 = double(eig(J2));
J3 = subs(J,{P,T,Pp,Tp},{P3,T3,Pp3,Tp3});
Lambda3 = double(eig(J3));
%% Steady State Values
Pss = [P1 P2 P3]
Tss = [T1 T2 T3]
Ppss = [Pp1 Pp2 Pp3]
Tpss = [Tp1 Tp2 Tp3]
Lambda1ss = [Lambda1(1) Lambda2(1) Lambda3(1)]
Lambda2ss = [Lambda1(2) Lambda2(2) Lambda3(2)]
Lambda3ss = [Lambda1(3) Lambda2(3) Lambda3(3)]
Lambda4ss = [Lambda1(4) Lambda2(4) Lambda3(4)]


%% Plot for First Steady state
tspan1 = [0 1500];
y0 = [0.09352 0.09350 690.4455 690.6074];
[t,y] = ode15s(@(t,y) odefcn(t,y), tspan1, y0);
figure(2)
subplot(3,2,1)
plot(t,y(:,1),'-+r',t,y(:,2),'k')
grid on
xlabel('Time')
ylabel('P and Pp')
title('P and Pp vs time for First Steady State')
legend('P','Pp')
subplot(3,2,2)
plot(t,y(:,3),'-+b',t,y(:,4),'m')
grid on
xlabel('Time')
ylabel('T and Tp')
title('T and Tp vs time for First Steady State')
legend('T','Tp')
%% Plot for Second Steady State
y0 = [0.06746 0.06694 758.3454 759.0693 ];
[t,y] = ode15s(@(t,y) odefcn(t,y), tspan1, y0);
subplot(3,2,3)
plot(t,y(:,1),'-+r',t,y(:,2),'k')
grid on
xlabel('Time')
ylabel('P and Pp')
title('P and Pp vs time for Second Steady State')
legend('P','Pp')
subplot(3,2,4)
plot(t,y(:,3),'-+b',t,y(:,4),'m')
grid on
xlabel('Time')
ylabel('T and Tp')
title('T and Tp vs time for Second Steady State')
legend('T','Tp')

%% Plot for Third Steady State
y03 = [0.06746 0.06694 758.3454 761.1693];
[t,y] = ode15s(@(t,y) odefcn(t,y), tspan1, y03);
subplot(3,2,5)
plot(t,y(:,1),'-+r',t,y(:,2),'k')
grid on
xlabel('Time')
ylabel('P and Pp')
title('P and Pp vs time for Third Steady State')
legend('P','Pp')
subplot(3,2,6)
plot(t,y(:,3),'-+b',t,y(:,4),'m')
grid on
xlabel('Time')
ylabel('T and Tp')
title('T and Tp vs time for Third Steady State')
legend('T','Tp')

%% function for ODEs
function dydt = odefcn(t,y)
dydt = zeros(4,1);
dydt(1) = 0.1 + 320*y(2)-321*y(1);
dydt(2) = -(1.12*exp(20.7-15000/y(4))+1866.76)*y(2)+1866.76*y(1);
dydt(3) = 1752 + 266.67*y(4)-269.267*y(3);
dydt(4) = 6.22*exp(20.7-15000/y(4))*y(2)+1.296*y(3)-1.296*y(4);
end

%% The Numerical values from the plot

%Parameters    % First Steady State    % Second Steady State   %Third Steady State
%P-Numerical =           0.09319            0.06692                0.006377
%T-Numerical =            692.2              759.2                   914.2
%Pp-Numerical=           0.09322            0.06692                0.006377
%Tp-Numerical=            692                759.2                   916.4
