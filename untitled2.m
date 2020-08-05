% Model parameters
beta = 5*10^-9; % rate of infection
gamma = 0.12; % rate of recovery (try also 0.07)
delta = 0.0; % rate of immunity loss
N = 6*10^7; % Total population N = S + I + R
I0 = 10; % initial number of infected
T = 300; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)
fprintf('Value of parameter R0 is %.2f',N*beta/gamma)
% Calculate the model
[S,I,R] = sir_model(beta,gamma,delta,N,I0,T,dt);
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2); grid on;
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R');
% Map
plot(I(1:(T/dt)-1),I(2:T/dt),"LineWidth",1,"Color",'r');
hold on; grid on;
plot(I(2),I(1),'ob','MarkerSize',4);
xlabel('Infected at time t'); ylabel('Infected at time t+1');
hold off;

delta = 1/60; % rate of immunity loss
% Calculate the model
[S,I,R] = sir_model(beta,gamma,delta,N,I0,T,dt);
% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2); grid on;
xlabel('Days');ylabel('Number of individuals');
legend('S','I','R');
