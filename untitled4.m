% Model parameters
Ym= 16380526;
Yw =16576010;
G =131045;
Im = 1951;
Iw = 7;
Ig = 5854;
D = 236162;

T = 27; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)

[Ym,Yw,G,Im,Iw,Ig,D] = HIV_model(Ym,Yw,G,Im,Iw,Ig,D,T,dt);

plot(tt,Ym,'b',tt,Yw,'r',tt,G,'g',Im,'y',Iw,'m',Ig,'c',D,'k','LineWidth',2); grid on;
xlabel('Days');ylabel('Number of individuals');
legend('Ym','Yw','G','Im','Iw','Ig','D');