function [Ym,Yw,G,Im,Iw,Ig,D] = HIV_model(Ym0,Yw0,G0,Im0,Iw0,Ig0,D0,T,dt)
    % if delta = 0 we assume a model without immunity loss
    Ym = zeros(1,T/dt);
    Ym(1) = Ym0;
    Yw = zeros(1,T/dt);
    Yw(1) = Yw0;
    G = zeros(1,T/dt);
    G(1) = G0;
    Im = zeros(1,T/dt);
    Im(1) = Im0;
    Iw = zeros(1,T/dt);
    Iw(1) = Iw0;
    Ig = zeros(1,T/dt);
    Ig(1) = Ig0;  
    D = zeros(1,T/dt);
    D(1) = D0; 
    
    tau1=0;
    tau2=0;
    tau3=0;
    zeta=0.09;
    gamma1=0;
    gamma2=0;
    omega1=0.00144;
    omega2=0.00288;
    omega4=0.04986;
    d1=0.007165862;
    d2=5.5047E-05;
    M=246538;
    W=257026;
    
    for tt = 1:(T/dt)-1
        % Equations of the model
        dYm = M-(1-tau3)*zeta*Ym(tt)+tau3*G(tt)-omega1*(1-gamma1)*(1-tau1)*((Iw(tt))/(Iw(tt)+Yw(tt)))*Ym(tt)-d1*Ym(tt);
    
        dYw =  W-omega2*(1-gamma1)*(1-tau1)*((Im(tt))/(Im(tt)+Ym(tt)))*Yw(tt)-d1*Yw(tt);
    
        dG =  (1-tau3)*zeta*Ym(tt)-tau3*G(tt)-(1-gamma2)*(1-tau2)*(omega4)*(Ig(tt)/(Ig(tt)+G(tt)))*G(tt)-d1*G(tt);
    
        dIm= omega1*(1-gamma1)*((Iw(tt))/(Iw(tt)+Yw(tt)))*Ym(tt)-d2*Im(tt);
    
        dIw = omega2*(1-gamma1)*(1-tau1)*((Im(tt))/(Im(tt)+Ym(tt)))*Yw(tt)-d2*Iw(tt);                                                              
    
        dIg = (1-gamma2)*(1-tau2)*(omega4)*(Ig(tt)/(Ig(tt)+G(tt)))*G(tt)-d2*Ig(tt);             
    
        dD = d1*Ym(tt)+d1*Yw(tt)+d1*G(tt)+d2*Im(tt)+d2*Iw(tt)+d2*Ig(tt);
        
        Ym(tt+1) = Ym(tt) + dYm;
        Yw(tt+1) = Yw(tt) + dYw;
        G(tt+1) = G(tt) + dG;
        Im(tt+1) = Im(tt) + dIm;
        Iw(tt+1) = Iw(tt) + dIw;
        Ig(tt+1) = Ig(tt) + dIg;
        D(tt+1) = D(tt) + dD;
 
    end
end