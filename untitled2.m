doc ode45

function dzdt = rhs(t,z,par)
% Right-hand side of differential equations
% par: Parameter struct
    % Unpack vector of dependent variables:
    
    N_M(0) = z(1);
    N_W(0) = z(2);
    N_G(0) = z(3);
    S_M(0) = z(4);
    S_W(0) = z(5);
    S_G(0) = z(6);
    ID_M(0) = z(7);
    ID_W(0) = z(8);
    ID_G(0) = z(9);
    UD_M(0) = z(10);
    UD_W(0) = z(11);
    UD_G(0) = z(12);
    A_M(0) = z(13);
    A_W(0) = z(14);
    A_G(0) = z(15);
    
    
    dN_Mdt = (N_M + N_G) * par.c(12) - par.c(1) * S_M + par.c(4) * (S_G + UD_G + ID_G + A_G) - par.c(14) * ( S_M + ID_M + UD_M + A_M) -  par.c(15) * A_M;
    
    dN_Wdt = N_W * par.c(12) - par.c(14) * ( S_W + ID_W + UD_W + A_W) -  par.c(15) * A_W;
    
    dN_Gdt = par.c(1) * S_M  - par.c(4) * (S_G + UD_G + ID_G + A_G) - par.c(14) * ( S_G + ID_G + UD_G + A_G) -  par.c(15) * A_G;
    
   
    dS_Mdt = (N_M + N_G) * par.c(12) - par.c(1) * S_M + par.c(4) * S_G - (1-par.c(2)) * (1-par.c(5)) * (1-par.c(7)) * par.c(10) * par.c(13) * UD_W * S_M /N_W - (1-par.c(2)) * (1-par.c(5)) * (1-par.c(7)) * par.c(10) * (1-par.c(13)) * UD_W * S_M /N_W - par.c(14) * S_M;
    
    dS_Wdt = N_W * par.c(12) - (1-par.c(2)) * (1-par.c(5)) * (1-par.c(8)) * par.c(11) * par.c(13) * UD_M * S_W /N_M - (1-par.c(2)) * (1-par.c(5)) * (1-par.c(8)) * par.c(10) * (1-par.c(13)) * UD_M * S_W /N_M - par.c(14) * S_W;
    
    dS_Gdt = par.c(1) * S_M - par.c(4) * S_G - (1-par.c(2)) * (1-par.c(3)) * (1-par.c(6)) * (1-par.c(9)) * par.c(12) * par.c(13) * UD_G * S_G /N_G - (1-par.c(2)) * (1-par.c(3)) * (1-par.c(6)) * (1-par.c(9)) * par.c(12) * (1-par.c(13)) * UD_G * S_G /N_G - par.c(14) * S_G;
    
    
    dID_Mdt = (1-par.c(2)) * (1-par.c(5)) * (1-par.c(7)) * par.c(10) * par.c(13) * UD_W * S_M /N_W + par.c(13) * UD_M + par.c(17) * A_M - par.c(16) * ID_M + par.c(4) * ID_G - par.c(14) * ID_M;
    
    dID_Wdt = (1-par.c(2)) * (1-par.c(5)) * (1-par.c(8)) * par.c(10) * par.c(13) * UD_M * S_W /N_M + par.c(13) * UD_W + par.c(17) * A_W - par.c(16) * ID_W - d1 * ID_M;
    
    dID_Gdt = (1-par.c(2)) * (1-par.c(3)) * (1-par.c(6)) * (1-par.c(8)) * par.c(11) * par.c(13) * UD_G * S_G /N_G + par.c(13) * UD_G + par.c(17) * A_G - par.c(16) * ID_M - par.c(4) * ID_G - par.c(14) * ID_G;
  
    
    dUD_Mdt = (1-par.c(2)) * (1-par.c(5)) * (1-par.c(7)) * par.c(10) * (1-par.c(13)) * UD_W * S_M /N_W - par.c(13) * UD_M - par.c(16) * ID_M + par.c(16) * UD_G - par.c(14) * UD_M;
    
    dUD_Wdt = (1-par.c(2)) * (1-par.c(5)) * (1-par.c(8)) * par.c(10) * (1-par.c(13)) * UD_M * S_W /N_M - par.c(13) * UD_W  - par.c(16) * ID_W - par.c(14) * UD_W;
    
    dUD_Gdt = (1-par.c(2)) * (1-par.c(3)) * (1-par.c(6)) * (1-par.c(8)) * par.c(11) * (1-par.c(13)) * UD_G * S_G /N_G - par.c(13) * UD_G  - par.c(16) * UD_G - par.c(4) * UD_G - par.c(14) * UD_G;
    
    
    dA_Mdt = par.c(16) * (UD_M + ID_M) - par.c(17) * A_M - (par.c(14)+ par.c(15)) * A_M + par.c(4) * A_G ;
    
    dA_Wdt = par.c(16) * (UD_W + ID_W) - par.c(17) * A_W - (par.c(14)+ par.c(15)) * A_W;
    
    dA_Gdt = par.c(16) * (UD_G + ID_G) - par.c(17) * A_G - (par.c(14)+ par.c(15)) * A_G - par.c(4) * A_G ;
    
    
    % In the same way for the remaining equation
    dzdt = [dN_Mdt;dN_Wdt;dN_Gdt;dS_Mdt;dS_Wdt;dS_Gdt;dID_Mdt;dID_Wdt;dID_Gdt;dUD_Mdt;dUD_Wdt;dUD_Gdt;dA_Mdt;dA_Wdt;dA_Gdt];
end


par.c = [];  

% zeta, tau1, tau2,tau3, gamma1, gamma2, omega1, omega2, omega3, eta1,eta2, kappa, phi, d1,d2, beta, theta




N_M0 = 0; N_W0 = 0; N_G0 = 0; S_M0 = 0; S_W0 = 0; S_G0 = 0; ID_M0 = 0; ID_W0 = 0; ID_G0 = 0; UD_M0 = 0; UD_W0 = 0; UD_G0 = 0; A_M0 = 0; A_W0 = 0; A_G0 = 0;
z0 = [N_M0; N_W0; N_G0; S_M0; S_W0; S_G0; ID_M0; ID_W0; ID_G0; UD_M0; UD_W0; UD_G0; A_M0; A_W0; A_G0];

tmax = 1;   % Start with a low end time for debugging purposes
[t,z] = ode45(fun,[0,tmax],z0);  % You may have to try alternative solvers

S = z(:,1);
I = z(:,2);
plot(t,S,t,I);









