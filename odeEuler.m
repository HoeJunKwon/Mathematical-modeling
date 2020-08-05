function dy = odeEuler(T)
% Parameters of the model
mi = 1.0 / 69.54; b = 2.1 * mi; beta = 1.6; 
etaC = 0.015; etaA = 1.3; fi = 1.0; ro = 0.1; 
alfa = 0.33; omega = 0.09; d = 1.0;

% Parameters of the Euler method
test = -1; deltaError = 0.001; M = 100; t = linspace(0,T,M+1); h = T / M;
S = zeros(1,M+1); I = zeros(1,M+1);
C = zeros(1,M+1); A = zeros(1,M+1);
% Initial conditions of the model
S(1) = 0.6; I(1) = 0.2; C(1) = 0.1; A(1) = 0.1;
% Iterations of the method 
while(test < 0)
        oldS = S; oldI = I; oldC = C; oldA = A;
        for i = 1:M
           % Differential equations of the model
           aux1 = beta * (I(i) + etaC * C(i) + etaA * A(i)) * S(i);
           aux2 = d * A(i);
           auxS = b * (1 - S(i)) - aux1 + aux2 * S(i);
           auxI = aux1 - (ro + fi + b - aux2) * I(i) + alfa * A(i) + omega * C(i); 
           auxC = fi * I(i) - (omega + b - aux2) * C(i);
           auxA = ro * I(i) - (alfa + b + d - aux2) * A(i);
           
           % Euler new approximation
           S(i+1) = S(i) + h * auxS;
           I(i+1) = I(i) + h * auxI;
           C(i+1) = C(i) + h * auxC;
           A(i+1) = A(i) + h * auxA;
        end
        
        % Absolute error for convergence
        temp1 = deltaError * sum(abs(S)) - sum(abs(oldS - S)); 
        temp2 = deltaError * sum(abs(I)) - sum(abs(oldI - I)); 
        temp3 = deltaError * sum(abs(C)) - sum(abs(oldC - C)); 
        temp4 = deltaError * sum(abs(A)) - sum(abs(oldA - A)); 
        test = min(temp1,min(temp2,min(temp3,temp4)));
end     
dy(1,:)  = t; dy(2,:)  = S; dy(3,:)  = I;
dy(4,:)  = C; dy(5,:)  = A;