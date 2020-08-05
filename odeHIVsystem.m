function dy = odeHIVsystem(t,y)
% Parameters of the model
mi = 1.0 / 69.54; b = 2.1 * mi; beta = 1.6; etaC = 0.015; etaA = 1.3; fi = 1.0; ro = 0.1;
alfa = 0.33; omega = 0.09; d = 1.0;
% Differential equations of the model
dy = zeros(4,1);
aux1 = beta * (y(2) + etaC * y(3) + etaA * y(4)) * y(1);
aux2 = d * y(4);
dy(1) = b * (1 - y(1)) - aux1 + aux2 * y(1);
dy(2) = aux1 - (ro + fi + b - aux2) * y(2) + alfa * y(4) + omega * y(3); dy(3) = fi * y(2) - (omega + b - aux2) * y(3);
dy(4) = ro * y(2) - (alfa + b + d - aux2) * y(4);