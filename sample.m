



syms N_M(t) N_W(t) N_G(t) S_M(t) S_W(t) S_G(t) ID_M(t) ID_W(t) ID_G(t) UD_M(t) UD_W(t) UD_G(t) A_M(t) A_W(t) A_G(t) 
ode = diff(N_M,t) == t*y;
ode = diff(N_W,t) == t*y;
ode = diff(N_G,t) == t*y;
ode = diff(S_M,t) == t*y;
ode = diff(S_W,t) == t*y;
ode = diff(S_G,t) == t*y;
ode = diff(ID_M,t) == t*y;
ode = diff(ID_W,t) == t*y;
ode = diff(ID_G,t) == t*y;
ode = diff(UD_M,t) == t*y;
ode = diff(UD_W,t) == t*y;
ode = diff(UD_G,t) == t*y;
ode = diff(A_M,t) == t*y;
ode = diff(A_W,t) == t*y;
ode = diff(A_G,t) == t*y;




cond = N_M(0) == 0;
cond = N_W(0) == 0;
cond = N_G(0) == 0;
cond = S_M(0) == 0;
cond = S_W(0) == 0;
cond = S_G(0) == 0;
cond = ID_M(0) == 0;
cond = ID_W(0) == 0;
cond = ID_G(0) == 0;
cond = UD_M(0) == 0;
cond = UD_W(0) == 0;
cond = UD_G(0) == 0;
cond = A_M(0) == 0;
cond = A_W(0) == 0;
cond = A_G(0) == 0;

ySol(t) = dsolve(ode,cond)