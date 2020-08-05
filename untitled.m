syms y(t) z(t)
eqns = [diff(y,t) == (0.0268/5)*z + 0.0927*y - 0.25*y-0.019*y,diff(z,t) == (0.0268*4/5)*z+0.25*y-0.019*z];
[ySol(t),zSol(t)] = dsolve(eqns)