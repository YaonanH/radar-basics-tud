sigma = 1;
P_fa = 0.2;
gamma = sigma * qfuncinv(P_fa)

gamma = 1;
P_fa = qfunc(gamma/sigma)

P_D = qfunc(-2)

P_MD = 1 - P_D
