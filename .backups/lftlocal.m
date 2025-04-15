function pdGcl = lftlocal(pdG,pdK,pv)

% to compute the closed loop model


G = ss(pdG,pv.points);
K = ss(pdK,pv.points);
Gcl = lft(G,K);

pdGcl = ppss(Gcl,pv);


