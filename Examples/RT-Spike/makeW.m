% for MOA analysis
% make w structure with initial values for MOA MCMC

load pout

w.model.theta=pout.model.theta;
w.model.betaU=pout.model.betaU;
w.model.betaV=pout.model.betaV;
w.model.lamVz=pout.model.lamVz;
w.model.lamUz=pout.model.lamUz;
w.model.lamWs=pout.model.lamWs;
w.model.lamWOs=pout.model.lamWOs;
w.model.lamOs=pout.model.lamOs;

w.mcmc.thetawidth=pout.mcmc.thetawidth;
w.mcmc.rhoUwidth=pout.mcmc.rhoUwidth;
w.mcmc.rhoVwidth=pout.mcmc.rhoVwidth;
w.mcmc.lamVzwidth=pout.mcmc.lamVzwidth;
w.mcmc.lamUzwidth=pout.mcmc.lamUzwidth;
w.mcmc.lamWswidth=pout.mcmc.lamWswidth;
w.mcmc.lamWOswidth=pout.mcmc.lamWOswidth;
w.mcmc.lamOswidth=pout.mcmc.lamOswidth;

save wRT w
