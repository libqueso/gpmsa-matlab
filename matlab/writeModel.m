function writeModel(params,filename)
% This writes out the data for a GPM/SA model, to be read by the C version
% of the likelihood and sampling code.

n=params.model.n;
m=params.model.m;
p=params.model.p;
q=params.model.q;
pv=params.model.pv;
pu=params.model.pu;

fp=fopen(filename,'w');
fprintf(fp,'%d %d %d %d %d %d\n',n,m,p,q,pv,pu);

fprintf(fp,'%f ',params.data.x);
fprintf(fp,'\n');

fprintf(fp,'%f ',[params.data.z params.data.t]);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.vu);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.w);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.lamOs);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.lamWOs);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.theta);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.betaV);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.lamVz);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.betaU);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.lamUz);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.lamWs);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.LamSim);
fprintf(fp,'\n');

fprintf(fp,'%f ',params.model.SigObs);
fprintf(fp,'\n');


% write prior and MCMC parameters
if n>0 
  varsA={'theta' 'rhoV'  'rhoU'  'lamVz' 'lamUz' 'lamWs' 'lamWOs' 'lamOs'};
  varsB={'theta' 'betaV' 'betaU' 'lamVz' 'lamUz' 'lamWs' 'lamWOs' 'lamOs'};
else
  varsA={'rhoU'  'lamUz'  'lamWs'  'lamWOs'};
  varsB={'betaU' 'lamUz'  'lamWs'  'lamWOs'};
end
for ii=1:length(varsA)
  fprintf(fp,'%f ',params.mcmc.([varsA{ii} 'width'])(1));
  fprintf(fp,'%f ',params.priors.(varsB{ii}).bLower);
  fprintf(fp,'%f ',params.priors.(varsB{ii}).bUpper);
  fprintf(fp,'%f ',params.priors.(varsA{ii}).params(1,:));
  fprintf(fp,'\n');
end


fclose(fp);
  
  

