function L=gThetaPrior(model,priors)

 mu=priors.mean; chCov=priors.chCov;
 ind=priors.ind;
 theta=(priors.xrange.*model(ind)+priors.xmin)';
 gamma=theta(3); theta(3)=-log(gamma);

 p1=(chCov')\(theta-mu);
 L=-0.5*(p1'*p1)+theta(3);

 ind=setxor(1:length(model),ind); 
 L = L+feval('gLogNormalPrior',model(ind), ...
             priors.default(ind,:));

end
