function [h,cv]=gXval(pout,pvec,varargin)
% function handles=gXval(pout,pvec,nreal,mode)
% leave-one-out cross-validation of eta emulator, using the pvals given in
%    the pout structure (i.e., no re-estimation of params)
%
% cv contains information for computing relative cross-validation
%    errors (computed in 'PCplot')
%    pred - CV predictions
%     cov - CV prediction covariance matrix
%
% mode is a selection of:
%   'PCplot' - (default) plot of each PC response.
%   'PCplotOrder' - plot of each PC response, in canonical order boxplot
%   'residErr' - prediction accuracy of multivariate response
%              (don't try this if it's highly multivariate)
%   'residSummary' - residual summary from each HO
%         (integrating over all multivariate responses, and all pvals)
% numSamp = number of samples to draw on xVal
% % nreal = number of realizations to draw of each point (default 1)
% figNum = figure number to use (default varies by plot)
% standardized = in native space, whether results are standardized or
%                on the original scale (default is 1 = standardized scale)
% labels - sometimes works

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: James R. Gattiker, Los Alamos National Laboratory
%
% This file was distributed as part of the GPM/SA software package
% Los Alamos Computer Code release LA-CC-06-079, C-06,114
%
% Copyright 2008.  Los Alamos National Security, LLC. This material 
% was produced under U.S. Government contract DE-AC52-06NA25396 for 
% Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
% National Security, LLC for the U.S. Department of Energy. The U.S. 
% Government has rights to use, reproduce, and distribute this software.  
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY 
% WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
% THIS SOFTWARE.  If software is modified to produce derivative works, 
% such modified software should be clearly marked, so as not to confuse 
% it with the version available from LANL.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simData=pout.simData;

% set up defaults
numSamp=size(simData.yStd,2);
mode='PCplot'; nreal=1; figNum=false; standardized=1; labels=[];
parseAssignVarargs({'mode','figNum','numSamp','standardized','labels'});  
cv.pred=[]; cv.cov=[];

pvals=pout.pvals(pvec);
npv=length(pvec);
pu=size(simData.Ksim,2);

p=size(simData.Ksim,1);

m=size(simData.yStd,2);
%msamp=gSample(m,numSamp);
%msamp=ilinspace(1,m,numSamp);
msamp=floor(linspace(1,m,numSamp));

% do predictions for all display modes
simLoo=simData;
fprintf('Predicting. '); counter('stime',1,numSamp,6,10);
for ii=1:numSamp
  counter(ii);
  simLoo.yStd=simData.yStd(:,setxor(1:m,msamp(ii)));
  simLoo.x=simData.x(setxor(1:m,msamp(ii)),:);
  pLoo=setupModel([],simLoo,[],'verbose',0);
    %pred(ii)=gPred(simData.x(msamp(ii),:),pvals,pout.model,pout.data,'etamod');
  pred(ii)=gPred(simData.x(msamp(ii),:),pvals,pLoo.model,pLoo.data,'etamod');
  wdat(ii,:)=pout.data.w(msamp(ii),:);
end
counter('end');;

% choose plot type
switch(mode)
case 'PCplot'
  if figNum; figure(figNum); else figure(61); end; clf
  cv.pred=zeros(numSamp,pu); cv.cov=zeros(numSamp,pu,pu);
  for ii=1:numSamp
    Myhat=pred(ii).Myhat; cv.pred(ii,:)=mean(Myhat);
    for jj=1:npv
      cv.cov(ii,:,:)=squeeze(cv.cov(ii,:,:))+pred(ii).Syhat{jj};
      delta=Myhat(jj,:)-cv.pred(ii,:);
      cv.cov(ii,:,:)=squeeze(cv.cov(ii,:,:))+delta'*delta;
    end
    cv.cov(ii,:,:)=squeeze(cv.cov(ii,:,:))./npv;
  end
  isize=ceil(sqrt(pu)); jsize=ceil(pu/isize);
  pcpred=zeros(npv,numSamp);
  for pc=1:pu
    for ii=1:numSamp
      pcpred(:,ii)=mean(squeeze(pred(ii).w(:,pc)),2);
    end
    h(pc)=gPackSubplot(isize,jsize,pc,0,0.1); hold on;
    plot(wdat(:,pc),mean(pcpred),'.');
    for ii=1:numSamp
      plot(wdat(ii,pc)*[1 1],gQuantile(pcpred(:,ii),[0.1 0.9]),...
           'linewidth',1.5);
      %plot(wdat(ii,pc)*[1 1],gQuantile(pcpred(:,ii),[0 1]),...
      %     'linewidth',1);
    end
    a=axis; a=min(a(1),a(3))*[1 0 1 0]+max(a(2),a(4))*[0 1 0 1]; axis(a);
    line(a([1 2]),a([3 4]));
    text(a([1 2])*[0.9 0.1]',a([3 4])*[0.1 0.9]',['PC' num2str(pc)]);
    drawnow;
  end
  
case 'PCplotOrder'
  if figNum; figure(figNum); else figure(62); end; clf
  pcpred=zeros(npv,m);
  for pc=1:pu
    for ii=1:numSamp
      pcpred(:,ii)=mean(squeeze(pred(ii).w(:,pc)),2);
    end
    gPackSubplot(pu,1,pc,1,0);
    gBoxPlot(pcpred);
    drawnow
  end
  
case 'residErr'
  if figNum; figure(figNum); else figure(63); end; clf
  if isscalar(simData.orig.ysd)
    ysd=simData.orig.ysd;
  else
    ysd=repmat(simData.orig.ysd,1,2);
  end
  ymean=repmat(simData.orig.ymean,1,2);
  if standardized; yorig=simData.yStd;
              else yorig=simData.orig.y;
  end
  mipred=zeros(numSamp,p,2);
  mripred=zeros(numSamp,p,2);
  for ii=1:numSamp
    mv=zeros(npv,pu);
    for jj=1:length(pred(ii).Syhat); mv(jj,:)=diag(pred(ii).Syhat{jj}); end
    micdf=(gGMICDF(pred(ii).Myhat',mv',[0.1 0.9])' * simData.Ksim')';
    mricdf=(gGMICDF(pred(ii).Myhat',mv' + pred(ii).residVar',[0.1 0.9])' * simData.Ksim')';
    if ~standardized
      micdf=(micdf.*ysd)+ymean;
      mricdf=(mricdf.*ysd)+ymean;
    end
    mipred(ii,:,:)=micdf;
    mripred(ii,:,:)=mricdf;
  end
  % now plot by response variable
  isize=ceil(sqrt(p)); jsize=ceil(p/isize);
  for ii=1:p
    h(ii)=gPackSubplot(isize,jsize,ii,0,0.1); hold on;
    %fprintf('Plot %d\n',ii);
    for jj=1:numSamp
      mi=squeeze(mipred(jj,ii,:));
      mri=squeeze(mripred(jj,ii,:));
      %fprintf('  '); fprintf('%10.6f ',mi);
      %fprintf(' | '); fprintf('%10.6f ',mri); fprintf('\n')
      plot(yorig(ii,msamp(jj)),mean(pred(jj).Myhat)*simData.Ksim(ii,:)','k+');
      plot(yorig(ii,msamp(jj))*[1 1],mri,'r','linewidth',1);
      plot(yorig(ii,msamp(jj))*[1 1],mi,'linewidth',1.5);
    end
    a=axis; a=min(a(1),a(3))*[1 0 1 0]+max(a(2),a(4))*[0 1 0 1]; axis(a);
    line(a([1 2]),a([3 4]));
    if ~isempty(labels); label=labels{ii}; 
    else label=['var ' num2str(ii)];
    end
    text(a([1 2])*[0.9 0.1]',a([3 4])*[0.1 0.9]',label);
    drawnow;
  end
  

case 'residSummary'  
  if figNum; figure(figNum); else figure(64); end; clf
  if isscalar(simData.orig.ysd)
    ysd=simData.orig.ysd;
  else
    ysd=repmat(simData.orig.ysd,1,npv);
  end
  ymean=repmat(simData.orig.ymean,1,npv);
  if standardized; yorig=simData.yStd;
              else yorig=simData.orig.y;
  end
  for ii=1:numSamp
    y=(pred(ii).w*simData.Ksim')';
    if ~standardized
      y=(y.*ysd)+ymean;
    end
    tres=(y-repmat(yorig(:,msamp(ii)),1,npv));
    yres(:,ii)=tres(:);
  end
  gBoxPlot(yres);
  tstr='residuals'; if standardized; tstr=['standardized ' tstr]; end
  title(tstr);
end



