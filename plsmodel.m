function model= plsmodel(XTrn,YTrn,plsparam,XTst,YTst)
% Input:
% XTrn:  spectral matrix in train set;
% YTrn:  response martix in train set;
% XTst:  spectral matrix in test set;
% YTst:  response matrix in test set;

%  Output:
%   model: the result of model;



%%%%%%%%%%%%%%%%% Extract parameters of PLS model %%%%%%%%%%%%%%%%%%%%%
%        LV: The maximum number of PLS;
LV=min([plsparam.LV,size(XTrn,2)]);
%   F_value: Parameter in F test to determine the number of LVs
F_value=plsparam.F_value;
%  plsGroup: the k-fold crossvalidation;
Group=plsparam.Group;
CV=plsparam.CV;
% cross validation of pls model building;
[RMSECV,PTcv]=crossvalidpls(XTrn,YTrn,LV,1,CV,Group);
% RMSECV=crossvalidpls(XTrn,YTrn,LV,1,'full');
%  find the optimal number of principal components;
LVopt=FindPC(RMSECV,length(YTrn),'Ftest_RMSECV',F_value);
if exist('XTst','var')&&exist('YTst','var')
    [B, options, Result] = pls2reg(XTrn, YTrn, LV, 1, XTst, YTst);
    model.Result.PTst=Result.PTst;
    model.RMSEtst=Result.RMSEtst(LVopt);
    model.Result.RMSEP=Result.RMSEtst;
    model.Result.Etst=Result.Etst;
    sse=sum(Result.Etst(:,LVopt).^2);
    sst=sum((YTst-ones(size(YTst,1),1)*mean(YTst)).^2);
    model.R2tst=1-sse/sst;
elseif exist('XTst','var')
    [B, options, Result] = pls2reg(XTrn, YTrn, LV, 1, XTst);
    model.Result.PTst=Result.PTst;
else
    [B, options, Result] = pls2reg(XTrn, YTrn, LV, 1);
end

model.type='pls';
model.Result.PTrn=Result.PTrn;
model.Result.RMSEC=Result.RMSEtrn;
model.RMSEtrn=Result.RMSEtrn(LVopt);
model.Result.PTcv=PTcv;
model.LVopt=LVopt;
model.Result.B=B;
model.RMSECVopt=RMSECV(LVopt);
model.Result.RMSECV=RMSECV;
model.Result.PVEtrn=Result.PVEtrn;
model.options=options;
model.Result.Etrn=Result.Etrn;
sst=sum((YTrn-ones(size(YTrn,1),1)*mean(YTrn)).^2);
sse=sum(Result.Etrn(:,LVopt).^2);
model.R2trn=1-sse/sst;




