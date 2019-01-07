function model=WMSCVS(Xtrnp,ytrn,opts)
%+++ WMSCVS: Weighted Multiplicative Scatter Correction with Variable Selectionn
%+++ Xtrnp: The data matrix of size n x p, spectra matrix after baseline removal
%+++ ytrn: The reponse vector of size n x 1
%+++ opts: The parameter settings of WMSCVS
%+++ Author: spotwu@163.com

% Initial settings
[~,p]=size(Xtrnp);
xm=mean(Xtrnp)';
nBoo=opts.nBoo;
nIt=opts.nIt;
ratio=opts.ratio;
plsopts=opts.plsopts;
nBes=floor(nBoo*ratio);
W=zeros(p,nIt+1);
W(:,1)=ones(p,1)/p;
RMSECVit=zeros(nIt,1);
RMSECit=zeros(nIt,1);
IndBesit=zeros(p,nIt);
IndBesGrouit=cell(nIt,1);
RMSECVGrouit=zeros(nBes,nIt);
hisBesit=zeros(p,nIt);
lvBesit=zeros(nIt,1);
retIt=zeros(nIt,2);
% Caculate EDF parameters
a=(p/2)^(1/(nIt-1));
b=log(p/2)/(nIt-1);
for i=1:nIt
    fprintf('Processing iteration %d\n',i)
    selMat=zeros(p,nBoo);
    rmsecv=zeros(nBoo,1);
    rmsec=zeros(nBoo,1);
    lv=zeros(nBoo,1);
    w=W(:,i);
    rati=a*exp(-b*i);
    retIt(i,1)=rati;
    nRet=ceil(p*rati);
    retIt(i,2)=nRet;
    for j=1:nBoo
        indsel=randsample(p,nRet,true,w);
        indsel=unique(indsel);
        selMat(indsel,j)=1;
        % MSC correction
        [~,coeftrn]=emsc(Xtrnp(:,indsel),xm(indsel),'slopeOnly');
        XtrnLc=Xtrnp./(coeftrn(:,1)*ones(1,p));
        % PLS model
        modelpls=plsmodel(XtrnLc,ytrn,plsopts);
        lv(j)=modelpls.LVopt;
        rmsecv(j)=modelpls.RMSECVopt;
        rmsec(j)=modelpls.RMSEtrn;
    end
    [~,indcvSort]=sort(rmsecv);
    RMSECVit(i)=rmsecv(indcvSort(1));
    RMSECit(i)=rmsec(indcvSort(1));
    lvBesit(i)=lv(indcvSort(1));
    IndBesit(:,i)=selMat(:,indcvSort(1));
    IndBesGrouit{i}=selMat(:,indcvSort(1:nBes));
    RMSECVGrouit(:,i)=rmsecv(indcvSort(1:nBes));
    hisBesit(:,i)=sum(selMat(:,indcvSort(1:nBes)),2);
    w=hisBesit(:,i);
    w=w/sum(w);
    W(:,i+1)=w;
end
model.W=W;
model.RMSECVit=RMSECVit;
model.RMSECit=RMSECit;
model.IndBesit=IndBesit;
model.hisBesit=hisBesit;
model.lvBesit=lvBesit;
model.retIt=retIt;
model.IndBesGrouit=IndBesGrouit;
model.RMSECVGrouit=RMSECVGrouit;
end