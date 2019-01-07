% mositure of Tecator data is used as an example
load Tecator_moi
% set the parameters of PLS modeling
plsopts.LV=15;
plsopts.Group=5;
plsopts.CV='syst123';
plsopts.F_value=0.25;
modelpls=plsmodel(Xtrn,ytrn,plsopts,Xtst,ytst);
% set the parameters of variable selection
wscopts.plsopts=plsopts;
wscopts.nBoo=500;
wscopts.nIt=50;
wscopts.ratio=0.1;
% baseline removal
[~,p]=size(Xtrn);
lamb=linspace(-1,1,p);
M=[ones(1,p);lamb;lamb.^2];
P=eye(p)-M'*((M*M')\M);
P=(P+P')/2;
Xtrnp=Xtrn*P;
xm=mean(Xtrnp)';
% perform variable selection
modelwsc=WMSCVS(Xtrnp,ytrn,wscopts);
[~,indCvMin]=min(modelwsc.RMSECVit);
varSelMsk=modelwsc.IndBesit(:,indCvMin);
vasSelInd=find(varSelMsk);
% MSC correction
[~,coeftrn]=emsc(Xtrnp(:,vasSelInd),xm(vasSelInd),'slopeOnly');
XtrnLc=Xtrnp./(coeftrn(:,1)*ones(1,p));
% performance on the test set
Xtstp=Xtst*P;
[~,coeftst]=emsc(Xtstp(:,vasSelInd),xm(vasSelInd),'slopeOnly');
XtstLc=Xtstp./(coeftst(:,1)*ones(1,p));
modelplsLc=plsmodel(XtrnLc,ytrn,plsopts,XtstLc,ytst);

% figures
% Result of variable selection
figure 
stackedbar = @(x, A) bar(x, A,1,'k');
prettyline = @(x, y) plot(x, y);
[AX,H1,H2]=plotyy(wn,varSelMsk,wn,Xtrnp,stackedbar,prettyline);
set(AX,'xlim',[850 1050])
set(AX(1),'ytick',[])
xlabel('Wavelength(nm)')
ylabel('Absorbance')
title('Result of variable selection')
% Process of variable selection
figure
subplot(211)
plot(modelwsc.retIt(:,2));
xlim([1 50])
ylim([0 120])
xlabel('Iterations')
ylabel('Number of sampled variables')
subplot(212)
plot(modelwsc.RMSECVit);
hold on
plot(repmat(indCvMin,1,10),linspace(0,4,10),'k*','LineWidth',0.5);
xlabel('Iterations')
ylabel('RMSECV')
% RMSEP vs LVs
figure
plot(modelpls.Result.RMSEP,'b*-')
hold on
plot(modelplsLc.Result.RMSEP,'r*-')
legend('RAW','WMSCVS')
xlabel('Number of LVs')
ylabel('RMSEP')
title('Moisture')