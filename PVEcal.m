function PVE= PVEcal(T,R,P,Q,X,Y,premethod )
if nargin<7
    premethod=1;
end
%% Meancenter
n= size(X,1);
if premethod==0
    X0=X;
    Y0=Y;
elseif premethod == 1
    X0  = X - ones(n,1)*mean(X);
    Y0  = Y - ones(n,1)*mean(Y);
    %% Standard
elseif premethod == 2
    X0  = (X - ones(n,1)*mean(X))./ (ones(n,1)*std(X));
    Y0  = (Y - ones(n,1)*mean(Y))./ (ones(n,1)*std(Y));
end
trX=trace(X0'*X0);
trY=trace(Y0'*Y0);
%% Compute the percent of variance explained for X and Y
for i=1:size(T,2);
    t=T(:,i);
    r=R(:,i);
    p=P(:,i);
    q=Q(:,i);
    PVEX(1,i)=(t'*t)*(p'*p)/trX*100;
    PVEY(1,i)=(r'*r)*(q'*q)/trY*100;
end
PVE=[PVEX;cumsum(PVEX);PVEY;cumsum(PVEY)];

