function [B, options,Result] = pls2reg(X, Y, PC, premethod, XTst, YTst)
% The code are written for sigle response variable ( Y is a column vector
% with only one component).

if nargin < 4
    premethod = 1;
end
if nargin < 3
    PC = 10;
end
%% Meancenter
ntrn= size(X,1);
if premethod == 1
    Xuse  = X - ones(ntrn,1)*mean(X);
    Yuse  = Y - ones(ntrn,1)*mean(Y);
    %% Standard
elseif premethod == 2
    Xuse  = (X - ones(ntrn,1)*mean(X))./ (ones(ntrn,1)*std(X));
    Yuse  = (Y - ones(ntrn,1)*mean(Y))./ (ones(ntrn,1)*std(Y));
end
%% Main function on PLS regression coefficient B;
PVE=zeros(2,PC);
for i = 1 : PC
    stopcond = 0;
    u = Yuse(:,1);
    count=1;
    while ~stopcond
        w = Xuse'*u;
        w = w / sqrt(w'*w);
        t = Xuse*w;
        q= (Yuse'*t)/(t'*t);
        q=q/sqrt(q'*q);
        unew = (Yuse*q)/(q'*q);
        diff = (u-unew)'*(u-unew);
        if diff < eps || count>500;
            stopcond = 1;
        end
        u = unew;
        count=count+1;
    end
    
    p = (Xuse'*t) / (t'*t);
    b(i)=(t'*u)/(t'*t);
    Xuse = Xuse - t*p';
    Yuse = Yuse - b(i)*t*q';
    W(:,i) = w;
    T(:,i) = t;
    P(:,i) = p;
    Q(:,i) = q;
    U(:,i) = u;
    R(:,i)=b(i)*t;
    B(:,i) = W(:, 1:i) * inv(P(:,1:i)'*W(:, 1:i)+10^(-6)*eye(i)) *diag(b(1:i))* Q(:,1:i)';
    
end
%% Exact result:
options.PC = PC;
options.premethod = premethod;
options.weight   = W;
options.xscore   = T;
options.yscore  =U;
options.xloading = P;
options.yloading = Q;
options.Regcoef  = B;
options.R=R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLS on the train set;
PTrn =zeros(ntrn,PC);
%% Meancenter
if premethod == 1
    XTrnm=X-ones(ntrn,1)*mean(X);
    PTrn=XTrnm*B+ones(ntrn,PC)*mean(Y);
    %% Standard
elseif premethod == 2
    XTrnms=(X-ones(ntrn,1)*mean(X))./(ones(ntrn,1)*std(X));
    PTrn=XTrnms*B*std(Y)+ones(ntrn,PC)*mean(Y);
end

Etrn=repmat(Y,1,PC) - PTrn;
RMSEtrn = sqrt( sum( Etrn.^2,1 ) /ntrn);
PVEtrn= PVEcal(T,R,P,Q,X,Y,premethod);
Result.Etrn=Etrn;
Result.PTrn=PTrn;
Result.RMSEtrn=RMSEtrn;
Result.PVEtrn=PVEtrn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLS on the test set;
if exist('XTst','var')
    ntst = size(XTst,1);
    PTst = zeros(ntst,PC);
    
    %% Meancenter
    if premethod == 1
        XTstm=XTst-ones(ntst,1)*mean(X);
        PTst= XTstm*B+ones(ntst,PC)*mean(Y);
        for i=1:PC;
            Ttst(:,i)=XTstm*W(:,i);
            XTstm=XTstm-Ttst(:,i)*P(:,i)';
            Utst(:,i)=b(i)*Ttst(:,i);
        end
        %% Standard
    elseif premethod == 2
        XTstms=(XTst-ones(ntst,1)*mean(X))./(ones(ntst,1)*std(X));
        PTst=XTstms*B*std(Y)+ones(ntst,PC)*mean(Y);
        for i=1:PC;
            Ttst(:,i)=XTstms*W(:,i);
            XTstms=XTstms-Ttst(:,i)*P(:,i)';
            Utst(:,i)=b(i)*Ttst(:,i);
        end
    end
    
    if exist('YTst','var')
        Etst=repmat(YTst,1,PC) - PTst;
        RMSEtst = sqrt( sum( Etst.^2 ,1) /ntst);
        Result.RMSEtst =RMSEtst;
        Result.Etst=Etst;
    end
    Result.PTst=PTst;
    options.Ttst   = Ttst;
    options.Utst   = Utst;
end








