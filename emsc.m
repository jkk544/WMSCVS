function [emscData,emscParam,Ref,Xmod] = emsc(data,Ref,type,MSCWeights, extraModelParam)

%  Performs Extended Multiplicative Signal Correction (EMSC)
%  Xorg = [1 Xref Xref^2 lambda lambda^2 xknown1 xknown2 ...]*b + e
%
% Input Variables:
%  data:        Data matrix to be transformed, m x n, m spectra, n channels
%  offset:      1 for including additive model term. 0 otherwise.
%  slope:       1 for including multiplicative model term. 0 otherwise.
%               The multiplicative term is chosen as the mean spectrum.
%  slopeSquared: 1 for including squared an additive model term
%                (the square of the mean spectrum). 0 otherwise.
%  xValues:     1 for including the x-variable values (ex. Wavelengths in spectra)
%               in the model. 0 otherwise. Additive.
%  xValuesSquared: 1 for including the square of the x-variable values.
%                0 otherwise. Additive
%  extraModelParam: Extra model Parameter matrix to include in the model.
%               The matrix must have the same number of columns as x-variables.
%               [] or 0 if not used.
%
% Output Variables:
%   emscData:    EMSC transformed data matrix
%   emscParam:   Correction parameters for the different model terms:
%                [offsetPar, xValuesPar,xValuesSqPar,slopeSqPar, extraModPar, slopePar]

% Default:  offset = 1; slope = 1; slopeSquared = 0; xValues = 1;xValuesSquared = 1;
%  extraModelParam = 0;
%  extraModelParam = data(2,:) - data(47,:);

X = data;
[M,N] = size(X);
if nargin < 5
    extraModelParam = 0;
end

if nargin <4
    MSCWeights = ones(1,N);
end


if isempty(Ref)
    Ref = mean(X)';
end

slope = 1;
switch type
    case 'regular'
        offset = 1;
        slopeSquared = 0;
        xValues = 1;
        xValuesSquared = 1;
    case 'slopeOnly'
        offset = 0;
        slopeSquared = 0;
        xValues = 0;
        xValuesSquared = 0;
end



CondNumber = 1e19;


[M2,N2] = size(extraModelParam);
Xmod = [];


%Offset
if offset == 1
    Xmod = [ones(N,1)];
end % if

%xValues
if xValues == 1
    lambda = linspace(-1,1,N);
    Xmod = [Xmod lambda'];
end % if

%xValuesSquared
if xValuesSquared == 1
    lambdaSq =lambda.^2;
    Xmod = [Xmod lambdaSq'];
end % if

%Slope Squared
if slopeSquared == 1
    mRef = mean(Ref);
    RefSq = (Ref - mRef).^2;
    Xmod = [Xmod RefSq];
end % if

% Extra Model Parameters
if extraModelParam == 0
    extraModelParam = [];
end % if

flagVal = ~isempty(extraModelParam);

if flagVal == 1
    Xmod = [Xmod extraModelParam'];
end % if

%Slope
if slope == 1
    Xmod = [Xmod Ref];
end % if

diagW = diag(MSCWeights);
XX = Xmod'*diagW*Xmod;

% Ridging if CondNumber is specified. Calculetion of inverse XX.

if CondNumber>0
    [u,S,v]=svd(XX);
    s=diag(S)';
    sCorrected=max(s, (s(1)/CondNumber));
    XXCorrected=u*(diag(sCorrected))*v';
    InvXX = inv(XXCorrected);
else
    InvXX=inv(XX);
end % if

MSCParam = InvXX*Xmod'*diagW*X';
MSCParam = MSCParam';

%estimating Xnew

Xnew = X;
j = 1;

%Offset
if offset == 1
    ai = MSCParam(:,j);
    Xnew = Xnew-(ai*ones(1,N));
    j = j+1;
end % if

% xValues
if xValues == 1
    ei = MSCParam(:,j);
    Xnew = Xnew-(ei*lambda);
    j = j+1;
end % if

% xValues Squared
if xValuesSquared == 1
    fi = MSCParam(:,j);
    Xnew = Xnew-(fi*lambdaSq);
    j = j+1;
end % if

% Slope Squared
if slopeSquared == 1
    b2i = MSCParam(:,j);
    Xnew = Xnew - (b2i*RefSq');
    j = j+1;
end % if

%Extra Model Parameters
if flagVal == 1
    H = MSCParam(:,j:j+(M2-1));
    j = j+M2;
end % if

%Slope
if slope == 1
    bi = MSCParam(:,j);
    Xnew = Xnew./(bi*ones(1,N));
    j = j+1;
end % if

emscData = Xnew;
emscParam = MSCParam;
