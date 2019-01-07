function [RMSECV,Ypre] = crossvalidpls(X, Y, h, premethod, valmethod, group, iter)
% Cross-Validation:  selecting the number of PLS latent variables (LV)
%
% Input:
%   X --- m x n matrix,  m spectra, n variables
%   Y --- m x 1,  univariate response
%   h --- the maximum number of PLS LV
%   premethod -  1,'meancenter', 2, 'standard'
%   valmethod - method of crossvalid
%    lout  --- leave-lout-out cross validation.
%    group ---Default: 1, or group >1
%    iter, --- scalar, number of iterations in Monte Carlo Cross-Validation
%
% Output:
%   RMSECV --- 1 x h, cross validation error corresponding to each LV


if nargin < 7
    iter = [];
end
if nargin < 6
    group = 5;  % leave-group-out, used in 'syst123' and 'syst111'
end
if nargin < 5
    valmethod = 'full';  % leave-one-out
end
if nargin < 4
    premethod = 1;
end
if nargin < 3
    h = 10;
end

m = size(X,1);      % size of the data
%% Monte Carlo Cross-Validation
if ~isempty(iter)
    for i = 1 : iter
        indm = randperm(m);
        indt = indm(1:group);         % indices of the test set
        indm(indt) = [];              % indices of the model set
        
        currXTrn= X(indm,:);  currYTrn= Y(indm,:);
        currXTst = X(indt,:);  currYTst = Y(indt,:);
        [~, ~, Result] = pls2reg(currXTrn, currYTrn, h, premethod, currXTst, currYTst);
        tstY = Result.PTst;
        Ypre(indt,:) = tstY;
        sumsqerr(i,:) = sum((repmat(currYTst,1,h) - tstY).^2);
    end
    RMSECV = sqrt( sum(sumsqerr) / (group*iter) );
    %% K-fold Cross-Validation;
else
    gpnum = floor(m/group);
    remainder = m - group*gpnum;
    No = ones(1,group)*gpnum;
    No(1:remainder) = No(1:remainder) + 1;
    
    switch valmethod
        case 'full'
            group = m;
    end
    
    for i  = 1 : group
        switch valmethod
            %% leave-one-out
            case 'full'
                indt = i;
            case 'syst123'
                indt = i : group: m;
            case 'syst111'
                if i == 1
                    indt = 1 : No(i);
                else
                    indt = sum(No(1:i-1)) + 1 : sum(No(1:i));
                end
        end
        indm = 1:m;
        indm(indt) = [];
        currXTrn= X(indm,:);  currYTrn= Y(indm,:);
        currXTst = X(indt,:);  currYTst = Y(indt,:);
        [~, ~, Result] = pls2reg(currXTrn, currYTrn, h, premethod, currXTst, currYTst);
        tstY = Result.PTst;
        Ypre(indt,:) = tstY;
        sumsqerr(i,:) = sum((repmat(currYTst,1,h) - tstY).^2,1);
    end
    RMSECV = sqrt( sum(sumsqerr) / m );
end
