function f = FindPC(Press, ntrn, method, alpha)
% deiermine the optimal latent variables (principal compononts) based on
% Press value. Most of cases, F-test is used. 
%
% When the given input Press is RMSECV, 
% using Ftest_RMSECV (or equivalently, Ftest with Press = RMSECV.^2*ntrn )
%
% Note: The freedom and siginificant alpha are not consistent in the papers
%           freedom = ntrn or ntrn-1 and  alpha or alpha/2 
%       In practice, the alpha is given by user. default is 0.25

if nargin < 4
    alpha = 0.25;
end
dper = 1 - alpha;

switch method
     case 'Ftest'  
        % Haaland and Thomas, PLS for Spectral Analyses, 1998, (Press)
        % here Press=RMSECV.^2*ntrn;
        Press=Press.^2*ntrn;
        [x, k] = min(Press);
        F = Press(1:k) / Press(k);
        dper = 1 - alpha/2;  % here, alpha/2 
        Fminus = F - ones(size(F))*finv(dper, ntrn, ntrn); % freedom = ntrn 
        k1 = find(Fminus<0);
        f = k1(1);
        
    case 'Ftest_RMSECV'
        % here Press = RMSECV
        [x,k] = min(Press);
        f = 1;
        while Press(f) > Press(k) * sqrt(finv(dper,ntrn-1,ntrn-1));
            f = f+1;
        end
%         f=k;
    case 'Ftest_RMSECV1'
        % based on RMSECV, freedom = ntrn - 1;  
        % Selection of individual variables versus intervals of variables in PLSR
        [x,k] = min(Press);
        nonsignincr = finv(dper,ntrn-1,ntrn-1);
        f = 1;
        while Press(f)-Press(k) > nonsignincr
            f = f+1;
        end
        
    case 'localmin'
        % find the first local minimum PRESS
        f = length(Press);
        for j = 2 : length(Press)-1
            if (Press(j) < Press(j-1)) && (Press(j) < Press(j+1))
                f = j;
                break;
            end
        end
        
    case 'minsingular'
        % find the minimum PRESS, and no 2% surpass
        [x,k] = min(Press);
        % nonsignincr = Press(k)*0.02;
        nonsignincr = ( max(Press) - min(Press) ) * 0.05;
        f = 1;
        while Press(f)-Press(k) > nonsignincr
            f = f+1;
        end
       
   
end