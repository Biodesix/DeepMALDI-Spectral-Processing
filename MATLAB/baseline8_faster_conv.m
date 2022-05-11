function [z,conv] = baseline8_faster_conv(y, lambda1, lambda2, p, iter)
%BASELINE Estimate baseline with asymmetric least squares (Eilers)
%
% author: HR
% version: 4/20/2021
% Generalized to asymmetric rigidity:
% lambda1 (higher) is rigidity for negative second derivative (z'' < 0)
% lambda2 (lower)  is rigidity for positive second derivative (z'' >= 0)
%
% Input:
%   y           data
%   lambda1     Eilers parameter
%   lambda2     Eilers parameter
%   p           Eilers parameter
%   iter        max number of iterations
% Output:
%   z           baseline estimator
%   conv        L^2 change in background


m = length(y);
D = diff(speye(m), 2);
w = ones(m, 1);
ww = 0.5*(lambda1 + lambda2) * ones(m-2, 1);

conv = zeros(iter-1,1);

% To improve convergence:
% First, do 10 iterations with lambda2a gradually changing (in log scale) 
% from lambda1 to lambda2

iter1 = min(iter, 10);
ic = 0;
for it = 1:iter1
    lambda2a = lambda1*(lambda2/lambda1)^((1.0*it)/iter1);
    W = spdiags(w, 0, m, m);
    WW = spdiags(ww, 0, m-2, m-2);

    
    if( it > 1)
        zpre = z;
    end
    z = (W + D' * WW * D) \ (w .* y);
    if( it > 1)
        ic = ic + 1;
        conv(ic) = norm((zpre-z) )/m;
    end
    w = p * (y > z) + (1 - p) * (y < z);
    ww = lambda2a*(D*z >= 0) + lambda1*(D*z < 0);
end

% Now we already have initialization for w and ww 
% left over from the loop above.

% Do remaining iterations with actual lambda1, lambda2.

iter2 = iter - iter1;

for it = 1:iter2
    W = spdiags(w, 0, m, m);
    WW = spdiags(ww, 0, m-2, m-2);

    
        zpre = z;
    z = (W + D' * WW * D) \ (w .* y);
        ic = ic + 1;
        conv(ic) = norm((zpre-z) )/m;
    w = p * (y > z) + (1 - p) * (y < z);
    ww = lambda2*(D*z >= 0) + lambda1*(D*z < 0);
end

end


