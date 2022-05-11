function [noisep] = noise2percentiles(noise,w,percentile)
% noise2percentiles retruns a local noise estimator noisep from the 
% detailed noise using percentiles percentile over
% a width w. Boundaries are extended by mirror b.c.s.
%
% author: HR
% version: 1.0.0 (4/20/2021)
%
% Input/output: see description
%
% dependencies:
%   pad_by_pts  (for mirror boundary conditions)
%
    if ( w > length(noise) )
        error( ' too large window in noise2percentiles ')
    end
    
    [np,i1,~] = pad_by_pts(noise,w+2);
    
    noisep = zeros( size(noise,1),size(noise,2) );
    for i = 1:length(noise)
        ip = i + i1;
        data = np( (ip-w):(ip+w) );
        tmp = prctile(data,percentile);
        noisep(i) = tmp;      
    end
    
end

