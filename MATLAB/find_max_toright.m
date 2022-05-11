function [right_max,idx] = find_max_toright(m0,m,maxima)
% find_max_toright: from an index of maxima. maxima, and a vector of m/z
% values, m, it finds the closest maximum to the right of a value m0.
%  
% author: HR
% version: 0.1.1 (4/21/2021)
%
% input: see description
% output:
%   right_max:  the value in m closest to the right(greater) than m0
%   idx:        the index of right_max in m
%

    max_mz = m(maxima);
    if (m0 > max_mz(end) )
        right_max = m(end);
        idx = length(m);
        return;
    end
    
    right_max = min(max_mz( max_mz >= m0 ) );
    
    if (size(right_max,1) ~= 1 | size(right_max,2) ~= 1)
        display([' something went wrong in find_max_toright m0 ' num2str(m0)])
        display([' max_mz(1) ' num2str(max_mz(1)) ' max_mz(N) ' num2str(max_mz(end)) ])
        error(' there is no maximum > input argument 1')
    end
    
    [~,idx] = min(abs(m - right_max));
    
end

