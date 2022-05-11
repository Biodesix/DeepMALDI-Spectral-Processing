function [left_max,idx] = find_max_toleft(m0,m,maxima)
% find_max_toleft: from an index of maxima. maxima, and a vector of m/z
% values, m, it finds the closest maximum to the left of a value m0.
%  
% author: HR
% version:  0.1.1 (4/21/2021)
%           1.0.0 (5/3/2021)  changed error to warning and retrun defaults
%
% input: see description
% output:
%   left_max:   the value in m closest to the left (smaller) than m0
%   idx:        the index of left_max in m
%
% Note: there are a couple of error conditions (likely unnecessary), and a
% try-catch block to enable debugging

    if ( isempty(maxima) )
        error (' in find_max_toleft: empty maxima' )
    end
    
    if ( isempty(m) )
        error (' in find_max_toleft: empty m' )
    end
    
    try
        max_mz = m(maxima);
        left_max = max(max_mz( max_mz <= m0 ) );
        [~,idx] = min(abs(m - left_max));
    catch
        idx =1;
        left_max = m(1);
        warning (' in find_max_toleft: somethings weird' )
    end
    
end

