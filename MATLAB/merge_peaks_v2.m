function [merged] = merge_peaks_v2(peaks,width_coeffs,COMBINE_FRACTION)
%
% Merges peaks if they are closer together than COMBINE_FRACTION*peakwidth.
%
% author: HR, MK
% version: 1.0.0 (4/21/2021)
%
% Input: 
%       peaks               peaks 
%       width_coeffs        3x3 table of peak shape parameters
%       COMBINE_FRACTION    parameter
%
% Output:
%       merged   vector of merged peaks
    
    if (length(peaks) == 1)
        merged = peaks;
        return
    end
    
    % set up peakshape as a function of m/z from Equation 2 of main text
    % Linear-to-quadratic cutoff (in paper as m_{int})
    m_int = width_coeffs{2,'Cutoff'};

    % Left and right HWHM coefficients from Equation 2 of main text
    sig_L        = @(x) (width_coeffs{2,'c0'} + width_coeffs{2,'c1'}*x + width_coeffs{2,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{2,'a0'} + width_coeffs{2,'a1'}*x) .* (x < m_int);
    sig_R        = @(x) (width_coeffs{3,'c0'} + width_coeffs{3,'c1'}*x + width_coeffs{3,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{3,'a0'} + width_coeffs{3,'a1'}*x) .* (x < m_int);
    % FWHM function
    pw          = @(x) sig_L(x) + sig_R(x) ;

    tmpout = peaks;
    reduction = 1;
    
    % Merge peaks 
    while (reduction > 0)
        tmp = zeros(length(tmpout),1);
        NP = 1;
        tmp(1) = tmpout(1);
        for i=2:length(tmpout)
            delta = tmpout(i) -tmpout(i-1);
            midp = 0.5*(tmpout(i) + tmpout(i-1));
            % If close, merge to mid point
            if ( delta < COMBINE_FRACTION*pw(midp) )
                tmp(NP) = midp;
            % Otherwise keep separate
            else
                NP = NP + 1;
                tmp(NP) = tmpout(i);
            end
        end
        reduction = length(tmpout) - NP;
        tmpout = tmp(1:NP);
    end
    
    merged = tmpout;
    
end

