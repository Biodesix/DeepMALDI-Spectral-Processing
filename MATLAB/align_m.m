function [ms] = align_m(m,align_ranges,m_coeffs)
%
% Given alignment coefficients and alignment ranges align_m rescales a
% vector m.
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% Input:
%       m                vector of m/z values
%       align_ranges     N_ranges x 2 array of alignment ranges
%       m_coeffs:        N_ranges x 3 array of m-scaling coefficients
%
% Output:
%       ms:              vector of scaled mzs
%
% Dependencies: none
%
% Note: fully vectorized
%
% usage: to generate a set of transformed m/z values use align_m, e.g.
%       mzs_new = align_m(mzs_old,align_ranges,m_coeffs);

    N_ranges = size(align_ranges,1);
    ms   = zeros(size(m,1),size(m,2));
    for ir = 1:N_ranges
        selection = (m > align_ranges(ir,1)) & (m < align_ranges(ir,2));
        tmp = @(x) m_coeffs(ir,1) + m_coeffs(ir,2)*x + m_coeffs(ir,3)*x.^2;
        ms(selection) = tmp(m(selection));
    end

end

