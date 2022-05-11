function [value] = aG(m,m0,sig_L,sig_R)
% aG generates the values of an asymmetric Gaussian defined by its
% position, m0, and the left and right halfwidths, sig_L and sig_R, at the
% positions in the vector m.
%
% author: HR
% version: 1.0.0 (4/21/2021)
% 
% input:        see description
% output:
%       value   arrray of the same size as m, values
%
% Note: this routine is fully vectorized
    
    m_shift     = m - m0;
    tmp         = m_shift <= 0;
    ml          = m_shift( tmp );
    mr          = m_shift( ~tmp );
    fl          = -1 * log(2)/sig_L^2;
    fr          = -1 * log(2)/sig_R^2;
    val_l       = exp( fl*ml.^2 );
    val_r       = exp( fr*mr.^2 );
    
    value       = zeros(size(m,1),size(m,2));
    value(1:length(ml)) = val_l;
    value( (length(ml)+1):end) = val_r;
end

