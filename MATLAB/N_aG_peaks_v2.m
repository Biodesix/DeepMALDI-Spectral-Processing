function [values] = N_aG_peaks_v2(b,m,NP,sig_L,sig_R)
%
% N_aG_peaks: calculates the values of NP asymmetric Gaussians with
% predefined positions and amplitudes common from an array usable for
% nlinfit.
%
% author: HR
% version:  1.0.0 (4/21/2021)
%           2.0.0 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region.
%
% input:
%   b       vector of size 2*NP containing positions (even index) and
%           amplitudes (odd index)
%   m       vector of m/z values
%   NP      number of asymmetric Gaussians
%   sig_L   function handle for left peakwidth
%   sig_R   function handle for right peakwidth
%
% dependencies:     aG
%
% note: fully vectorized
    
    % loop over peaks
    values = zeros(size(m,1),size(m,2));
    for ip = 1:NP
        m0 = b(2*(ip-1)+1);
        values = values + b(2*ip)*aG(m,m0,sig_L(m0),sig_R(m0));
    end
end

