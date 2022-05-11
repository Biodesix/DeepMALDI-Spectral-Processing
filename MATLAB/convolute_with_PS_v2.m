function [conv] = convolute_with_PS_v2(m,s,sig_L,sig_R,peak_range)
% convolute_with_PS performs a convolution of a signal s on points m with a
% peakshape function defined by an asymetric Gaussian with with FWHM sig1
% and sig2 defined by function handles. the signal is mirror padded. It
% uses the trapezoidal rule for the convolution integral. The convolution
% range is a peak_range * peakwidth.
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% input:
%       m           vector of m/z values
%       s           vector of intensities
%       sig_L       function handle for left halfwidth
%       sig_R       function handle for right halfwidth
%       peak_range  parameter value defining convolution range in units of
%                   the peakwidth
%
% output:
%       conv
%
% dependencies:
%               aG
%               pad
%
% Note: this takes a while

    %% extend range by mirroring at the ends by an amount defined by
    [m_pad,s_pad,i11,i22] = pad(m,s,sig_L,sig_R,peak_range); 
    
    %% convolution
    conv = zeros(size(s,1),size(s,2));
    fact = 2*log(2)/sqrt(pi);
    for i = i11:i22
        % get left and right intervals
        m0 = m_pad(i);
        sl = sig_L(m0);
        sr = sig_R(m0);
        pw = peak_range*( sl + sr );
        [~,il] = min(abs(m_pad - (m0 - pw) ));
        [~,ir] = min(abs(m_pad - (m0 + pw) ));
        m_range = m_pad(il:ir+1);
        s_range = s_pad(il:ir+1);
        % get peakshape kernel
        ps = aG(m_range,m0,sl,sr);
        % calculate function ot integrate
        k_range = s_range .* ps;
        % calculate delta xs
        d_ms = m_range(2:end) - m_range(1: (length(k_range)-1));
        % calculate delta f's
        d_fs = k_range(2:end) + k_range(1: (length(k_range)-1));
        % concolution integral
        conv( i-i11 +1) = 0.5 * sum ( d_ms .* d_fs) * fact/pw;
    end

end

