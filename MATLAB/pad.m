function [m_pad,s_pad,i11,i22] = pad(m,s,sig_L,sig_R,peak_range)
% mirror pads m/z and intensity by a range defined by peak_range *
% peak_width.
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
%       m_pad       padded m/z values
%       s_pad       padded intensities
%       i11         points to the first point of m,s in m_pad,s_pad
%       i22         points to the last point of m,s in m_pad,s_pad
%

    %% extend range by mirroring at the ends by an amount defined by
    %  peak_range
    % on the left
    w1 = 1.2*peak_range*( sig_L(m(1)) + sig_R(m(1)) );    % 1.2 to have some grace
    [~,i1] = min( abs( m -(m(1) + w1)));
    s2pad_l = fliplr(s(1:i1));
    m2pad_l = zeros(length(s2pad_l),1);
    mi = m(1);
    dm = m(2) - m(1);
    for i = length(s2pad_l):-1:1
        mi = mi - dm;
        m2pad_l(i) = mi;
    end
    %on the right
    w1 = 1.2*peak_range*( sig_L(m(end)) + sig_R(m(end)) );    % 1.2 to have some grace
    [~,i1] = min( abs( m -(m(end) - w1)));
    s2pad_r = fliplr(s(i1:end));
    m2pad_r = zeros(length(s2pad_r),1);
    mi = m(end);
    dm = m(end) - m((length(m) -1));
    for i = 1:length(s2pad_r)
        mi = mi + dm;
        m2pad_r(i) = mi;
    end
    
    s_pad = zeros( length(s) + length(s2pad_l) + length(s2pad_r) , 1);
    m_pad = zeros( length(s) + length(s2pad_l) + length(s2pad_r) , 1);
    s_pad(1:length(s2pad_l)) = s2pad_l;
    s_pad( (length(s2pad_l)+1): (length(s2pad_l)+length(s))) = s;
    s_pad( (length(s2pad_l)+length(s)+1):end) = s2pad_r;
    m_pad(1:length(s2pad_l)) = m2pad_l;
    m_pad( (length(s2pad_l)+1): (length(s2pad_l)+length(s))) = m;
    m_pad( (length(s2pad_l)+length(s)+1):end) = m2pad_r;
    
    i11 = (length(s2pad_l)+1);
    i22 = (length(s2pad_l)+length(s));
    
end

