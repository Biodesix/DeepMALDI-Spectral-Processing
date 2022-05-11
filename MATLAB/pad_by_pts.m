function [s_pad, i11,i22] = pad_by_pts(s,w)
% pads by w points by mirroring to the left and the right
%
% author: HR
% version: 1.0.0 (4/20/2021)

%   input: 
%       s           data to pad
%       w           number of data points to mirror
%
%   output:
%       s_pad       padded signal (extended by 2*w pts)
%       i11         index to beginning of s in s_pad
%       i22         index to end of s in s_pad
%
    
    if ( w > length(s) )
        error(' cannot pad signal with more points than its length ')
    end
%% extend range by mirroring at the ends by an amount defined by w
    % on the left
    s2pad_l = fliplr(s(1:w));
    %on the right
    w1 = length(s) - w + 1;
    s2pad_r = fliplr(s(w1:end));
    
    s_pad = zeros( length(s) + length(s2pad_l) + length(s2pad_r) , 1);
    s_pad(1:length(s2pad_l)) = s2pad_l;
    s_pad( (length(s2pad_l)+1): (length(s2pad_l)+length(s))) = s;
    s_pad( (length(s2pad_l)+length(s)+1):end) = s2pad_r;

    i11 = length(s2pad_l)+1;
    i22 = length(s2pad_l)+1 + length(s);
    
end

