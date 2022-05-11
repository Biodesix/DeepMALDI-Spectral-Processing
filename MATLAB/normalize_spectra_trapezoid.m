function [sum_int] = normalize_spectra_trapezoid(mzs,its,left_mz,right_mz)
% normalize_spectra takes an incoming spectra m/z (mzs) and intensity (its)
% and normalizes the intensity based on a trapezoidal integration
%
% author: MK
% version: 1.0.0 (6/2/2021)
%
% input:
%       mzs         vector of m/z values
%       its         vector of intensities
%       left_mz     start of the integral in m/z
%       right_mz    end of the integral in m/z
%
% output:
%       sum_int     integrated intensity of its from left_mz to right_mz
%
% dependencies:     none 
%              


%Find the left and right index to integrate over
[~,il] = min(abs(mzs - left_mz ));
[~,ir] = min(abs(mzs - right_mz ));

% Run conditions if il or ir is >max
if il < 1
    il = 1;
end
if ir > length(mzs)
    ir = length(mzs);
end

%Integrate using the trapezoid rule
d_mzs = mzs( (il+1):ir ) - mzs(il:(ir-1));
s_its = its( (il+1):ir ) + its(il:(ir-1));

sum_int = 0.5 * sum(d_mzs .* s_its);


end

