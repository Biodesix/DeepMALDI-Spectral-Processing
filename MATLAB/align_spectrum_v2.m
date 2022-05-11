function [success, m_coeffs, ...
    matched_peaks_s, matched_peaks_r, ...
    pre_alignment_errors, ...
    alignment_errors] = align_spectrum_v2(peaks,amplitudes,...
    ref_peaks, ...
    width_coeffs, ...
    align_ranges, ...
    ALIGN_A_PRCTILE, ...
    ALIGN_MATCH_FRACTION, ...
    ALIGN_MIN_PEAKS)
% assuming possibly multiple alignment ranges align_spectrum calculates the
% coefficients in a quadratic polynomial m_new = a + b*m_old + c*m_old^2
% from minimizing the leat squares error between the peaks of a spectrum
% and a set of predescribed references peaks.
%
% author: HR
% version:  1.0.0 (4/21/2021)
%           1.0.1 (5/3/2021) HR: added success flag
%           2.0.1 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region.
%
% Input: 
%       peaks                   peaks of the spectrum from fit
%       amplitudes              Amplitudes of a spectrum from fit
%       ref_peaks               reference peaks
%       width_coeffs            table of peak shapes 3x3
%       align_ranges            N_ranges x 2 array of alignment ranges
%       ALIGN_A_PRCTILE         percentile of amplitudes per range def = 80
%       ALIGN_MATCH_FRACTION    fraction of peakwidth to define a match to
%                               reference peaks def = 0.5
%       ALIGN_MIN_PEAKS         mimimum number of peaks required per range
%                               for alignment def = 5
%
% Output:
%       success:        flag if alignment is successful = 1
%       m_coeffs:        N_ranges x 3 array of m-scaling coefficients
%       matched_peaks_s: structure array containing the matched sample peaks
%                       per range
%       matched_peaks_r: structure array containing the matched ref peaks
%                       per range
%       pre_alignment_errors: N_ranges x 1 containg the misalingment of
%                       matched peaks as a mean of the absolute deviations
%       alignment_errors:   N_ranges x 1 containg the misalingment of
%                       matched peaks as a mean of the absolute deviations
%                       after alignment
%
% Dependencies: none
%
% Note: throws an error if, in a range, the number of matches is less than
%       ALIGN_MIN_PEAKS
%       may give a warning in the solution of the linear equations; it is
%       save to ignore.
%       Currently it does not check whether alignment points are uniformly
%       distributed in a range which can lead to an accumulation of
%       alignment points around large peak clusters. This should be fixed,
%       e.g. by only picking the largest point in an interval of say 100
%       Da.
%
% usage: to generate a set of transformed m/z values use align_m, e.g.
%       mzs_new = align_m(mzs_old,align_ranges,m_coeffs);
        
    success = 1;

    N_ranges = size(align_ranges,1);
    sample_peaks    = peaks;
    sample_SN       = amplitudes;
    pds_peaks       = ref_peaks;
    

    % Peak shape as function of m/z based on Equation 2 of main text
    % Intersection of linear and quadratic regions
    m_int = width_coeffs{2,'Cutoff'};
    % Left HWHM
    sig_L        = @(x) (width_coeffs{2,'c0'} + width_coeffs{2,'c1'}*x + width_coeffs{2,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{2,'a0'} + width_coeffs{2,'a1'}*x) .* (x < m_int);
    % Right HWHM
    sig_R        = @(x) (width_coeffs{3,'c0'} + width_coeffs{3,'c1'}*x + width_coeffs{3,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{3,'a0'} + width_coeffs{3,'a1'}*x) .* (x < m_int);
    
    % FWHM function
    pw          = @(x) sig_L(x) + sig_R(x);
    
    % loop over ranges
    pre_alignment_errors    = zeros(N_ranges,1);
    alignment_errors        = zeros(N_ranges,1);
    m_coeffs                = zeros(N_ranges,3);
    for ir = 1:N_ranges
        
        i_range = ir;

        % get sub sets
        peaks_ref = pds_peaks( (pds_peaks > align_ranges(i_range,1)) & (pds_peaks < align_ranges(i_range,2)) );
        
        % reduce sample peaks by percentile cutoff within range
        selection       = (sample_peaks > align_ranges(i_range,1)) & (sample_peaks < align_ranges(i_range,2));
        sn_inrange      = sample_SN(selection);
        SN_cutoff       = prctile(sn_inrange,ALIGN_A_PRCTILE);
        SN_selection    = sn_inrange > SN_cutoff;
        tmp             = sample_peaks( (sample_peaks > align_ranges(i_range,1)) & (sample_peaks < align_ranges(i_range,2))  );
        peaks_sam       = tmp(SN_selection);

        % make matches
        matched_sample_peaks    = zeros(length(peaks_ref),1);
        alignment_peaks         = zeros(length(peaks_ref),1);
        delta                   = zeros(length(alignment_peaks),1);
        ic = 0;
        for i=1:length(peaks_ref)
            tmp = peaks_ref(i);
            [~,idx] = min(abs(peaks_sam - tmp));
            delta_tmp = abs(peaks_sam(idx) - tmp);
            if (delta_tmp < ALIGN_MATCH_FRACTION*pw(tmp))
                ic                          = ic + 1;
                matched_sample_peaks(ic)    = peaks_sam(idx);
                alignment_peaks(ic)         = tmp;
                delta(ic)                   = peaks_sam(idx) - tmp;
            end
        end

        N_match = ic;
        
        % Fail and return if does not have at least ALIGN_MIN_PEAKS to
        % align to
        if (N_match < ALIGN_MIN_PEAKS) 
            success         = 0;
            m_coeffs        = [];
            matched_peaks_r = [];
            matched_peaks_s = [];
            pre_alignment_errors = [];
            alignment_errors = [];
            
            warning([' did not find sufficient peaks for range ' num2str(ir)])
            return
        end
        alignment_peaks = alignment_peaks(1:N_match);
        matched_sample_peaks = matched_sample_peaks(1:N_match);
        pre_alignment_errors(ir) = mean(abs(delta(1:N_match)));

        matched_peaks_r(ir).peaks = alignment_peaks; 
        matched_peaks_s(ir).peaks = matched_sample_peaks;
        
        % Calculate alignment with a quadratic linear regression
        fit_matrix = zeros(3,3);
        rhs = zeros(3,1);
        alignment_weigts = ones(N_match,1) ./ (alignment_peaks .* alignment_peaks);
        % Explicit quadratic linear regression matrix 
        for i = 1:length(alignment_peaks)
            fit_matrix(1,1) = fit_matrix(1,1) + alignment_weigts(i);
            fit_matrix(1,2) = fit_matrix(1,2) + alignment_weigts(i)*matched_sample_peaks(i);
            fit_matrix(1,3) = fit_matrix(1,3) + alignment_weigts(i)*matched_sample_peaks(i)^2;
            fit_matrix(2,1) = fit_matrix(2,1) + alignment_weigts(i)*matched_sample_peaks(i);
            fit_matrix(2,2) = fit_matrix(2,2) + alignment_weigts(i)*matched_sample_peaks(i)^2;
            fit_matrix(2,3) = fit_matrix(2,3) + alignment_weigts(i)*matched_sample_peaks(i)^3;
            fit_matrix(3,1) = fit_matrix(3,1) + alignment_weigts(i)*matched_sample_peaks(i)^2;
            fit_matrix(3,2) = fit_matrix(3,2) + alignment_weigts(i)*matched_sample_peaks(i)^3;
            fit_matrix(3,3) = fit_matrix(3,3) + alignment_weigts(i)*matched_sample_peaks(i)^4;
            rhs(1)          = rhs(1) + alignment_weigts(i)*alignment_peaks(i);
            rhs(2)          = rhs(2) + alignment_weigts(i)*alignment_peaks(i)*matched_sample_peaks(i);
            rhs(3)          = rhs(3) + alignment_weigts(i)*alignment_peaks(i)*matched_sample_peaks(i)^2;
        end

        % Linear algebra for quadratic linear regression
        coeffs = fit_matrix \ rhs;
        m_coeffs(ir,:) = coeffs;
        new_m = @(x) coeffs(1) + coeffs(2)*x + coeffs(3)*x.^2; %Alignment coefficients
        fitted = new_m(matched_sample_peaks(1:N_match));
        diff2 = (fitted - alignment_peaks);
        alignment_errors(ir) = mean(abs(diff2));
    end
    
end

