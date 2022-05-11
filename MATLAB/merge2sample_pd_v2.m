function [all_peaks, all_mults] = merge2sample_pd_v2(...
    ref_peaks, ref_mult,...
    check_peaks, check_mult,...
    width_coeffs,COMBINE_FRACTION)
%
% from two sets of definitions ( peaks and multiplicities) generates a
% combined list of definitions. It first finds isolated single peaks and
% matches them. matches are if the 2 peaks are closer than
% COMBINE_FRACTION* peakwidth. It then makes a list of multiple peaks (just
% takes all from both lists) and defines boundaries at the end of the same
% multiplicities. This list is then matched internally, and new list of
% clusters (single or multiple peaks is generated). These are then merged
% if they are too close (at the beginning and end). Finally this cluster list
% converted back to peaks and their multiplicities.
%
% author: HR
% version:  1.0.0 (4/21/2021)
%           2.0.1 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region.
%
% Input: 
%       ref_peaks           peaks from set 1
%       ref_mult            corresponding multiplicities
%       check_peaks         peaks from set 2
%       check_mult          corresponding multiplicities
%       width_coeffs        3x3 table of peak shape parameters
%       COMBINE_FRACTION    fraction of peakwidth for matching def 
%
% Output:
%       all_peaks   the combined peak list
%       all_mults   corresponding multiplicities
%
% Dependencies: 
%           merge_peaks
%           make_clusters
%           merge_clusters
%           make_flag_from_mult
%
% Notes: a check that everything works is to check whether all_peaks is
% sorted. 

    % Setup initial parameters
    N_ref = length(ref_peaks);
    N_check = length(check_peaks);
    minimum_ref_distance = min(diff(ref_peaks,1));

    matched_ref     = zeros(N_ref,1);
    matched_check   = zeros(N_check,1);
    match_index     = zeros(N_ref,1);
    match_dist      = zeros(N_ref,1);
    matches         = zeros(2*N_ref,1);
    match_mult      = zeros(2*N_ref,1);
    N_match = 0;
    
    % find all single peaks that match within mimimum_ref_distance
    for i=1:N_ref
        pr = ref_peaks(i);
        [~,closest_idx] = min(abs(check_peaks-pr));
        dist = abs(pr - check_peaks(closest_idx));
        multsum = ref_mult(i) + check_mult(closest_idx);
        if ( dist < minimum_ref_distance && ~matched_check(closest_idx) && multsum == 2)
            matched_ref(i)              = 1;
            matched_check(closest_idx)  = 1;
            N_match                     = N_match + 1;
            match_index(N_match)        = closest_idx;
            match_dist(N_match)         = dist;
            matches(N_match)            = pr;
            match_mult(N_match)         = 1;
        end
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

    % add all singles from check have no corresponding peak in ref within +-
    % one peakwidth; mark correspondingly
    ic = 0;
    for i = 1:N_check
        if (~logical(matched_check(i)) && check_mult(i) == 1)   % only go over unmatched singles
            pr      = check_peaks(i);
            left    = pr - pw(pr);
            right   = pr + pw(i);
            test    = (ref_peaks > left) & (ref_peaks < right);  % only 1 if there are peaks in range
            if (sum(test) == 1 && ~logical(matched_ref(test)) && (ref_mult(test)==1) )
                N_match                     = N_match + 1;
                matches(N_match)            = pr;
                match_mult(N_match)         = 1;
                matched_check(i)            = 1;
                matched_ref(test)           = 1;
                ic                          = ic + 1;
            end
        end
    end
    
    % reduce to a list only of multiples
    peaks_mult_ref      = ref_peaks(~logical(matched_ref));
    mult_mult_ref       = ref_mult(~logical(matched_ref));
    peaks_mult_check    = check_peaks(~logical(matched_check));
    mult_mult_check     = check_mult(~logical(matched_check));

    % get the (end) flags of the multiples
    flag_mult_ref       = make_flag_from_mult(mult_mult_ref);
    flag_mult_check     = make_flag_from_mult(mult_mult_check);

    % combine into one list sorted by m, first column: m, second column:
    % multiplicity ref, third column: multiplicity check
    comb_mult_list      = nan(length(peaks_mult_ref)+length(peaks_mult_check),5);
    comb_mult_list(:,1) = [peaks_mult_ref ; peaks_mult_check];
    comb_mult_list(1:length(peaks_mult_ref),2) = flag_mult_ref;
    comb_mult_list(1:length(peaks_mult_ref),3) = mult_mult_ref;
    comb_mult_list(length(peaks_mult_ref)+1:end,4) = flag_mult_check;
    comb_mult_list(length(peaks_mult_ref)+1:end,5) = mult_mult_check;
    [~,idx] = sort(comb_mult_list(:,1));
    comb_mult_list = comb_mult_list(idx,:);

    % scan over list and get boundary flag
    ref_state   = nan(length(peaks_mult_ref)+length(peaks_mult_check),1);
    check_state = nan(length(peaks_mult_ref)+length(peaks_mult_check),1);
    final_state = nan(length(peaks_mult_ref)+length(peaks_mult_check),1);
    sw_r = 1;
    sw_c = 1;
    for i=1:length(ref_state)
        if ( ~isnan(comb_mult_list(i,2)))
            ref_state(i) = comb_mult_list(i,2);
            sw_r = comb_mult_list(i,2);
        else
            ref_state(i) = sw_r;
        end
        if ( ~isnan(comb_mult_list(i,4)))
            check_state(i) = comb_mult_list(i,4);
            sw_c = comb_mult_list(i,4);
        else
            check_state(i) = sw_c;
        end
        final_state(i) = ref_state(i) + check_state(i) == 2;
    end

    merged_multiple = [ comb_mult_list ref_state check_state final_state];

    % reduce list by combining peaks within super clusters that are less than
    % COMBINE_FRACTION * peakwidth apart and make a list of peakpos and
    % multiples
    final_multiple_peaks = zeros(size(merged_multiple,1),1);
    final_multiple_multiplicites = zeros(size(merged_multiple,1),1);
    ic  = 1;
    isc = 1;
    np = 0;
    for i = 1:size(merged_multiple,1)    
        %gather the peaks until next 1
        if (final_state(i) == 0)        % keep counting
            np = np + 1;
        else                            % new supercluster
            np = np + 1;
            isc = isc + 1;
            tmp_p = merged_multiple((i-np+1):i,1);
            % now merge peaks 
            tmp_merged = merge_peaks_v2(tmp_p,width_coeffs,COMBINE_FRACTION);
            final_multiple_peaks(ic:(ic-1+length(tmp_merged))) = tmp_merged;
            final_multiple_multiplicites(ic:(ic-1+length(tmp_merged))) = length(tmp_merged);
            ic = ic + length(tmp_merged);  
            np = 0;
        end
    end

    tmp1 = [matches(1:N_match) ; final_multiple_peaks(1:(ic-1))];
    tmp2 = [match_mult(1:N_match) ; final_multiple_multiplicites(1:(ic-1))];
    [~,idx]         = sort(tmp1);
    all_peaks_tmp1  = tmp1(idx);
    all_mults_tmp1  = tmp2(idx);
    

    % make cluster memebership from mults
    [~,clusters] = make_clusters(all_peaks_tmp1, all_mults_tmp1);  
    final_clusters = merge_clusters(clusters,COMBINE_FRACTION, pw);
  
    %convert back to peaks and mults
    npeaks = 0;
    for i=1:length(final_clusters)
        npeaks = npeaks + length(final_clusters(i).pos);
    end
   
    all_peaks = zeros(npeaks,1);
    all_mults = zeros(npeaks,1);
    ic = 0;
    for i = 1:length(final_clusters)
        poses = final_clusters(i).pos;
        all_peaks(ic+1:ic + length(poses)) = poses;
        all_mults(ic+1:ic + length(poses)) = length(poses);
        ic = ic + length(poses);
    end
    
end

