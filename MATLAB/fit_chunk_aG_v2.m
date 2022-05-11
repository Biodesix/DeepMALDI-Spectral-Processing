function [peak_summary_chunk,recon] = fit_chunk_aG_v2( ...
    s, ...                      % signal 
    m, ...                      % mz
    c, ...                      % convolution
    sig_L, ...                  % handle for left HWHM
    sig_R, ...                  % handle for right HWHM
    MIN_SEP_FRACTION , ...      % mimum separation of mimima in convolution
    SN_cutoff_for_signal, ...   % order of 10
    influence_cutoff, ...       % if the intersection peak value for ocmbining peaks to clusters
    min_peak_separation)        % if peaks are closer than this times the peakwidth they are considered to overlap

% fit_chunk_aG analyses a chunk. 
% it finds peak candidates from the 2nd (smoothed) derivative of the
% convolution with the peakshape. 
% peaks are deemed reasonable candidates by comparing their S/N value (in
% undefined units) to a cut-off, and bracketed related to the local peak
% width.
% clusters of of peaks are defined via their overlap and their sepration.
% Sets of asymmetric Gaussians (set size defined by overlap) are fitted
% using nlinfit (warnings off!).
% 
% author: HR
% version:  1.0.2 (4/20/2021)
%           2.0.0 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region.
%
% input:     (see function signature)
%
% output:
%   peak_summary_chunk      table
%                   variables:
%                       {'index'}   not really used index into pk list
%                       {'m0'}      position of a peaks
%                       {'A'}       peak amplitudes
%                       {'int'}     peak integral
%                       {'SN'}      signal/noise estimator (not scaled!)
%                       {'mult'}    multiplicity of a peak (i.e. size of
%                                   the cluster a peaks belongs to
%
%   recon           the reconstructed signal using the fit parameters 
%                   returns s if there are no peaks in chunk
% 
% dependencies:
%       noise2percentiles
%       find_max_toleft
%       find_max_toright
%       aG
%       n_in_cluster
%       N_aG_peaks
%
% notes:
%       uses smoothdata and islocalmin, islocalmax MATLAB functions; any
%       rewrite needs to follow these very closely as parameters are tricky
%
%       throws a warning if there is nonsense coming out of nlinfit (does
%       not happen more than a couple of times)
%
%       throws an error of the peak m0s are not sorted (should never
%       happen)

%% make guess for peakwidths from mid_point of chunk
% Because chunks are relatively small, the peak width is not going to
% change significantly over the range. If larger ranges are used, then this
% method should be refined.
    midpoint = 0.5*(m(1)+m(end));
    PEAKWIDTH = sig_L(midpoint) + sig_R(midpoint);
    
    % Calculate peak width in data points
    if ( midpoint < (m(1)+ PEAKWIDTH ) )
        midpoint = m(1)+ PEAKWIDTH;
    end
    
    if ( midpoint > (m(length(m)) - PEAKWIDTH) )
        midpoint = m(length(m)) - PEAKWIDTH;
    end
    
    [~,idxl] = min( abs(m-midpoint + PEAKWIDTH/2) );
    [~,idxr] = min( abs(m-midpoint - PEAKWIDTH/2) );
    
    PEAKWIDTH_in_pts = idxr - idxl;

    
%% get diagnostic information on chunk
    % make a smoothed version of signal
    s_smoothed = smoothdata(s,'gaussian',round(PEAKWIDTH_in_pts/3));

    % generate 2nd der of conv and smoothed 2nd der of conv
    chunk_2nd_der = [0 ; diff(c,2) ./ diff(m,2) ; 0];
    chunk_2nd_der_smoothed = smoothdata(chunk_2nd_der,'gaussian',round(PEAKWIDTH_in_pts/2));

    % get list of mimima of smoothed 2nd der of conv
    chunk_2nd_der_smoothed_minima = ... 
        find(islocalmin(chunk_2nd_der_smoothed,'ProminenceWindow',PEAKWIDTH_in_pts,'MinSeparation',round(PEAKWIDTH_in_pts/MIN_SEP_FRACTION)));
    
    % get list of maxima of smoothed 2nd der of conv
    chunk_2nd_der_smoothed_maxima = ... 
        find(islocalmax(chunk_2nd_der_smoothed,'ProminenceWindow',PEAKWIDTH_in_pts,'MinSeparation',round(PEAKWIDTH_in_pts/MIN_SEP_FRACTION)));

    % make a noise estimator by smoothing signal and subtracting from signal 
    noise = s - smoothdata(s,'gaussian',PEAKWIDTH_in_pts);
    % running percentile of abs(noise)
    noisep = noise2percentiles(abs(noise),round(PEAKWIDTH_in_pts),50);

    % estimated S/N ratio using smoothed
    int_sn = abs(s_smoothed ./ noisep) ;
    
%% identify peaks
    % make a list of peaks with boundaries
    n_try = length(chunk_2nd_der_smoothed_minima);
    % initialize peak candidate storage
    cand_mz         = zeros(n_try,1);
    cand_idx        = zeros(n_try,1);
    cand_left       = zeros(n_try,1);
    cand_right      = zeros(n_try,1);
    cand_left_idx   = zeros(n_try,1);
    cand_right_idx  = zeros(n_try,1);
    cand_sn         = zeros(n_try,1);
    n_cand          = 0;
    
    % Loop over number of potential peaks to fit
    for i = 1:n_try
        current = chunk_2nd_der_smoothed_minima(i); 
        % Only save peaks with SNR > cutoff
        if ( ~(int_sn(current)<SN_cutoff_for_signal) )

           % Find the indeces of the peaks
           % Finds the closest maxima towards lower mass
           % If multiple maxima in the chunk
           if ( length(chunk_2nd_der_smoothed_maxima) > 1)
               
                [~,left_idx] = find_max_toleft(m(current),m,chunk_2nd_der_smoothed_maxima);
           % Otherwise go half a peak width away
           else
                left_idx = (current - round(PEAKWIDTH_in_pts/2));
           end
           left_idx = min( (current - round(PEAKWIDTH_in_pts/2)), left_idx);  
           left_idx = max(1,left_idx);
           left = m(left_idx);
           
           % Finds the closest maxima towards higher mass
           % If multiple maxima in the chunk
           if ( length(chunk_2nd_der_smoothed_maxima) > 1)
                [~,right_idx] = find_max_toright(m(current),m,chunk_2nd_der_smoothed_maxima);
           % Otherwise go half a peak width away
           else
                right_idx = (current + round(PEAKWIDTH_in_pts/2));
           end
           right_idx = max( (current + round(PEAKWIDTH_in_pts/2)), right_idx);  
           right_idx = min(right_idx,length(m));
           right = m(right_idx);
           
           % store away peak candidates
           n_cand = n_cand + 1;
           cand_mz(n_cand)          = m(current);
           cand_idx(n_cand)         = current;
           cand_right(n_cand)       = right;
           cand_left(n_cand)        = left;
           cand_right_idx(n_cand)   = right_idx;
           cand_left_idx(n_cand)    = left_idx;
           cand_sn(n_cand)          = int_sn(current);
           
        end
    end
    
%% stop if no peaks

    if (n_cand == 0 )
        recon = s ;
        peak_summary_chunk = [];
        return;
    end

%% get peak to peak influence to define clusters
    peak2nextpeak_influence = zeros(n_cand-1,1);
    nextpeak2peak_influence = zeros(n_cand-1,1);
    too_close = zeros(n_cand-1,1);
    % Check for all peaks in chunk if they are influenced
    for i = 1:(n_cand-1)
        % Peaks are too close if their centers are within
        % min_peak_separation*peak width
        too_close(i) = (cand_idx(i+1)-cand_idx(i)) < min_peak_separation*PEAKWIDTH_in_pts;
        
        %If peaks are far enough away check for influence (one is much
        %larger than the other and affects potential fitting of single
        %peak)
        if ( ~too_close(i))
            m1          = cand_mz(i);
            m2          = cand_mz(i+1);
            w1_l        = sig_L(m1);
            w1_r        = sig_R(m1);
            w2_l        = sig_L(m2);
            w2_r        = sig_R(m2);
            A1          = s_smoothed(cand_idx(i));
            A2          = s_smoothed(cand_idx(i+1));

            first_mz    = m1;
            second_mz   = m2;

            % Find where the two peaks intersect in intensity 
            %(peak1 - peak2)
            diff_p      = @(x) A2*aG(x,m2,w2_l,w2_r) - A1*aG(x,m1,w1_l,w1_r);
            
            % if peaks are similar intensity/spaced far apart
            if (sign( diff_p(first_mz)*diff_p(second_mz)) <0)
                intercept = fzero(diff_p,[first_mz second_mz]);
                intercept_value = A2*aG(intercept,m2,w2_l,w2_r);
                % now influence functions (peaks may be influenced based on
                % cutoff)
                nextpeak2peak_influence(i) = intercept_value/A1;
                peak2nextpeak_influence(i) = intercept_value/A2;
                
            % Otherwise, one peak is too large and tail of one peak is
            % greater in intensity than the other peak's max. 
            else
                nextpeak2peak_influence(i) = 1; %Peaks are influenced
                peak2nextpeak_influence(i) = 1; %Peaks are influenced
            end
            
        %Peaks too close and are thus influenced
        else
            nextpeak2peak_influence(i) = 1;
            peak2nextpeak_influence(i) = 1;
        end
    end
    tmp1 = nextpeak2peak_influence > influence_cutoff;
    tmp2 = peak2nextpeak_influence > influence_cutoff;
    test_nextpeak = tmp1 | tmp2;  % if either peak influencs the other
    test_nextpeak = test_nextpeak | too_close ; 
    
%% Now to fit peaks to determine exact m/z positions and amplitudes
    i_cand = 1;
    ic = 0;
    ip = 0;
    peak_mul    = zeros(n_cand,1);
    peak_m0     = zeros(n_cand,1);
    peak_A      = zeros(n_cand,1);
    peak_SN     = zeros(n_cand,1);
    peak_int    = zeros(n_cand,1);
    peak_idx    = zeros(n_cand,1);
    
    %Fit total number of candidate peaks in list. 
    % Isolated peaks are fit to single gaussian, influenced peaks are fit
    % to N gaussians based on multiplicity.
    while (i_cand < n_cand +1)
        n_peaks_in_cluster = n_in_cluster(i_cand,test_nextpeak);
        ic = ic + 1;
        % find boundaries of cluster    
        bleft = cand_left(i_cand);
        bright = cand_right(i_cand + n_peaks_in_cluster-1);
        bleft_idx = cand_left_idx(i_cand);
        bright_idx = cand_right_idx(i_cand+n_peaks_in_cluster-1);

        % get the signal and m range 
        s_cluster = s(bleft_idx:bright_idx);
        m_cluster = m(bleft_idx:bright_idx);
        
        % make the initial conditions
        b = zeros(2*n_peaks_in_cluster,1);
        for i = 1:n_peaks_in_cluster
            b(2*(i-1)+1)    = cand_mz(i_cand + (i-1));
            b(2*i)          = s_smoothed( cand_idx(i_cand + (i-1)));
        end
        
        % make a handle for fitting and fit
        fit_handle = @(b,x) N_aG_peaks_v2(b,x,n_peaks_in_cluster,sig_L,sig_R);%%%%%%%%
        warning('off')
        b_out = nlinfit(m_cluster,s_cluster,fit_handle,b);
        warning('on')
        
        % add some checking as ignoring warnings
        % only add to list if the following is fulfilled
        % are the peaks sorted and are they within the bounds [bleft,bright]
            % make a tmp list of peaks
            tmp = zeros(n_peaks_in_cluster,1);
            for i = 1:n_peaks_in_cluster
                tmp(i) = b_out( 2*(i-1) + 1);
            end
            fail = ~issorted(tmp);
            if (~fail)
                if (tmp(1) < bleft || tmp(end) > bright)
                    fail = 1;
                end
            end
        if (~fail)    
            for i = 1:n_peaks_in_cluster
                ip = ip + 1;
                peak_mul(ip)    = n_peaks_in_cluster;
                peak_m0(ip)     = b_out( 2*(i-1) + 1);
                peak_A(ip)      = b_out( 2*i); 
                [~,tmp] = min(abs(m-peak_m0(ip)));
                peak_idx(ip)    = tmp;
                peak_SN(ip)     = int_sn(tmp);
                peak_int(ip)    = 0.5*peak_A(ip)*(sig_L(peak_m0(ip))+sig_R(peak_m0(ip)))*sqrt(pi/log(2));
            end
        else      
            warning(['in fit_chunk_aG: skipped peak cluster for nonsense fit from nlinfit: ' num2str(bleft) ' - ' num2str(bright)])
        end
        i_cand = i_cand + n_peaks_in_cluster;
        if (i_cand > n_cand)
            break;
        end
    end
    
    % resize peak info arrays
    peak_m0     = peak_m0(1:ip);
    peak_mul    = peak_mul(1:ip);
    peak_A      = peak_A(1:ip);
    peak_idx    = peak_idx(1:ip);
    peak_SN     = peak_SN(1:ip);
    peak_int    = peak_int(1:ip);

    % reconstructed signal
    b = zeros(2*length(peak_m0),1);
    for i = 1:length(peak_m0)
        b( 2*(i-1) + 1) = peak_m0(i);
        b( 2*i )        = peak_A(i);
    end
    
    % Reconstructed spectra based on fitting
    recon = N_aG_peaks_v2(b,m,length(peak_m0),sig_L,sig_R);

    %check that peaks are sorted
    if ( ~issorted(peak_m0) )
        error(['fit_chunk_aG produces a peak list that is not sorted for chunk ' num2str(m(1)) ' - ' num2str(m(end)) ])
    end
    
    peak_summary_chunk = table(peak_idx,peak_m0,peak_A,peak_int,peak_SN,peak_mul, ...
    'VariableNames',{ 'index' 'm0' 'A' 'int' 'SN' 'mult'});
    
end

