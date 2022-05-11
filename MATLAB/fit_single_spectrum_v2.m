function [current_peak_summary,total_recon,BG1,BGSUB,FINE,BUMPS,CHUNK_indices,CONV] = ...
    fit_single_spectrum_v2( ...
            MZS, ...
            ITS, ...
            width_coeffs, ...
            PARMS)
% fit_single_spectrum first estimates a flat background, generates 
% background subtracted data, estimates BUMPS and the fine structure. The
% BUMPS are used (via their 2nd derivative to seperate the m/z range into
% separate intervals (chunks).
% Looping over chunks it accumulates fits to clusters of asymmetric
% Gaussian peaks calling fit_chunk_aG.
%
% author: HR, MK
% version:  1.0.1 (4/20/2021)
%           1.0.2 (6/3/2021) HR: occasional issues with VERTCAT if tables
%           have no contents; added clause to if and a try catch block
%           around line 125
%           1.0.3 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region. 
% Input:
%           MZS             (Nx1) vector of m/z
%           ITS             (Nx1) vector of intensities
%           width_coeffs    table of peak_shape parameters 
%           PARMS           structure of parameters with fields:
%                           (default values are in 2nd line)
%                 BG_LAMBDA1 Eilers for background                
%                   default 10^11
%                 BG_LAMBDA2 Eilers for background                
%                   default 10^4
%                 BG_P       Eilers for background                 
%                   default 0.001
%                 FINE_LAMBDA1 Eilers for fine structure              
%                   default 10^6
%                 FINE_LAMBDA2 Eilers for fine structure            
%                   default 10^2
%                 FINE_P       Eilers for fine structure             
%                   default 0.001;
%                 BUMP_MINPROMINENCE  to detect chunk boundaries      
%                   default 4*10^(-4)
%                 PEAK_RANGE_FOR_CONVOLUTION how far away to convolute with
%                                           peak shape
%                   default 1.2
%                 MIN_SEP_FRACTION  mimum separation of mimima in convolution        
%                   default 4
%                 SIGNAL_SN_CUTOFF  SN_cutoff_for_signal        
%                   default 10
%                 INFLUENCE_CUTOFF  if the intersection peak value is greater than this * values consider the peaks to be overlapping        
%                   default 0.1
%                 MIN_PEAK_SEPARATION  if peaks are closer than this times the peakwidth they are conisdered to overlap      
%                   default 0.5
%
% Output:
%           current_peak_summary    Table with peaks
%                   variables:
%                       {'index'}   not really used index into pk list
%                       {'m0'}      position of a peaks
%                       {'A'}       peak amplitudes
%                       {'int'}     peak integral
%                       {'SN'}      signal/noise estimator (not scaled!)
%                       {'mult'}    multiplicity of a peak (i.e. size of
%                                   the cluster a peaks belongs to
%           total_recon         the sum of all asymmetric Gaussian peaks
%                               evaluated at the values of MZS
%           BG                  the flat background
%           BGSUB               the background subtracted intensities
%           FINE                the fine structure (BGSUB - BUMPS)
%           BUMPS               a background closer to peaks, i.e. broad
%                               structures varying on a lengthscale much
%                               larger than the peakwidth
%           CHUNK_indices       index into MZS defining the chunks
%           CONV                the convolution of FINE with the peakshape
%
% Dependencies:
%           baseline8_faster_conv
%           convolute_with_PS_v2
%           fit_chunk_aG
%           normalize_spectra_trapezoid
%
% NOTE:     
%           ITER for EILER's is hardcoded to be 40
% Updates
    % Number of iterations to run BG
    iter        = 40;
    
    % estimate background (Eilers' method)
    [BG1,~]      = baseline8_faster_conv(ITS, PARMS.BG_LAMBDA1, PARMS.BG_LAMBDA2, PARMS.BG_P, iter);

    % Calculate Fine and bums with elastic Eilers BG estimation
    [BG2,~]  = baseline8_faster_conv(ITS, PARMS.FINE_LAMBDA1, PARMS.FINE_LAMBDA2, PARMS.FINE_P, iter);
    BGSUB       = ITS - BG1; % Background subtracted spectra
    BUMPS       = BG2 - BG1; %Bumps structure
    FINE        = ITS - BG2; % Fine structure
    
    % generate chunks (fit peaks along natural discontinuities in spectra)
    % normalize BG2 for more consistent computational results
    integrated_zbias = normalize_spectra_trapezoid(MZS,BG2,MZS(1),MZS(end));         
    tmp_zbias  = length(BG2)*BG2/integrated_zbias; 
    % 2nd Derivative of BG2 to find bump discontinuities
    z_dif       = 10000*diff(tmp_zbias,2);                                              
    CHUNK_indices          = find(islocalmax(z_dif,'MinProminence',PARMS.BUMP_MINPROMINENCE));      
    
    % set up peakshape as a function of m/z from Equation 2 of main text
    % Linear quadratic intersection
    m_int = width_coeffs{2,'Cutoff'};
    % Left HWHM
    sig_L        = @(x) (width_coeffs{2,'c0'} + width_coeffs{2,'c1'}*x + width_coeffs{2,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{2,'a0'} + width_coeffs{2,'a1'}*x) .* (x < m_int);
    % Right HWHM
    sig_R        = @(x) (width_coeffs{3,'c0'} + width_coeffs{3,'c1'}*x + width_coeffs{3,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{3,'a0'} + width_coeffs{3,'a1'}*x) .* (x < m_int);
    
    % convolute signal with peakshape
    CONV        = convolute_with_PS_v2(MZS,FINE,sig_L,sig_R,PARMS.PEAK_RANGE_FOR_CONVOLUTION);
    
    % loop over chunks
    N_chunks    = length(CHUNK_indices);
    current_peak_summary = [];
    total_recon = zeros(size(FINE,1),size(FINE,2));
    for ic = 1:N_chunks - 1
        chunk   = [CHUNK_indices(ic) CHUNK_indices(ic+1)];
        mc      = MZS(chunk(1):chunk(2));
        FINE_c  = FINE(chunk(1):chunk(2));
        CONV_c  = CONV(chunk(1):chunk(2));
        
        %check that chunk is wide enough to fit
        chunk_range     = mc(end)-mc(1);
        mz_midp         = 0.5*(mc(end)+mc(1));
        chunk_pw        =  sig_L(mz_midp) + sig_R(mz_midp);
        
        % Fit if chunk is larger than a peak width, otherwise do nothing
        if ( chunk_range > chunk_pw) 
            [peak_summary_chunk,recon] = fit_chunk_aG_v2( ... 
            FINE_c, ...                 % signal 
            mc, ...                     % mz
            CONV_c, ...                 % convolution
            sig_L, ...                   % handle for sig1
            sig_R, ...                   % handle for sig2
            PARMS.MIN_SEP_FRACTION , ...      % mimum separation of mimima in convolution 
            PARMS.SIGNAL_SN_CUTOFF, ...   % order of 10 SN_cutoff_for_signal
            PARMS.INFLUENCE_CUTOFF, ...       % if the intersection peak value for ocmbining peaks to clusters influence_cutoff
            PARMS.MIN_PEAK_SEPARATION);        % if peaks are closer than this tiems the peakwidth they are conisdered to overlap min_peak_separation
            total_recon(chunk(1):chunk(2)) = recon;
        else
            total_recon(chunk(1):chunk(2)) = zeros(length(FINE_c),1);
            peak_summary_chunk = [];
        end
        
        if (isempty( current_peak_summary) || isempty(current_peak_summary{:,'index'}))
              current_peak_summary = peak_summary_chunk;
        else
            try
              current_peak_summary = [current_peak_summary ; peak_summary_chunk];
            catch
                warning(['not updating current_peak_summary with ' num2str(length(current_peak_summary{:,'index'})) ' peaks' ])
                fff = 0;
            end
        end
        if ( ~isempty(peak_summary_chunk))
            clear recon
        end
        clear FINE_c mc CONV_c
    end

end

