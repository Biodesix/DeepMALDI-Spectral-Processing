function [feature_table,recon] = get_featuretable_v2(MZS,ITS,width_coeffs,peaks,mults)
% Applies a set of feature definitions (from peak positions and
% multiplicities) to a spectrum. The peak positions remain fixed, and only
% the amplitudes are allowed to vary depending on the spectrum. The result
% is a table with positions, amplitudes, and integrals over a peak.
%
% author: HR, MK
% version:  1.0.0 (4/21/2021)
%           2.0.1 (6/16/2021) MK: Changing the peakshape function to be
%           piecewise fit. Linear for low mass region and quadratic for
%           high mass region.
%           2.0.2 (6/21/2021) MK: Added a condition to remove antipeak 
%           fits when fitting the spectrum to determine the features of 
%           interest (lines 109-156)
%
% Input: 
%       MZS             vector of a spectrum's mz's
%       ITS             vector of a spectrum's intensities
%       width_coeffs    table of peakshape parameters
%       peaks           vector of peaks
%       mults           vector of corresponding multiplicites 
%
% Output:
%       feature_table:  table with variables:
%                       'm0'        peak position
%                       'height'    peak amplitude
%                       'area'      peak area
%       recon: 
% dependencies
%       SG_filter_ignore_boundaries 
%       make_clustermembership_from_mult
%       m2pts
%       aG
%
        
    N = length(MZS);
    NP = length(peaks);
    RECON_RANGE = 4;    % in units of peak width
    recon = zeros(size(ITS,1),size(ITS,2));
    
    % use a SG smoothed intensity
    SM = SG_filter_ignore_boundaries(ITS,65,16);
    SM = SM';
    
    % make the clusters
    [~,nc] = make_clustermembership_from_mult(mults);
    

    % set up peakshape as a function of m/z (Equation 2 of main text)
    m_int = width_coeffs{2,'Cutoff'};
    % Left and right HWHM
    sig_L        = @(x) (width_coeffs{2,'c0'} + width_coeffs{2,'c1'}*x + width_coeffs{2,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{2,'a0'} + width_coeffs{2,'a1'}*x) .* (x < m_int);
    sig_R        = @(x) (width_coeffs{3,'c0'} + width_coeffs{3,'c1'}*x + width_coeffs{3,'c2'}*x.^2) .* (x>=m_int) ...
        + (width_coeffs{3,'a0'} + width_coeffs{3,'a1'}*x) .* (x < m_int);

    % FWHM
    pw          = @(x) sig_L(x) + sig_R(x) ;
    
    %preset table columns
    m0      = zeros(NP,1);
    height  = zeros(NP,1);
    p_int   = zeros(NP,1);
    i0      = 1;                % initialize first cluster index
    for i = 1:nc
        i1 = i0;
        i2 = i0 + mults(i0) - 1;
        
        % peaks positions are pre-defined based on input
        m0(i1:i2)       = peaks(i1:i2);
        
        % amplitudes need some fitting for amplitudes (only linear though)
        % first find the range of the cluster
        if( i1 == 1)
            left = peaks(i1) - pw(peaks(i1));
        else
            left = min(peaks(i1) - pw(peaks(i1)),0.5*(peaks(i1)+peaks(i1-1)));
        end
        if( i2 == NP)
            right = peaks(i2) + pw(peaks(i2));
        else
            right = min(peaks(i2) + pw(peaks(i2)),0.5*(peaks(i2+1)+peaks(i2)));
        end

        
        [~,idx_l]   = min(abs(MZS-left));
        [~,idx_r]   = min(abs(MZS-right));
        ms          = MZS(idx_l:idx_r);
        vals        = SM(idx_l:idx_r);
        % set up the linear system for the amplitudes
        MM      = mults(i0);
        G_mat   = zeros(MM,MM);
        SG_vect = zeros(MM,1); %MK added the ,1
        gk_s    = zeros(length(ms),MM);
        for l = 1:MM
            gk_s(:,l) = aG(ms,m0(l+i1-1),sig_L(m0(l+i1-1)),sig_R(m0(l+i1-1)));
        end
        for l = 1:MM
            SG_vect(l) = dot(vals,gk_s(:,l));
            for k = 1:MM
                G_mat(l,k) = dot(gk_s(:,k),gk_s(:,l));
            end
        end
        As = G_mat \ SG_vect ;


        % Refit if original fit gives an antipeak (negative amplitudes)
        refit = 1;
        bad_loc = false(length(As),1);
        while min(As)<0 && refit ==1
            refit = 0;

            %Create a reconstruction over the region for the fit
            m = MZS(idx_l:idx_r);
            tmp_recon = zeros(size(m,1),size(m,2));
            for l =1:MM
                fit     = As(l) * aG(m,m0(l+i1-1),sig_L(m0(l+i1-1)),sig_R(m0(l+i1-1)));
                tmp_recon = tmp_recon + fit;
            end
            
            %Find Local minima in the reconstruction
            diff_tmp_recon = [0; diff(tmp_recon,2) ./ diff(m,2); 0];
            midpoint = 0.5*(m(1)+m(end));
            PEAKWIDTH = sig_L(midpoint) + sig_R(midpoint);
            PEAKWIDTH_in_pts = m2pts(m,midpoint,PEAKWIDTH);
            max_diff_tmp_recon = find(islocalmax(diff_tmp_recon,'MinSeparation',PEAKWIDTH_in_pts));
            min_tmp_recon = m(max_diff_tmp_recon);
            
            %Test the negative fits if they are local minima
            ind_As_bad = find(As<0);
            for itest = 1:length(ind_As_bad)
                irefit = ind_As_bad(itest);
                m_bad = m0(i1-1+irefit);
                [~,bad_idx] =  min(abs(min_tmp_recon-m_bad));
                if abs(min_tmp_recon(bad_idx)-m_bad) < pw(m_bad)/2 
                    refit = 1;
                    bad_loc(irefit) = 1;
                end
            end
           
            %Remove any anti-peaks found above and refit remaining peaks
            if refit == 1
                % Redefine the matrices removing the anti-peaks
                SG_vect_refit = SG_vect(~bad_loc);
                G_mat_refit = G_mat(~bad_loc,~bad_loc);
                As_refit = G_mat_refit \ SG_vect_refit ;

                % Set anti-peaks to 0 and change others to new fits
                As(bad_loc) = zeros(sum(bad_loc),1);
                As(~bad_loc) = As_refit;
            end
        end        


        for l = 1:MM
            height(l+i1-1)   = As(l);
            p_int(l+i1-1)    = 0.5*height(l+i1-1)*(pw(m0(l+i1-1)))*sqrt(pi/log(2));
        end       
        
        % now reconstruction
        left        = max(left - RECON_RANGE*pw(m0(i)),MZS(1));
        right       = min(right + RECON_RANGE*pw(m0(i)),MZS(N));
        [~,idx_l]   = min(abs(MZS - left)); 
        [~,idx_r]   = min(abs(MZS - right));
        m       = MZS(idx_l:idx_r);
        for l =1:MM
            fit     = height(l+i1-1) * aG(m,m0(l+i1-1),sig_L(m0(l+i1-1)),sig_R(m0(l+i1-1)));
            recon(idx_l:idx_r) = recon(idx_l:idx_r)+fit;
        end
        i0 = i2 + 1;            % beginning of the next cluster
    end
    
    feature_table = table(m0,height,p_int,'VariableNames',{ 'm0' 'height' 'area'});
    
end

