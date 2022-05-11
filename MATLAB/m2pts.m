function [w_in_pts] = m2pts(m,m0,w_in_mz)
%m2pts calculates the width in points of a peak centered at m0 with width
%w_in_mz over the range m.
    
    if ( m0 < (m(1)+ w_in_mz ) )
        m0 = m(1)+ w_in_mz;
    end
    
    if ( m0 > (m(length(m)) - w_in_mz) )
        m0 = m(length(m)) - w_in_mz;
    end
    
    [~,idxl] = min( abs(m-m0 + w_in_mz/2) );
    [~,idxr] = min( abs(m-m0 - w_in_mz/2) );
    
    w_in_pts = idxr - idxl;
end

