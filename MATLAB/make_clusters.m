function [c_member,clusters] = make_clusters(peaks,mult)
%
% makes clusters for peaks from their multiplicities. A cluster is a set of
% peaks with the same multiplicity. It assumes that peaks is sorted (not
% checked, but could be done). 
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% Input: 
%       peaks           peaks 
%       mult            corresponding multiplicities
%
% Output:
%       c_member   the cluster number of a peak
%       clusters   structure array of number of clusters with field pos;
%           pos is a vector of the peaks in that cluster
%
% Notes: there are soem debugging variables kept e.g. ic1  

    c_member        = zeros(size(mult,1),size(mult,2));
    cluster_size    = 0;
    cluster_count   = 0;
    nc              = 0;
    ic1             = 0;
    pos(1)          = peaks(1);
    for i=1:length(mult)
        cluster_count = cluster_count + 1;
        if ( cluster_count <= cluster_size )  % within current cluster
            c_member(i)     = nc;
            pos(cluster_count) = peaks(i);
        else                            % new cluster  
            nc                  = nc + 1;
            %save finished cluster away 
            if (nc > 1 )
                clusters(nc-1).pos = pos(1:cluster_size);
                ic1 = ic1 + cluster_size;
            end
            
            % initiate new cluster 
            c_member(i)         = nc;
            cluster_size        = mult(i);
            cluster_count       = 1;
            pos(cluster_count)   = peaks(i);            
        end
        %display([ 'i ' num2str(i) ' cluster_cnt ' num2str(cluster_count) ' cluster_size ' num2str(cluster_size) ' nc ' num2str(nc)])
    end
    % and the last
    clusters(nc).pos = pos(1:cluster_size);
    ic1 = ic1 + cluster_size;
    ic = 0;
    for i=1:length(clusters)
        ic = ic +length(clusters(i).pos);
    end
    
end

