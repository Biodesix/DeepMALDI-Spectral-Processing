function [c_member,nc] = make_clustermembership_from_mult(mult)
    % Determines the number of cluseters from the input multiplicity list
    c_member        = zeros(size(mult,1),size(mult,2));
    cluster_size    = 0;
    cluster_count   = 0;
    nc              = 0;
    for i=1:length(mult)
        cluster_count = cluster_count + 1;
        if ( cluster_count <= cluster_size )  % within current cluster
            c_member(i)     = nc;
        else                            % new cluster
            nc                  = nc + 1;
            c_member(i)         = nc;
            cluster_size        = mult(i);
            cluster_count       = 1;
        end
    end
end

