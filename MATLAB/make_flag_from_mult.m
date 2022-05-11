function [c_flag,nc] = make_flag_from_mult(mult)
%
% loops backwards over multiplicities and generates a logical array with a
% 1 indicating the end of a cluster. also returns the number of clusters.
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% Input:
%       mult        multiplicities of a feature definition
%
% Output:
%       c_flag      vector of logicals (see description)
%       nc          number of clusters
%  
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
    c_flag = -17*ones(size(mult,1),size(mult,2));
    c_flag(length(c_flag)) = 1;
    for i = (length(c_flag)-1):-1:1
        if (c_member(i+1) == c_member(i) )
            c_flag(i) = 0;
        else
            c_flag(i) = 1;
        end
    end
end

