function [merged] = merge_clusters(clusters,COMBINE_FRACTION,pw)
%
% merges clusters if its lasr memeber is too close, < COMBINE_RACTION *
% peakwidth, to the first memeber of the next.
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% Input:
%           clusters        
%               structure array of clusters with field pos
%           COMBINE_FRACTION
%               fraction of pw at which to combine clusters
%           pw
%               function handle returning the peakwidth
%
% Output:
%           merged
%               merged cluster, structure array of clusters with field pos
%
% Notes: left a display on for merging loops; throws error if more than
% MAXITER merges


    tmpout = clusters;
    reduction = 1;
    MAXITER = 100;
    iter = 0;
    while (reduction > 0)
            iter = iter + 1;
            if (iter > MAXITER)
                error(['merge_clusters: too many iterations > ' num2str(MAXITER) ])
            end
            fc = 0;
            NC = length(tmpout);
            used = zeros(NC,1);
            for ic = 1:NC - 1
                cluster_i = tmpout(ic).pos;
                cluster_ip = tmpout(ic+1).pos;
                last = cluster_i(end);
                first = cluster_ip(1);
                delta = first - last;
                mid = 0.5*(first + last);
                if ( delta < COMBINE_FRACTION*pw(mid) && ~logical(used(ic)))      % merge
                    merged = [cluster_i(1:length(cluster_i)-1) mid cluster_ip(2:end)];
                    fc = fc + 1;
                    used(ic) = 1;
                    used(ic + 1)  = 1;
                    final_clusters(fc).pos = merged;
                elseif ( ~logical(used(ic)) )                   % copy
                    fc = fc + 1;
                    used(ic) = 1;
                    final_clusters(fc).pos = cluster_i;
                end
            end
            % add the last cluster if not used
            if ( ~logical(used(NC)))  
                fc = fc + 1;
                final_clusters(fc).pos = tmpout(NC).pos;
            end
    
            reduction = length(tmpout) - fc;
            tmpout = final_clusters;
            %disp([' iter: = ' num2str(iter) ' fc = ' num2str(fc)])
            clear final_clusters
    end
    
    %disp([' merge_clusters used ' num2str(iter) ' iterations '])
    merged = tmpout;
    
end

