function [L] = n_in_cluster(i,test)
%
% Given an array of logicals , test, and an index into that array, i,
% n_in_cluster returns either 1 (if the index points to a 0) or the number
% of consecutiove 1s starting at i.
%
% author: HR
% version: 1.0.0 (4/21/2021)
%
% input/output: see description
%
% Note: uses recursion
%
    L = 1;
    if ( i == ( length(test) + 1) )
        return ;
    else
        if ( test(i) )
            L = 1 + n_in_cluster(i+1,test);
        else
            return ;
        end
    end
    
end

