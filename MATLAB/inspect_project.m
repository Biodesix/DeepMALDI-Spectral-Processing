function [project_table] = inspect_project(data_path)
%inspect_project reads through the input spectra directory for spectra to
% Input:
%       data_path:          directory containing folder of spectra to be
%                           read in
% Output:
%       project_table:      table with file name and location for all
%                           spectra
% Dependencies:     none
% Revisions:        
%               v1.0, HR 3/11/2021

    ref_moniker = 'SerumP'; 
    % Samples with the name "SerumP..." is how we determine if a sample is
    % a reference (run 2x at beginning and 2x and end of a given batch) to
    % determine any batch effects on the samples tested.
    batches = dir(data_path); %Reads the number of batches in the data_path directory
    batches(1:2) = []; % removes "\." and "\.." items 
    project_table = [];
    
    % Loop over number of batches to analyze
    for ib=1:length(batches)
        batch_path = [data_path '\' batches(ib).name]; 
        samples_per_batch = dir(batch_path); % reads files within batch directory
        samples_per_batch(1:2) = []; % removes "\." and "\.." items 
        n_in_batch = length(samples_per_batch);
        tmp = struct2table(samples_per_batch);
        sampleID = tmp{:,'name'};
        is_ref = zeros(n_in_batch,1); %Checks if it is a reference sample
        BatchID = strings(n_in_batch,1);
        spectrum_path = strings(n_in_batch,1);
        for i = 1:n_in_batch
            BatchID(i) = convertCharsToStrings(batches(ib).name);
            if contains(sampleID(i),ref_moniker)
                is_ref(i) = 1;
            end
            tmp_list = dir([batch_path '\' sampleID{i}]);
            % Creates the full file name to be read in
            spectrum_path(i) = convertCharsToStrings([batch_path '\' sampleID{i} '\' tmp_list(3).name]);
        end
        % Creates output table for a single batch
        tmp = table(sampleID,BatchID,is_ref,spectrum_path);
        % Appends table to final output table
        project_table = [project_table ; tmp];
    end

end

