function [spectrum] = read_spectrum(spectrum_file)
% Reads a Biodesix spectrum file
%   limits input to < 30kDa
% Input:
%       spectrum_file:      fully qualified filename
% Output:
%       spectrum:           structure with elements mz and its, both vectors 
% Usage:
%       SpectrumFile = 'C:\Projects\GoldStandard\SP5\Data\SP5_avg_perm_1_NextgenProcessedSpectra_00001.txt';
%       spectrum = read_spectrum(SpectrumFile);
% Dependencies:     none
% Revisions:        
%               v1.0, HR 3/11/2021

    fid = fopen(spectrum_file);
    tline = fgetl(fid);
    start_read_switch = 0;
    ic = 0;
    
    % Reads through Biodesix file to extract m/z and intensity
    while ischar(tline)
        % If we have passed the header, read spectra and output
        if (start_read_switch )
            ic = ic + 1;
            vals = sscanf(tline,'%g',2);
            if (vals(1) > 30000)        % don't look at > 30kDa
                break;
            end
            spectrum.mz(ic) = vals(1); % m/z location
            spectrum.its(ic) = vals(2); % intensity
        % If still in the header, skip all header lines
        else
            if tline(1) == '#' % Character that signals we have reached the spectral data
                start_read_switch = 1;
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
end

