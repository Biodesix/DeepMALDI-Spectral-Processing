function [shifted] = SG_filter_ignore_boundaries(x,window_length,order)
% Does Savatsky-Golay smoothing of a signal; ignores boundaries, i.e. is
% really only defined in the range ((window_length-1)/2+1:(length(x)-(window_length-1)/2
% modified from MAtlab file exchange at:
% https://www.mathworks.com/matlabcentral/answers/335433-how-to-implement-savitzky-golay-filter-without-using-inbuilt-functions
%
% Input:
%       x:              input signal
%       window_length:  length of smoothing window, needs to be odd
%       order:          integer, order of the SG polynomial
% Output:
%       shifted:        smoothed signal, the first (window_length-1)/2
%                       values are 0
% Usage:
%       for example to estimate the noise via x - SG_filter_ignore_boundaries(x,window_length,order)

% Dependencies:     none
% Revisions:        
%               v1.0, HR 3/16/2021


    N1 = length(x);% length of signal
    n = (1:N1); % time vector
    
    %%construction of Savitzky-Golay Filter
    WinL = window_length;   %in samples
    Ord = order; % order of the filter
    shiftL = 1; % hop size in samples
    nFr  = round(length(x)/shiftL); %no., of frames
    WIND = zeros(WinL,nFr);
    for c = 1:nFr - round(WinL/shiftL)
      FB = (c-1)*shiftL+1; % beginning of the frame in samples
      FE = FB + WinL -1;   % ending of the frame in samples
      WIND(:,c) = x(FB:FE);
    end
    for c = 1:nFr - round(WinL/shiftL) % computing no., of frames into windows
      FB = (c-1)*shiftL+1; % beginning of the frame in samples
      FE = FB + WinL -1;   % ending of the frame in samples
      N(:,c) = n(FB:FE);
    end
    %adj = zeros(WinL,size(WIND,2)-size(N,2));
    WIND(:,[size(N,2)+1:size(WIND,2)]) = [];
    polcoeff = zeros(Ord+1,size(N,2)); % coefficients of the polynomial
    polvalues = zeros(WinL,size(N,2)); % value of the function with 'p' polynomial coefficient
    for c = 1:size(N,2)
        t = N(:,c);
        [p,s,mu] = polyfit(t,WIND(:,c),Ord);
        polcoeff(:,c) = p;
        polvalues(:,c) = polyval(p,t(round(WinL/2)),s,mu);
    end
    polvalues(2:WinL,:) = [];
    polvalues = [polvalues,zeros(1,WinL)];
    
    shift = (window_length-1)/2;
    shifted = zeros(1,length(x));
    shifted((shift+1):end) = polvalues(1:(N1-shift));
    %     plot(spectrum.mz(chunk(1):chunk(2)),x(chunk(1):chunk(2)),'-r',spectrum.mz((chunk(1)):(chunk(2))),shifted(chunk(1):chunk(2)),'-b');
    %     ylabel('Amplitude'),xlabel('Number of samples');title('Original signal + Noise');
    %     %axis([0 4096 -1.5 1.5]);
    %     subplot(312)
    %     SG_noise = x - shifted;
    %     plot(spectrum.mz(chunk(1):chunk(2)),SG_noise(chunk(1):chunk(2)));
    %     ylabel('Amplitude'),xlabel('Number of samples');title('Noise');

end

