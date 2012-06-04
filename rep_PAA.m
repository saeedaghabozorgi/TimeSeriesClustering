
function [PAA] =  rep_PAA(data, N, n, alphabet_size)

if nargin < 4
    disp('usage: sax_modified(data, window_len, num_segment, alphabet_size, [numerosity_reduction_option]');
    return;
end

%if (N/n - floor(N/n))                               % N/n must be an integer.
%    disp('N/n must be an integer. Aborting '); , return;
%end;

if alphabet_size > 32
    disp('Currently alphabet_size cannot be larger than 20.  Please update the breakpoint table if you wish to do so');
    return;
end

if nargin < 5
    NR_opt = 2;
end

win_size = floor(N/n);                              % win_size is the number of data points on the raw time series that will be mapped to a single symbol

pointers         = [];                                                  % Initialize pointers,
symbolic_data = zeros(1,n);                                             % Initialize symbolic_data with a void string, it will be removed later.
all_string = zeros(length(data)-N+1,n);


% Scan accross the time series extract sub sequences, and converting them to strings.
for i = 1 : length(data) - (N -1)
    
    if mod(i, 1000) == 0
        disp(num2str(i));
    end
    
    % Remove the current subsection.
    sub_section = data(i:i + N -1);
    
    % Z normalize it.
    sub_section = (sub_section - mean(sub_section))/std(sub_section);
    
    % take care of the special case where there is no dimensionality reduction
    if N == n
        PAA = sub_section;
        
        % N is not dividable by n
    else
        if (N/n - floor(N/n))
            temp = zeros(n, N);
            for j = 1 : n
                temp(j, :) = sub_section;
            end
            expanded_sub_section = reshape(temp, 1, N*n);
            PAA = [mean(reshape(expanded_sub_section, N, n))];
            % N is dividable by n
        else
            PAA = [mean(reshape(sub_section,win_size,n))];
        end
        % Convert to PAA.
        %else
        %     PAA = [mean(reshape(sub_section,win_size,n))];
    end
end

