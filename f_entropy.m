function [H] = f_entropy(data_binned)
% Returns the (joint) entropy of an 1-to-any-dimensional discrete (binned) frequency distribution
% Note: This is done fast as only non-zero bin occupations are considered
% Input
% - data_binned: [num_data, num_dim] matrix, where num_data is the number of data tuples (= sample size),
%   and num_dim is the number of dimensions of the data set (= number of variables)
%   Note
%   - data_binned must be NaN-free
%   - values in data_binned are positive integers, indicating the bin number into which the original data were classified
%   - values in columns of data_binned must be in the same order, data tuples are linked by the same row number
% Output
% - H: [1,1] entropy in [bit]
% Dependencies
% - none
% Version
% - 2022/02/24 Uwe Ehret: Renamed (old name was f_entropy_anyd_fast)
% - 2021/07/12 Uwe Ehret: initial version

% check if 'data' is NaN-free
    if ~isempty(find(isnan(data_binned)))
        error('data contains NaNs')
    end

% find unique data tuples
    [C,ia,ic] = unique(data_binned,'rows','stable');
    % If [C,ia,ic] = unique(A,'rows','stable'), then
    % C: list of unique rows in A, in order of appearance in A. Unique rows rather than unique colum values is ensured by parameter 'rows'
    % ia: C = A(ia), i.e. ia is the index, where in A each row in C occurs the first time (ensured by parameter 'stable') --> ia has the same length as C
    % ic: A = C(ic), i.e ic contains for each row in A the index of the unique row in C --> ic has the same length as A. 

% create probability distribution of unique data tuples in data_binned    
    binlabels = [1:length(ia)];     % as unique data tuples receive labels 1...num_tuples, 'binlabels' is the list of all unique data tuples 
    fs = hist(ic,binlabels);        % calculate frequencies of all unique data tuples
                                    % Note: All values in 'fs' are > 0 (no empty bins)
    ps = fs/sum(fs);                % normalize frequencies to probabilites                

% calculate entropy
    H = sum(ps .* -log2(ps));       

end

