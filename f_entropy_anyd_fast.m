function [H] = f_entropy_anyd_fast(data)
% Returns the (joint) entropy of an 1-to-any-dimensional discrete (binned) frequency distribution
% Note: This is done very fast as only non-zero bin occupations are considered
% Input
% - data: [num_data, num_dims] matrix, where num_data is the number of data tuples, and num_dims is the number of dimensions of the data set
%   - data must be NaN-free
%   - values in data are positive integers, indicating the bin number into which the original data were classified
% Output
% - H: [1,1] entropy in [bit]
% Version
% - 2021/07/12 Uwe Ehret: initial version

% check if 'data' is NaN-free
if ~isempty(find(isnan(data)))
    error('data contains NaNs')
end

% find unique data tuples in data
[C,ia,ic] = unique(data,'rows','stable');
% C: list of unique rows in 'data', in order of appearance in 'data'
% ia: C = A(ia), i.e. ia is the index, where in 'data' each row in C occurs the first time (NOT all occurences) --> ia has the same length as C
% ic: same length as 'data'. contains for each row in 'data' the index of the unique row in C
% 'rows': finds and returns unique rows (not single values) in C
% 'stable': order of unique rows in C is the same as in 'data'

bincenters = [1:length(ia)];    % as unique data tuples receive labels 1...num_tuples, 'bincenters' is the list of all data tuple labels 
fs = hist(ic,bincenters);       % calculate frequencies of all unique data tuple
                                % Note: All values in 'fs' are > 0 (no empty bins)
ps = fs/sum(fs);                % normalize frequencies to probabilites

H = sum(ps .* -log2(ps));       % calculate entropy

end

