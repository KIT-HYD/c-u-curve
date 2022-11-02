function [data_binned] = f_binme(data,edges)
% Returns values in 'data' classified/binned/discretized by the bins provided in 'edges'
% Note: This task can be done by using 'histcounts' or 'discretize'. According to Matlab help, use histcounts to find the number of elements in each bin, 
%       and discretize to find which bin each element belongs to (which is the case here). Computation time differences are small, though. 
% Input
% - data: [num_data, num_dim, num_ens] matrix, where num_data is the number of data tuples (= sample size), 
%   num_dim is the number of dimensions of the data set (= number of variables), 
%   and num_ens is the number of ensemble members (= number of values for each variable in a particular data tuple)
%   - data must be NaN-free
%   - the third dimension is optional (no problem if absent)
% - edges [1,num_dim] cell array, with arrays of bin edges for each dimension
%   - the bins of each dimension must completely cover the corresponding values in 'data'
% Output
% - data_binned [num_data, num_dim, num_ens] matrix (same as 'data'), with strictly positive integers, indicating the bin number into which the original data were classified
% Dependencies
% - none
% Version
% - 2022/03/29 Uwe Ehret: Added handling the 'ensemble' dimension in 'data'
% - 2022/03/09 Uwe Ehret: Initial version

% get dimensions and initialize output
    [num_data, num_dim, num_ens] = size(data);  % number of time steps, number of variables, number of ensemble members        
    data_binned = NaN(num_data, num_dim, num_ens); % initialize output

% reshape the data (remove the ensemble dimension)
% Note: This was used in an older version, but not needed any more
%     dummy = permute(data,[1 3 2]);        % swap the 2nd dimension (variables) and 3d dimension (ensemble members)
%     data_reshaped = reshape(dummy,[],num_dim);   % remove the ensemble dimension, glue them to the lower end of the 2-d matrix (time steps, variables)
    
% check input data for NaN
    if any(isnan(data),'all')   
        error('input data contain NaN');
    end

% check if input data fall outside the bin edges
    
    for i = 1 : num_dim     % loop over all variables
        if min(data(:,i,:),[],'all') < edges{i}(1) 
            error(strcat("there are values in dim ", num2str(i), " < lowermost bin edge"));
        end
        if max(data(:,i,:),[],'all') > edges{i}(end) 
            error(strcat("there are values in dim ", num2str(i), " > uppermost bin edge"));
        end
    end

% discretize the values (replace values by the number of the bin they fall into)
    
    for i = 1 : num_dim     % loop over all variables
        data_binned(:,i,:) = discretize(data(:,i,:),edges{i});
    end

end