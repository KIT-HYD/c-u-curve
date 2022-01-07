function [uncs,comps,ns,all_uncs] = f_c_u_curve(vals, edges_vals, edges_entropy, slice_widths)
% Returns uncertainty and complexity for dynamical system characterization 
% Input
% - vals: [nt,ndim,nens] array with data values. nt is number of time steps, ndim is number of variables, nens is number of ensemble members
%   Note: 
%   - must be NaN-free
% - edges_vals: [1,ndim] cell array, with a [1,nvb+1] array of bin edges for each variable. nvb is number of bins for each variable
%   Note: 
%   - These edges are used to calculate histograms of the data
%   - For each variable, the edges must completely cover its entire value range, for all ensemble members
% - edges_entropy: [1,1] cell array, with a [1,neb+1] array of bin edges for entropy values. neb is number of bins for entropy
%   Note:
%   - These edges are used to calculate histograms of the entropy values calculated for all time slices
%   - The edges must completely cover the entire range of entropy values
% - slice_widths: [nss,1] array with all time slice widths to be examined. nss is number of time slicing schemes
%   Note:
%   - order is ascending, minimum possibe value is 1, maximum possible value is nt
% Output
% - uncs: [nss,1] array, with mean entropy (=uncertainty) of all time slices for each of the time slicing schemes in 'slice_widths'
% - comps: [nss,1] array, with entropy of entropies (=complexity) for each of the time slicing schemes in 'slice_widths'
% - ns: [nss,1] array, with number of time slices the time series was split up into, for each of the time slicing schemes in 'slice_widths'
% - all_uncs: [1,nss] cell array, with [1,ns] array of entropies in each time slice for each of the time slicing schemes in 'slice_widths'
% Note
% - Calculation of entropies, and the treatment of ensemble data
%   - Entropy of a set of values is the same as expected Kullback-Leibler divergence (KLD) of each single value in the set
%     (single deterministic value as a reference) and the distrubution of all values (a probabilistic model for the value)
%   - This also applies if for a given time step and variable, not only a single value exists (in case 'nens'=1),
%     but when several values are provided ('nens'>1). In this case, the reference also is a pdf. Each of these ensemble values contributes with weight '1' 
%     to the overall pdf of the distribution, and as before the expected KLD of all values from their joint distribution is equal to the entropy of that distribution.
%   - The only difference is the lower bound of KLD, KLDmin:
%     - for a deterministic data set ('nens'=1), for the smallest possible time slice (1 time step = 1 time slice), 
%       KLD will by definition be zero, as the within-slice pdf is formed by a single value for each variable.
%     - for an ensemble data set ('nens'>1), for the smallest possible time slice (1 time step = 1 time slice),
%       KLD will only be zero if all ensemble members show the same value, but typically it will be >0.
% Required Matlab products
% - Matlab 9.9
% Version
% - 2022/01/07 Uwe Ehret: initial version
% Reference
% - For further explanations, please see the related publication below. Please cite this publication when applying the method.
% REF

% get dimensions
    [nt, ndim, nens] = size(vals);  % number of time steps, number of variables, number of ensemble members        
    nss = length(slice_widths);     % number of time-slicing schemes
    
% check if there are too-small time-slices
    if any(slice_widths < 1)
        error('there are time slices < 1!')
    end
    
% check if there are too-large time-slices
    if any(slice_widths > nt)
        error('there are time slices > nt!')
    end
    
% check input data for NaN
    if any(isnan(vals(:)))   
        error('input data contain NaN');
    end
    
% check if input data fall outside the bin edges
    mins = min(min(vals,[],1),[],3);   % smallest value per variable (across all time steps and ensemble members)
    maxs = min(max(vals,[],1),[],3);   % smallest value per variable (across all time steps and ensemble members)

    % loop over all variables
    for d = 1 : ndim 
        if mins(d) < edges_vals{d}(1) % smallest value is < leftmost edge
            error('input data < leftmost edge');
        elseif maxs(d) > edges_vals{d}(end)  % largest value is > rightmost edge
            error('input data > rightmost edge');
        end
    end

% check if binning is adequate and issue a warning if not
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    neb = length(edges_entropy{1});                 % number of entropy bins
    % loop over all time-slicing schemes
    for ss = 1 : length(slice_widths)
        if (nt/(neb*3)) < slice_widths(ss)
            bla = strcat ("Warning: slice width ",num2str(slice_widths(ss))," too wide to get sufficient population of distribution of entropies (3 values per bin on average)!");
            disp (bla);
        elseif (tot_nvb*3) > slice_widths(ss)
            bla = strcat ("Warning: slice width ",num2str(slice_widths(ss))," too narrow to get sufficient population of within-slice distribution of values (3 values per bin on average)!");
            disp (bla);
        end
    end
   
% initialize output
    uncs = NaN(nss,1);
    comps = NaN(nss,1);
    ns = NaN(nss,1);
    all_uncs = cell(1,nss);
        
% loop over all time-slicing schemes
for ss = 1 : nss
    
    sw = slice_widths(ss);      % width of a time slice
    num_slices = fix(nt/sw);    % number of time slices
    dummy_entropies = NaN(num_slices,1);    % dummy container for entropies
    
    % loop over all time slices of the current time slicing scheme
    for s = 1 : num_slices
        
        from = 1 + (s-1) * sw;                  % starting time of the current time slice
        vals_slice = vals(from:from+sw-1,:,:);  % extract data in the time slice
        
        % reshape the time slice (remove the ensemble dimension)
        dummy = permute(vals_slice,[1 3 2]);        % swap the 2nd dimension (variables) and 3d dimension (ensemble members)
        vals_reshaped = reshape(dummy,[],ndim,1);   % remove the ensemble dimension, glue them to the lower end of the 2-d matrix (time steps, variables)
        
        % discretize the values (replace values by the number of the bin they fall into)
            vals_discretized = NaN(size(vals_reshaped)); % create the container
        
            % loop over all variables
            for d = 1 : ndim
                vals_discretized(:,d) = discretize(vals_reshaped(:,d),edges_vals{d}); % replace values by the number of the bin they fall into
            end
        
            dummy_entropies(s) = f_entropy_anyd_fast(vals_discretized); % calculate within-slice varibility (uncertainty)
        
    end

    % calculate uncertainty
    all_uncs{ss} = dummy_entropies;     % save all time slice entropies of the the current time slicing scheme
    uncs(ss) = mean(dummy_entropies);   % save the mean of all time slice entropies (uncertainty)
    
    % calculate complexity
    entropies_discretized = discretize(dummy_entropies,edges_entropy{1});   % replace entropy values by the number of the entropy bin the fall into
    comps(ss) = f_entropy_anyd_fast(entropies_discretized);                 % calculate 'entropy of entropies' (complexity)   
    
    ns(ss) = num_slices;    % save number of time slices
end


end

