% Uwe Ehret, 2022/01/07
% This provides test application for the c-u-curve method
% Required Matlab products: Matlab 9.9
% For further explanations, please see the related publication below.Please cite this publication when applying the method.
% REF

clearvars
clear all
close all
clc

%% Test data set: 1-d, deterministic

% settings
    nt = 4096;      % number of rows (=time steps) in the data set
    ndim = 1;       % number of colums (=variables) in the data set
    nens = 1;       % number of ensemble members in the data set
    nvb = 10;       % [1,ndim] array with number of equal-size bins for the value range of the data values of each variable
    neb = 12;       % number of equal-size bins for the value range of entropies
    vals_min = 0;   % minimum value (used to re-scale data values)
    vals_max = 1;   % maximum value (used to re-scale data values)
    % array with all time slice widths to be examined
    % - size is [nss,1], with nss being total number of time slicing schemes to be examined
    % - order is ascending, minimum possible value is 1, maximum possible value is nt
    slice_widths = [1 2 4 8 16 32 64 128 256 512 1024 2048 nt]'; 

% create test data set (random normal)
    dummy = 2*randn(nt,1);                      % create data set
    vals = rescale(dummy,vals_min,vals_max);    % normalize to [0,1] range for convenient binning
    
% create edges of value bins
% - [1,ndim] cell array, with a [1,nvb+1] array of bin edges for each variable
    edges_vals = cell(1,ndim);
    edges_vals{1} = linspace(vals_min,vals_max,nvb+1);

% create edges of entropy bins
% - [1,1] cell array, with a [1,neb+1] array of bin edges for the entropy values
% - the possible value range for entropy values is always [0,log2(total number of value bins for the entire ndim-dimensional space of values)]
    tot_nvb = nvb;  % total number of value bins for the entire ndim-dimensional space of values
    edges_entropy = cell(1,1);
    edges_entropy{1} = linspace(0,log2(tot_nvb),neb+1); 

% calculate uncertainty and complexity
    [uncs,comps,ns,all_uncs] = f_c_u_curve(vals, edges_vals, edges_entropy, slice_widths);    
    
% calculate maximum and minimum possible uncertainty and complexity
% - this is useful to put actual uncertainties and complexities into perspective
% - maximum possible uncertainty is reached if the entire ndim-dimensional space of values
%   is filled by a uniform distribution, whose entropy is log2(total number of value bins for the entire ndim-dimensional space of values)
% - maximum possible complexity is reached if the entropies of all time slices of a time slicing scheme are uniformly distributed.
%   In that case, entropy is log2(number of entropy bins).
% - minimum possible value both for uncertainty and complexity is zero
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    max_possible_unc = log2(tot_nvb);   % maximum possible uncertainty
    max_possible_comp = log2(neb);      % maximum posible complexity
    min_possible_unce = 0;              % minimum possible uncertainty
    min_possible_comp = 0;              % minimum posible complexity

%% Test data set: 2-d, deterministic

% settings
    nt = 4096;      % number of rows (=time steps) in the data set
    ndim = 2;       % number of colums (=variables) in the data set
    nens = 1;       % number of ensemble members in the data set
    nvb = [10 10];  % [1,ndim] array with number of equal-size bins for the value range of the data values of each variable
    neb = 12;       % number of equal-size bins for the value range of entropies
    vals_min = 0;   % minimum value (used to re-scale data values)
    vals_max = 1;   % maximum value (used to re-scale data values
    % array with all time-slice widhts to be examined
    % - size is [nss,1], with nss being total number of time slicing schemes to be examined
    % - order is ascending, minimum possible value is 1, maximum possible value is nt
    slice_widths = [1 2 4 8 16 32 64 128 256 512 1024 2048 nt]'; 
    
% create test data set (random uniform)
    dummy = 2*randn(nt,ndim);                   % create data set
    vals = rescale(dummy,vals_min,vals_max);    % normalize to [0,1] range for convenient binning    
    
% create edges of value bins
% - [1,ndim] cell array, with a [1,nvb+1] array of bin edges for each variable
    edges_vals = cell(1,ndim);
    for i = 1 : ndim
        edges_vals{i} = linspace(vals_min,vals_max,nvb(i)+1); 
    end

% create edges of entropy bins
% - [1,1] cell array, with a [1,neb+1] array of bin edges for the entropy values
% - the possible value range for entropy values is always [0,log2(total number of value bins for the entire ndim-dimensional space of values)]
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    edges_entropy = cell(1,1);
    edges_entropy{1} = linspace(0,log2(tot_nvb),neb+1); 
    
% calculate uncertainty and complexity
    [uncs,comps,ns,all_uncs] = f_c_u_curve(vals, edges_vals, edges_entropy, slice_widths);    

% calculate maximum and minimum possible uncertainty and complexity
% - this is useful to put actual uncertainties and complexities into perspective
% - maximum possible uncertainty is reached if the entire ndim-dimensional space of values
%   is filled by a uniform distribution, whose entropy is log2(total number of value bins for the entire ndim-dimensional space of values)
% - maximum possible complexity is reached if the entropies of all time slices of a time slicing scheme are uniformly distributed.
%   In that case, entropy is log2(number of entropy bins).
% - minimum possible value both for uncertainty and complexity is zero
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    max_possible_unc = log2(tot_nvb);       % maximum possible uncertainty
    max_possible_comp = log2(neb);          % maximum posible complexity
    min_possible_unce = 0;                  % minimum possible uncertainty
    min_possible_comp = 0;                  % minimum posible complexity

%% Test data set: 2-d, 3-member ensemble

% settings
    nt = 4096;      % number of rows (=time steps) in the data set
    ndim = 2;       % number of colums (=variables) in the data set
    nens = 3;       % number of ensemble members in the data set
    nvb = [10 10];  % [1,ndim] array with number of equal-size bins for the value range of the data values of each variable
    neb = 12;       % number of equal-size bins for the value range of entropies
    vals_min = 0;   % minimum value (used to re-scale data values)
    vals_max = 1;   % maximum value (used to re-scale data values
    % array with all time-slice widhts to be examined
    % - size is [nss,1], with nss being total number of time-slicing schemes to be examined
    % - order is ascending, minimum possible value is 1, maximum possible value is nt
    slice_widths = [1 2 4 8 16 32 64 128 256 512 1024 2048 nt]';  
    
% create test data set (random uniform)
    dummy = 2*randn(nt,ndim,nens);                 % create data set
    vals = rescale(dummy,vals_min,vals_max);    % normalize to [0,1] range for convenient binning 
    
% create edges of value bins
% - [1,ndim] cell array, with a [1,nvb+1] array of bin edges for each variable
    edges_vals = cell(1,ndim);
    for i = 1 : ndim
        edges_vals{i} = linspace(vals_min,vals_max,nvb(i)+1); 
    end
    
% create edges of entropy bins
% - [1,1] cell array, with a [1,neb+1] array of bin edges for the entropy values
% - the possible value range for entropy values is always [0,log2(total number of value bins for the entire ndim-dimensional space of values)]
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    edges_entropy = cell(1,1);
    edges_entropy{1} = linspace(0,log2(tot_nvb),neb+1); 
    
% calculate uncertainty and complexity
    [uncs,comps,ns,all_uncs] =f_c_u_curve(vals, edges_vals, edges_entropy, slice_widths);    

% calculate maximum and minimum possible uncertainty and complexity
% - this is useful to put actual uncertainties and complexities into perspective
% - maximum possible uncertainty is reached if the entire ndim-dimensional space of values
%   is filled by a uniform distribution, whose entropy is log2(total number of value bins for the entire ndim-dimensional space of values)
% - maximum possible complexity is reached if the entropies of all time slices of a time-slicing scheme are uniformly distributed.
%   In that case, entropy is log2(number of entropy bins).
% - minimum possible value both for uncertainty and complexity is zero
    tot_nvb = prod(cellfun(@length,edges_vals)-1);  % total number of value bins for the entire ndim-dimensional space of values
    max_possible_unc = log2(tot_nvb);   % maximum possible uncertainty
    max_possible_comp = log2(neb);      % maximum posible complexity
    min_possible_unce = 0;              % minimum possible uncertainty
    min_possible_comp = 0;              % minimum posible complexity
     
     
