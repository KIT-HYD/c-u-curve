function [Hs] = f_maxEnt_known_mean(states, means)
% Returns the maximum possible entropy of a discrete distribution for which the mean is known. The set of discrete values the data can take has to be known.
% The function applies the semianalytical solution provided by Keith Conrads "Probability distributions and maximum entropy", Example 5.13, based on Theorem 5.12
% - https://kconrad.math.uconn.edu/blurbs/analysis/entropypost.pdf, downloaded 2022/04/29
% - Example solution for checking: means = 4.7, states = 1,2,...6 --> beta = -0.463; Hs = 2.208
% Input
% - states: [1,num_states] ordered set of discrete values (states) that the data coming from the distribution can take
% - means: [num_means,1] array of mean values for each of which the corresponding maximum possible entropy should be found
% Output
% - Hs: [num_means,1] array with the maximum entropy values for each mean value in 'means'
% Dependencies
% - none
% Version
% - 2022/11/01 Uwe Ehret: initial version

% settings
betas = (-100:0.001:100);   % range of tested beta exponents.
ntries = length(betas);     % number of iteration steps
num_means = length(means);      % number of mean values to be used
Hs = NaN(num_means,1);      % container for maximum entropies

% loop over all mean values in 'means'
for z = 1 : num_means
    
    % initialize vector with deviations from the mean
    dev = NaN(ntries,1);    
    
    % loop over all tested beta exponents
    for i = 1 : ntries
        val = sum(states.*exp(-betas(i).*states))/sum(exp(-betas(i).*states));   % calculate estimated mean for a chosen beta. This uses Eq. 5.5 in Conrads paper  
                                                                                 % to calculate bin probabilities, and then multiplies them with the bin values (states)
                                                                                 % to get the mean value of the distribution
        dev(i) = abs(val-means(z));                                              % calculate deviation of the estimated and the known mean
    end

    % find the optimal beta (the one leading to a mean clostes to the known mean)
    indx = find(dev==min(dev),1);   % find the smallest deviation
    beta = betas(indx);             % get the beta for this
    ps = (exp(-beta.*states))/sum(exp(-beta.*states));  % calculate maximum entropy histogram
    Hs(z) = sum(ps.*-log2(ps));                         % calculate entropy of the histogram

end

end