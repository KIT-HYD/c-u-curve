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
num_means = length(means);  % number of mean values to be used
Hs = NaN(num_means,1);      % container for maximum entropies

% loop over all mean values in 'means'
for z = 1 : num_means

    mean = means(z);            % pick the current mean
    prob = optimproblem;        % Initialize optimization problem
    beta = optimvar("beta");    % Define optimization variable
    prob.Objective =(sum(states.*exp(-beta*states))/sum(exp(-beta*states))-mean)^2; % Define the optimization fuction
    x0.beta = 0;                % define initial value to start iteration
    sol = solve(prob, x0);      % solve optimization problem
    beta = sol.beta;            % get the current optimized beta
    ps = (exp(-beta.*states))/sum(exp(-beta.*states));  % calculate maximum entropy histogram
    Hs(z) = sum(ps.*-log2(ps));                         % calculate entropy of the histogram
    
end

end