% Uwe Ehret, 2022/11/01
% Script to test the function 'f_maxEnt_known_mean'

clear all
close all
clc

%% simple example for a 6-sided dice
states = [1 2 3 4 5 6];
means = [4.7; 1; 6; 3.5];
Hmax = f_maxEnt_known_mean(states,means);

%% example for using it to plot the upper complexity bound in the c-u-curve
num_bins = 10;      % number of bins for histrograms of values
num_bins_kld = 10;  % number of bins for histograms of entropies

states = linspace(0,log2(num_bins),num_bins_kld); % discrete values the entropy distribution can take

means = (0:0.01:log2(num_bins));    % candidate mean values, covering the uncertainty range
                                    % uncertainty is limited to [0,log2(num_bin)]

Hmax = f_maxEnt_known_mean(states,means);

% plot entropy vs. mean value 
figure
plot(means,Hmax)

% calculate area under the Hmax curve
A_Hmax = trapz(means,Hmax);