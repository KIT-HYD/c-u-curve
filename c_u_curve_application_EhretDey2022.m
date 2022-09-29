% Uwe Ehret, 2022/09/29
% Creates all plots for Ehret and Dey (2022)
% Reference:
% Ehret, U., and P. Dey (2022), Technical note: c-u-curve: A method to analyse, classify and compare dynamical systems by uncertainty and complexity, Hydrol. Earth Syst. Sci. Discuss., 2022, 1-12.

clearvars
clc
close all

%% settings

num_data_arti = 30000;  % for artificial data: number of rows (=time steps) in the test data set 
num_data_real = 12418;  % for real-world data: number of rows (=time steps) in the test data set
num_dim = 1;            % number of colums (=variables) in the test data set
num_ens = 1;            % number of ensemble members
vals_min = 0;           % minimum value in all test data sets
vals_max = 1;           % maximum value in all test data sets
num_bins = 10;          % number of equal-size bins the value range of the data will be split in to calculate histograms
num_bins_kld = 10;      % number of equal-size bins the kld range will be split in to calculate the entropy of the kld distribution

% bin_edges: [1,num_dim] cell array, with a [1,num_bins+1] array of bin edges for each dimension inside
bin_edges = cell(1,num_dim);
bin_edges{1} = linspace(vals_min,vals_max,num_bins+1);

% bin_edges_kld: [1,1] cell array, with a [1,num_bins_kld+1] array of bin edges for the 1-d histogram of kld values
% - the possible value range for kld is always
%   - for a 1-d data set: [0,log2(num_bins)] 
%   - for a 2-d data set: [0,log2(num_bins_in_first_dimension * num_bins_in_second_dimension)]
%   --> set the smallest bin edge 0, the upper to the upper value of the value range
bin_edges_kld = cell(1,1);
bin_edges_kld{1} = linspace(0,log2(num_bins),num_bins_kld+1);

% array with all time-slice widhts to be examined (for which uncertainty and complexity are to be calculated)
% - size is [nss,1], with nss being total number of time-slicing schemes to be examined
% - order is ascending, minimum possible value is 1, maximum possible value is nt
slice_widths_arti = [1 30:10:90 100:50:200 300:100:500 1000 30000]; % for artificial data
slice_widths_real = [1 7 14 21 30 60 91 182 365 730 6209 12418]';   % for real-world data

%% create test data sets

% artifical data
    num_data = num_data_arti;

    % horizontal line (1-d)
        data_line = zeros(num_data,1) + mean([vals_min vals_max]); 
        
    % white noise (uniform) (1-d)
        y = rand(num_data,1);
        data_whiteu = rescale(y,vals_min,vals_max);        
        
    % lorenz attractor  
        [X Y Z] = lorenz(28, 10, 8/3,[0 1 1.05],[0 190],0.000001);
        % [X Y Z] = lorenz(28, 10, 8/3,[0 1 1.05],[0 50],0.000001);
        % function [x,y,z] = lorenz(rho, sigma, beta, initV, T, eps)
        %       X, Y, Z - output vectors of the strange attactor trajectories
        %       RHO     - Rayleigh number
        %       SIGMA   - Prandtl number
        %       BETA    - parameter
        %       INITV   - initial point
        %       T       - time interval
        %       EPS     - ode solver precision
        X = X(1:num_data);
        data_lorenz = rescale(X,vals_min,vals_max);

% observed data
    num_data = num_data_real;
    load c_u_curve_application_EhretDey2022.mat
    y = wet_dataset(1:num_data,1);
    data_datetime = datetime(y,'convertfrom','juliandate');
    d = '1-Oct-1980 00:00:00';	
    t = datetime(d,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    dts = days(1:num_data);
    data_datetime = t + dts;
    y = wet_dataset(1:num_data,2);
    data_p_wet = rescale(y,vals_min,vals_max); 
    y = wet_dataset(1:num_data,3);
    data_q_wet = rescale(y,vals_min,vals_max); 
    y = snow_dataset(1:num_data,3);
    data_q_snow = rescale(y,vals_min,vals_max);       

%% plot the test data sets

fsize = 12; % font size
lw = 1;   % line width
figure('units','normalized','outerposition',[0 0 0.8 1])

% line
    x = subplot(2,3,1)
    plot(data_line,'LineWidth',lw,'color',rgb('chocolate'))
    xlim([0 300]);
    ylim([0 1]);
    ylabel('Normalized value [0,1]');
    xlabel('Time step [-]')
    title('Line')

% random uniform
    x = subplot(2,3,2)
    plot(data_whiteu(1:800),'LineWidth',lw,'color',rgb('hotpink'))
    set(x,'YTick',zeros(1,0));
    xlim([0 300]);
    ylim([0 1]);
    xlabel('Time step [-]')
    title('Random noise')

% lorenz
    x = subplot(2,3,3)
    plot(data_lorenz(2000:5000),'LineWidth',lw,'color',rgb('darkviolet'))
    set(x,'YTick',zeros(1,0));
    xlim([0 3000]);
    ylim([0 1]);
    xlabel('Time step [-]')
    title('Lorenz attractor')

% p wet
    from_plot = 4748;
    to_plot = 6209;
    data_datetime_plot = data_datetime(from_plot:to_plot);
    
    x = subplot(2,3,4)
    data_p_wet_plot = data_p_wet(from_plot:to_plot);
    data_p_wet_plot = rescale(data_p_wet_plot,vals_min,vals_max); 
    xx = plot(data_datetime_plot,data_p_wet_plot,'LineWidth',lw,'color',rgb('skyblue'))
    ylim([0 1]);
    ylabel('Normalized value [0,1]');
    ax=xx.Parent;
    set(ax, 'XTick', [data_datetime_plot(93) data_datetime_plot(458) data_datetime_plot(823) data_datetime_plot(1189)]); 
    xticklabels({'1994','1995','1996','1997'})
    set(gca,'TickDir','out'); 
    title('Precipitation STR')

% q wet
    x = subplot(2,3,5)
    data_q_wet_plot = data_q_wet(from_plot:to_plot);
    data_q_wet_plot = rescale(data_q_wet_plot,vals_min,vals_max); 
    xx = plot(data_datetime_plot,data_q_wet_plot,'LineWidth',lw,'color',rgb('steelblue'))
    set(x,'YTick',zeros(1,0));
    ylim([0 1]);
    ax=xx.Parent;
    set(ax, 'XTick', [data_datetime_plot(93) data_datetime_plot(458) data_datetime_plot(823) data_datetime_plot(1189)]); 
    xticklabels({'1994','1995','1996','1997'})
    set(gca,'TickDir','out'); 
    title('Streamflow STR')

% q snow
    x = subplot(2,3,6)
    data_q_snow_plot = data_q_snow(from_plot:to_plot);
    data_q_snow_plot = rescale(data_q_snow_plot,vals_min,vals_max); 
    xx = plot(data_datetime_plot,data_q_snow_plot,'LineWidth',lw,'color',rgb('cadetblue'))
    set(x,'YTick',zeros(1,0));
    %xlim([0 1461]);
    ylim([0 1]);
    %xlabel('Time step [d]')
    ax=xx.Parent;
    set(ax, 'XTick', [data_datetime_plot(93) data_datetime_plot(458) data_datetime_plot(823) data_datetime_plot(1189)]); 
    xticklabels({'1994','1995','1996','1997'})
    set(gca,'TickDir','out'); 
    title('Streamflow GR')

annotation('textbox',[0.13 0.88 0.038 0.043],'String',{'(a)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
annotation('textbox',[0.41 0.88 0.038 0.043],'String',{'(b)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
annotation('textbox',[0.69 0.88 0.038 0.043],'String',{'(c)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
annotation('textbox',[0.13 0.405 0.038 0.043],'String',{'(d)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
annotation('textbox',[0.41 0.405 0.038 0.043],'String',{'(e)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
annotation('textbox',[0.69 0.405 0.038 0.043],'String',{'(f)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');

% set(findobj(gcf,'type','axes'),'FontSize',fsize,'FontWeight','Bold', 'LineWidth', 0.5);
set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
% print ('fig_timeseries', '-dpng', '-r600');

%% characterize the systems

% artifical data
    slice_widths = slice_widths_arti;

    % horizontal line (1-d)
    [kld_mean_line,kld_entropy_line,~,~] = f_c_u_curve(data_line, bin_edges, bin_edges_kld, slice_widths);
%     [kld_mean_line,kld_entropy_line,~,~] = f_uncertainty_complexity_fast(data_line, bin_edges, bin_edges_kld, slice_widths);

    % white noise (uniform) (1-d)    
    [kld_mean_whiteu,kld_entropy_whiteu,~,~] = f_c_u_curve(data_whiteu, bin_edges, bin_edges_kld, slice_widths);

    % lorenz system (1-d)
    [kld_mean_lorenz,kld_entropy_lorenz,~,~] = f_c_u_curve(data_lorenz, bin_edges, bin_edges_kld, slice_widths);

% real-world data
    slice_widths = slice_widths_real;

    % p data wet (1-d)
    [kld_mean_p_wet,kld_entropy_p_wet,~,~] = f_c_u_curve(data_p_wet, bin_edges, bin_edges_kld, slice_widths);

    % q data wet (1-d)
    [kld_mean_q_wet,kld_entropy_q_wet,~,~] = f_c_u_curve(data_q_wet, bin_edges, bin_edges_kld, slice_widths);

    % q data snow (1-d)
    [kld_mean_q_snow,kld_entropy_q_snow,~,~] = f_c_u_curve(data_q_snow, bin_edges, bin_edges_kld, slice_widths);

%% calculate mean c and u values

u_line = mean(kld_mean_line); c_line = mean(kld_entropy_line);
u_whiteu = mean(kld_mean_whiteu); c_whiteu = mean(kld_entropy_whiteu);
u_lorenz = mean(kld_mean_lorenz); c_lorenz = mean(kld_entropy_lorenz);
u_p_wet = mean(kld_mean_p_wet); c_p_wet = mean(kld_entropy_p_wet);
u_q_wet = mean(kld_mean_q_wet); c_q_wet = mean(kld_entropy_q_wet);
u_q_snow = mean(kld_mean_q_snow); c_q_snow = mean(kld_entropy_q_snow);

%% plot c-u-curves for artifical curves

fsize = 12;     % font size
lw = 1.5;         % line width
m_size = 30;    % marker size
t_size = 12;
inc = 0.025;     % offset increment

% get total number of bins used to calculate limits
num_bins = prod(cellfun(@length,bin_edges)-1);  % total number of bins for uncertainty
num_bins_kld = prod(cellfun(@length,bin_edges_kld)-1);    % total number of bins for complexity
max_uncertainty = log2(num_bins);
max_complexity = log2(num_bins_kld);

dummy_str = string(slice_widths_arti);
    
figure('units','normalized','outerposition',[0 0 0.8 1])   
hold on

plot(kld_mean_line,kld_entropy_line,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('chocolate'))  
plot(kld_mean_whiteu,kld_entropy_whiteu,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('hotpink'))
plot(kld_mean_lorenz,kld_entropy_lorenz,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('darkviolet'))

plot(u_line,c_line,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('chocolate'))
plot(u_whiteu,c_whiteu,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('hotpink'))
plot(u_lorenz,c_lorenz,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('darkviolet'))  

% show some or all
indx = ([2 3 4]);
text(kld_mean_whiteu(indx)+inc,kld_entropy_whiteu(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('hotpink'));  
indx = (1:1:length(slice_widths_arti));
text(kld_mean_lorenz(indx)+inc,kld_entropy_lorenz(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('darkviolet'));

xline(max_uncertainty,'-k','max. Uncertainty','LabelVerticalAlignment','bottom','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);
yline(max_complexity,'-k','max. Complexity','LabelHorizontalAlignment','left','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);   

xlim([0 max_uncertainty + 0.2])
ylim([0 max_complexity + 0.2])
xlabel('Uncertainty [bit]')
ylabel('Complexity [bit]')

h = legend('Line','Random noise','Lorenz attractor');

set(h,'Location','northeast');
set(h,'Position',[0.81 0.76 0.07 0.11]);
legend boxoff
set(gca,'FontSize',fsize,'FontWeight','bold')
set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
hold off

% print ('fig_c-u-curve_artifical', '-dpng', '-r600');

%% plot c-u-curves for real-world data

fsize = 12;     % font size
lw = 1.5;         % line width
m_size = 30;    % marker size
t_size = 12;
inc = 0.015;     % offset increment

% get total number of bins used to calculate limits
num_bins = prod(cellfun(@length,bin_edges)-1);  % total number of bins for uncertainty
num_bins_kld = prod(cellfun(@length,bin_edges_kld)-1);    % total number of bins for complexity
max_uncertainty = log2(num_bins);
max_complexity = log2(num_bins_kld);

dummy_str = string(slice_widths_real);
    
figure('units','normalized','outerposition',[0 0 0.8 1])   
hold on

plot(kld_mean_p_wet,kld_entropy_p_wet,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('skyblue'))
plot(kld_mean_q_wet,kld_entropy_q_wet,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('steelblue'))
plot(kld_mean_q_snow,kld_entropy_q_snow,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('cadetblue'))

plot(u_p_wet,c_p_wet,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('skyblue'))
plot(u_q_wet,c_q_wet,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('steelblue'))   
plot(u_q_snow,c_q_snow,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('cadetblue'))       

% show some
%indx = (1:1:length(slice_widths_real));
indx = ([1:1:10 12]);
text(kld_mean_p_wet(indx)+inc,kld_entropy_p_wet(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('skyblue'));

indx = ([1 3 12]);
text(kld_mean_q_wet(indx)+inc,kld_entropy_q_wet(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('steelblue'));

indx = ([1:1:10 12]);
text(kld_mean_q_snow(indx)+inc,kld_entropy_q_snow(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('cadetblue'));

yline(max_complexity,'-k','max. Complexity','LabelHorizontalAlignment','left','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);   

xlim([0 1])
ylim([0 max_complexity + 0.2])
xlabel('Uncertainty [bit]')
ylabel('Complexity [bit]')

h = legend('Precipitation STR','Streamflow STR','Streamflow GR');

set(h,'Location','northeast');
set(h,'Position',[0.81 0.76 0.07 0.11]);
legend boxoff
set(gca,'FontSize',fsize,'FontWeight','bold')
set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
hold off

% print ('fig_c-u-curve_real', '-dpng', '-r600');
