% Uwe Ehret, 2022/11/01
% This script creates all plots for 
% Ehret, U., and Dey, P.: Technical note: c-u-curve: A method to analyse, classify and compare dynamical systems by uncertainty and complexity, 
% Hydrol. Earth Syst. Sci. Discuss., 2022, 1-12, 10.5194/hess-2022-16, 2022.
% Required Matlab products: Matlab 9.9

clearvars
clc
close all

%% settings

    nt_arti = 30000;    % for artificial data: number of rows (=time steps) in the test data set 
    nt_real = 12418;    % for real-world data: number of rows (=time steps) in the test data set
    ndim = 1;           % number of colums (=variables) in the test data set
    nens = 1;           % number of ensemble members
    vals_min = 0;       % minimum value in all test data sets
    vals_max = 1;       % maximum value in all test data sets
    nvb = 10;           % number of equal-size bins the value range of the data will be split in to calculate histograms
    neb = 10;           % number of equal-size bins the kld range will be split in to calculate the entropy of the kld distribution

% create edges of value bins
    % - [1,ndim] cell array, with a [1,nvb+1] array of bin edges for each dimension inside
    edges_vals = cell(1,ndim);
    edges_vals{1} = linspace(vals_min,vals_max,nvb+1);

% create edges of entropy bins    
% - [1,1] cell array, with a [1,neb+1] array of bin edges for the 1-d histogram of entropy values
% - the possible value range for entropy values is always [0,log2(total number of value bins for the entire ndim-dimensional space of values)]
%   - for a 1-d data set: [0,log2(nvb)] 
%   - for a 2-d data set: [0,log2(nvb_in_first_dimension * nvb_in_second_dimension)]
%   --> set the smallest bin edge 0, the upper to the upper value of the value range
    edges_entropy = cell(1,1);
    edges_entropy{1} = linspace(0,log2(nvb),neb+1);

% array with all time-slice widhts to be examined (for which uncertainty and complexity are to be calculated)
% - size is [nss,1], with nss being total number of time-slicing schemes to be examined
% - order is ascending, minimum possible value is 1, maximum possible value is nt
    sw_arti = [1 30:10:90 100:50:200 300:100:500 1000 30000];   % for artificial data
    sw_real = [1 7 14 21 30 60 91 182 365 730 6209 12418]';     % for real-world data

%% create test data sets

% artifical data
    nt = nt_arti;

    % horizontal line (1-d)
        data_line = zeros(nt,1) + mean([vals_min vals_max]); 
        
    % white noise (uniform) (1-d)
        y = rand(nt,1);
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
        X = X(1:nt);
        data_lorenz = rescale(X,vals_min,vals_max);

% observed data (1-d)
    nt = nt_real
    load Dataset.mat
    y = wet_dataset(1:nt,1);
    data_datetime = datetime(y,'convertfrom','juliandate');
    d = '1-Oct-1980 00:00:00';	
    t = datetime(d,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    dts = days(1:nt);
    data_datetime = t + dts;
    y = wet_dataset(1:nt,2);
    data_p_wet = rescale(y,vals_min,vals_max); 
    y = wet_dataset(1:nt,3);
    data_q_wet = rescale(y,vals_min,vals_max); 
    y = snow_dataset(1:nt,3);
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
    %xlim([0 1461]);
    ylim([0 1]);
    ylabel('Normalized value [0,1]');
    ax=xx.Parent;
    set(ax, 'XTick', [data_datetime_plot(93) data_datetime_plot(458) data_datetime_plot(823) data_datetime_plot(1189)]); 
    xticklabels({'1994','1995','1996','1997'})
    set(gca,'TickDir','out'); 
    %xlabel('Time step [d]')
    title('Precipitation South Toe River STR')
    
    % q wet
    x = subplot(2,3,5)
    data_q_wet_plot = data_q_wet(from_plot:to_plot);
    data_q_wet_plot = rescale(data_q_wet_plot,vals_min,vals_max); 
    xx = plot(data_datetime_plot,data_q_wet_plot,'LineWidth',lw,'color',rgb('steelblue'))
    set(x,'YTick',zeros(1,0));
    %xlim([0 1461]);
    ylim([0 1]);
    %xlabel('Time step [d]')
    ax=xx.Parent;
    set(ax, 'XTick', [data_datetime_plot(93) data_datetime_plot(458) data_datetime_plot(823) data_datetime_plot(1189)]); 
    xticklabels({'1994','1995','1996','1997'})
    set(gca,'TickDir','out'); 
    title('Streamflow South Toe River STR')
    
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
    title('Streamflow Green River GR')
    
    annotation('textbox',[0.13 0.88 0.038 0.043],'String',{'(a)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.41 0.88 0.038 0.043],'String',{'(b)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.69 0.88 0.038 0.043],'String',{'(c)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.13 0.405 0.038 0.043],'String',{'(d)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.41 0.405 0.038 0.043],'String',{'(e)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.69 0.405 0.038 0.043],'String',{'(f)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    
    % set(findobj(gcf,'type','axes'),'FontSize',fsize,'FontWeight','Bold', 'LineWidth', 0.5);
    set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
    print ('Fig01', '-dpng', '-r600');

%% characterize the systems

% artifical data
    slice_widths = sw_arti;

    % horizontal line (1-d)
    [uncs_line,comps_line,~,~] = f_c_u_curve(data_line, edges_vals, edges_entropy, slice_widths,0);

    % white noise (uniform) (1-d)    
    [uncs_whiteu,comps_whiteu,~,~] = f_c_u_curve(data_whiteu, edges_vals, edges_entropy, slice_widths,0);

    % lorenz system (1-d)
    [uncs_lorenz,comps_lorenz,~,~] = f_c_u_curve(data_lorenz, edges_vals, edges_entropy, slice_widths,0);

% real-world data
    slice_widths = sw_real;

    % p data wet (1-d)
    [uncs_p_wet,comps_p_wet,~,~] = f_c_u_curve(data_p_wet, edges_vals, edges_entropy, slice_widths,0);

    % q data wet (1-d)
    [uncs_q_wet,comps_q_wet,~,~] = f_c_u_curve(data_q_wet, edges_vals, edges_entropy, slice_widths,0);

    % q data snow (1-d)
    [uncs_q_snow,comps_q_snow,~,~] =f_c_u_curve(data_q_snow, edges_vals, edges_entropy, slice_widths,0);

%% calculate mean values of uncertainty and complexity

    u_line = mean(uncs_line); c_line = mean(comps_line);
    u_whiteu = mean(uncs_whiteu); c_whiteu = mean(comps_whiteu);
    u_lorenz = mean(uncs_lorenz); c_lorenz = mean(comps_lorenz);
    u_p_wet = mean(uncs_p_wet); c_p_wet = mean(comps_p_wet);
    u_q_wet = mean(uncs_q_wet); c_q_wet = mean(comps_q_wet);
    u_q_snow = mean(uncs_q_snow); c_q_snow = mean(comps_q_snow);

%% calculate upper bounds of complexity as a function of uncertainty

    states = linspace(0,log2(nvb),neb); % discrete (binned) values the entropy distribution can take
    means = (0:0.01:log2(nvb));    % array of mean values, covering the uncertainty range [0,log2(nvb)]
    Hmax = f_maxEnt_known_mean(states,means);

%% plot Fig. 02: c-u-curves for artifical curves

    fsize = 12;     % font size
    lw = 1.5;       % line width
    m_size = 30;    % marker size
    t_size = 12;
    inc = 0.025;    % offset increment

% get total number of bins used to calculate limits
    nvb = prod(cellfun(@length,edges_vals)-1);      % total number of bins for uncertainty
    neb = prod(cellfun(@length,edges_entropy)-1);   % total number of bins for complexity
    max_uncertainty = log2(nvb);
    max_complexity = log2(neb);

    dummy_str = string(sw_arti);

% plot    
    figure('units','normalized','outerposition',[0 0 0.8 1])   
    hold on
    
    % plot c-u-curves of artifical data
    plot(uncs_line,comps_line,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('chocolate'))  
    plot(uncs_whiteu,comps_whiteu,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('hotpink'))
    plot(uncs_lorenz,comps_lorenz,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('darkviolet'))
  
    % plot c-u-curves of real-world data    
    plot(uncs_p_wet,comps_p_wet,'Linestyle','-','Marker','none','Markersize',m_size,'LineWidth',lw/2,'color',rgb('skyblue'))
    plot(uncs_q_wet,comps_q_wet,'Linestyle','-','Marker','none','Markersize',m_size,'LineWidth',lw/2,'color',rgb('steelblue'))
    plot(uncs_q_snow,comps_q_snow,'Linestyle','-','Marker','none','Markersize',m_size,'LineWidth',lw/2,'color',rgb('cadetblue'))
    
    % plot mean values of uncertainty and complexity
    plot(u_line,c_line,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('chocolate'))
    plot(u_whiteu,c_whiteu,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('hotpink'))
    plot(u_lorenz,c_lorenz,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('darkviolet'))  
    
    % plot labels
    indx = ([2 3 4]);
    text(uncs_whiteu(indx)+inc,comps_whiteu(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('hotpink'));  
    indx = (1:1:length(sw_arti));
    text(uncs_lorenz(indx)+inc,comps_lorenz(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('darkviolet'));
    
    % plot upper bounds for uncertainty and complexity
    xline(max_uncertainty,'-k','max. Uncertainty','LabelVerticalAlignment','bottom','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);
    yline(max_complexity,'-k','max. Complexity','LabelHorizontalAlignment','left','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);   
    
    % plot upper bound for complexity as plot(means,Hmax) as a function of uncertainty
    plot(means,Hmax,'Linestyle','-','Marker','none','Markersize',m_size,'LineWidth',lw,'color',rgb('black'))
    annotation('textbox',[0.15 0.807 0.119 0.043],'String',{'max. Complexity'},'EdgeColor','none','FontSize',fsize,'FitBoxToText','on');

    xlim([0 max_uncertainty + 0.2])
    ylim([0 max_complexity + 0.2])
    xlabel('Uncertainty [bit]')
    ylabel('Complexity [bit]')
        
    h = legend('Line','Random noise','Lorenz attractor','Precipitation STR','Streamflow STR','Streamflow GR');
    set(h,'Location','northeast');
    set(h,'Position',[0.81 0.76 0.07 0.11]);
    legend boxoff
    
    set(gca,'FontSize',fsize,'FontWeight','bold')
    set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
    hold off
    
    print ('Fig02', '-dpng', '-r600');

%% plot Fig. 03: c-u-curves for real-world data
    
    fsize = 12;     % font size
    lw = 1.5;       % line width
    m_size = 30;    % marker size
    t_size = 12;
    inc = 0.015;    % offset increment

% get total number of bins used to calculate limits
    nvb = prod(cellfun(@length,edges_vals)-1);      % total number of bins for uncertainty
    neb = prod(cellfun(@length,edges_entropy)-1);   % total number of bins for complexity
    max_uncertainty = log2(nvb);
    max_complexity = log2(neb);

    dummy_str = string(sw_real);

% plot    
    figure('units','normalized','outerposition',[0 0 0.8 1])   
    hold on

    % plot c-u-curves of real-world data   
    plot(uncs_p_wet,comps_p_wet,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('skyblue'))
    plot(uncs_q_wet,comps_q_wet,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('steelblue'))
    plot(uncs_q_snow,comps_q_snow,'Linestyle','-','Marker','.','Markersize',m_size,'LineWidth',lw,'color',rgb('cadetblue'))
    
    % plot mean values of uncertainty and complexity 
    plot(u_p_wet,c_p_wet,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('skyblue'))
    plot(u_q_wet,c_q_wet,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('steelblue'))   
    plot(u_q_snow,c_q_snow,'Linestyle','none','Marker','p','Markersize',m_size,'MarkerFaceColor',rgb('cadetblue'))       
        
    % plot labels
    indx = ([1:1:10 12]);
    text(uncs_p_wet(indx)+inc,comps_p_wet(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('skyblue'));
    indx = ([1 3 12]);
    text(uncs_q_wet(indx)+inc,comps_q_wet(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('steelblue'));
    indx = ([1:1:10 12]);
    text(uncs_q_snow(indx)+inc,comps_q_snow(indx)+inc,dummy_str(indx),'FontSize',t_size,'FontWeight','Bold','color',rgb('cadetblue'));
    
    % plot upper bounds for uncertainty and complexity
    %xline(max_uncertainty,'-k','max. Uncertainty','LabelVerticalAlignment','bottom','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);
    yline(max_complexity,'-k','max. Complexity','LabelHorizontalAlignment','left','LineWidth',lw,'color',rgb('black'),'FontSize',fsize);   
    
    % plot upper bound for complexity as a function of uncertainty
    plot(means,Hmax,'Linestyle','-','Marker','none','Markersize',m_size,'LineWidth',lw,'color',rgb('black'))
    annotation('textbox',[0.8 0.84 0.1 0.043],'String',{'max. Complexity'},'EdgeColor','none','FontSize',fsize,'FitBoxToText','on');

    xlim([0 1])
    ylim([0 max_complexity + 0.2])
    xlabel('Uncertainty [bit]')
    ylabel('Complexity [bit]')
    
    h = legend('Precipitation STR','Streamflow STR','Streamflow GR');
    
    set(h,'Location','northwest');
    set(h,'Position',[0.15 0.76 0.07 0.11]);
    legend boxoff
    set(gca,'FontSize',fsize,'FontWeight','bold')
    set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
    hold off
    
    print ('Fig03', '-dpng', '-r600');

%% plot Fig A1: time series and distribution of values

    fsize = 12; % font size
    lw = 2;     % line width
    figure('units','normalized','outerposition',[0 0 0.8 1])
      
    % slice 1 (H=0 = minH)
        from = 1; to = 60;
        p_uncs = histcounts(data_q_snow(from:to),edges_vals{1},'Normalization', 'probability');
        x = subplot(2,3,1)
        xx = plot(data_datetime(from:to),data_q_snow(from:to),'LineWidth',lw,'color',rgb('cadetblue'))
        ylim([0 1]);
        ylabel('Normalized value [0,1]');
        
        x = subplot(2,3,4)
        bar(p_uncs,0.85,'FaceColor',rgb('cadetblue'))
        ylim([0 1]);
        xlabel('Normalized value (binned) [0,1]')
        ylabel('Probability [-]');
        set(gca, 'XTick', (1:11) - 0.5);
        set(gca, 'XTickLabel', {'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'});
        title('Entropy = 0 bit')

    % slice 62 (H=0.61 ~ Hmean=0.6027)
        from = 3661; to = 3720;
        p_uncs = histcounts(data_q_snow(from:to),edges_vals{1},'Normalization', 'probability');
        x = subplot(2,3,2)
        xx = plot(data_datetime(from:to),data_q_snow(from:to),'LineWidth',lw,'color',rgb('cadetblue'))
        ylim([0 1]);

        x = subplot(2,3,5)
        bar(p_uncs,0.85,'FaceColor',rgb('cadetblue'))
        ylim([0 1]);
        xlabel('Normalized value (binned) [0,1]')
        set(gca, 'XTick', (1:11) - 0.5);
        set(gca, 'XTickLabel', {'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'});
        title('Entropy = 0.61 bit')

    % slice 168 (H=2.27 = maxH)
        from = 10021; to = 10080;
        p_uncs = histcounts(data_q_snow(from:to),edges_vals{1},'Normalization', 'probability');
        x = subplot(2,3,3)
        xx = plot(data_datetime(from:to),data_q_snow(from:to),'LineWidth',lw,'color',rgb('cadetblue'))
        ylim([0 1]);

        x = subplot(2,3,6)
        bar(p_uncs,0.85,'FaceColor',rgb('cadetblue'))
        ylim([0 1]);
        xlabel('Normalized value (binned) [0,1]')
        set(gca, 'XTick', (1:11) - 0.5);
        set(gca, 'XTickLabel', {'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'});
        title('Entropy = 2.27 bit')

    annotation('textbox',[0.13 0.88 0.038 0.043],'String',{'(a)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.41 0.88 0.038 0.043],'String',{'(b)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.69 0.88 0.038 0.043],'String',{'(c)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.13 0.405 0.038 0.043],'String',{'(d)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.41 0.405 0.038 0.043],'String',{'(e)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    annotation('textbox',[0.69 0.405 0.038 0.043],'String',{'(f)'},'FontSize',fsize,'FontWeight','bold','EdgeColor','none');
    
    % set(findobj(gcf,'type','axes'),'FontSize',fsize,'FontWeight','Bold', 'LineWidth', 0.5);
    set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
    print ('FigA1', '-dpng', '-r600');

%% plot Fig A2: distribution of entropies  

    sw = 60; % slice width (60 days)

    fsize = 12; % font size
    lw = 2;   % line width
    figure('units','normalized','outerposition',[0 0 0.8 1])
    [unc,comp,ns,all_uncs] = f_c_u_curve(data_q_snow, edges_vals, edges_entropy, sw,0);
    p_uncs = histcounts(all_uncs{1,1},edges_entropy{1},'Normalization', 'probability');
    bar(p_uncs,0.85,'FaceColor',rgb('cadetblue'))
    ylim([0 1]);
    ylabel('Probability [-]');
    xlabel('Uncertainty (binned) [bit]')
    set(gca, 'XTick', (1:11) - 0.5);
    set(gca, 'XTickLabel', {'0' '0.33' '0.66' '0.99' '1.33' '1.66' '1.99' '2.32' '2.66' '2.99' 'log(10)=3.32'});
    title('Entropy = 2.33 bit')

    set(gca,'LooseInset',get(gca,'TightInset')); % erase unnnecesary outside whitespace
    print ('FigA2', '-dpng', '-r600');
