function Growth_data_analysis(data_set, with_pie, rates_normalized, filter_suppressors)
% This function loads processed growth curve data (derived from raw plate
% reader data) a turns this into growth rate estimated for strain/medium
% combinations. Relevant combinations are plotted.
% This file uses a overdispersion correction (Birge ratio) from [1].
%
% Inputs:
% data_set          Either 1 (for data main text) or 2 (data supplementary
%                   information)
% with_pie          true/false whether marker symbols in fitness plot
%                   should contain pie diagrams on growth success
% rates_normalized  true/false whether fitness in plot should be normalized
%                   to WT fitness in the same medium
% filter_suppressors true/false whether growth in certain measurements
%                   should be discarded when emergence of a suppressor
%                   mutant is suspected
%
% Outputs:
% A .mat file with processed data and a fitness plot with relevant strains
%
% Bibliography:
%
% [1] Bodnar O, Elster C. On the adjustment of inconsistent data using the
% Birge ratio. Metrologia. 2014 Oct 1;51(5):516–21.
%

switch data_set
    case 1
        List_GUI_output_files   = { './Processed plate reader outputs/2017_12_22.mat'; ...
                                    './Processed plate reader outputs/2017_12_25.mat'; ...
                                    './Processed plate reader outputs/2018_01_26.mat'; ...
                                    './Processed plate reader outputs/2018_01_28.mat'};
    case 2
        List_GUI_output_files   = { './Processed plate reader outputs/13022021.mat'; ...
                                    './Processed plate reader outputs/21022021.mat'; ...
                                    './Processed plate reader outputs/30042021.mat'; ...
                                    './Processed plate reader outputs/10052021_1.mat'; ...
                                    './Processed plate reader outputs/10052021_2.mat'; ...
                                    './Processed plate reader outputs/21052021_1.mat'; ...
                                    './Processed plate reader outputs/21052021_2.mat'; ...
                                    './Processed plate reader outputs/06082021.mat'; ...
                                    './Processed plate reader outputs/24092021.mat'};
end

%% Loading data into one variable
% Load into cell
Outputs_cell                        = cell(numel(List_GUI_output_files), 1);
for file_no = 1 : numel(List_GUI_output_files)
    Output_temp                     = load( List_GUI_output_files{file_no, 1} );
    if ~isempty(regexp(List_GUI_output_files{file_no, 1}, '2018_01_28', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(1 : 88);
        % Remove contaminated last column from further analysis
    elseif ~isempty(regexp(List_GUI_output_files{file_no, 1}, '13022021', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(setdiff(1 : 96, ...
            [18 : 24 : 72 19 : 24 : end 30 : 24 : end 31 : 24 : end]));
        % Last line to remove unreliable YWKD073
    elseif ~isempty(regexp(List_GUI_output_files{file_no, 1}, '21022021', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(setdiff(1 : 96, [41 : 48 73 : 88]));
        % Last line to remove mislabelled YWKD065b and unreliable YWKD073
    elseif ~isempty(regexp(List_GUI_output_files{file_no, 1}, '30042021', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(setdiff(1 : 96, 65 : 80));
        % Last line to remove unreliable YWKD073
    elseif ~isempty(regexp(List_GUI_output_files{file_no, 1}, '(10052021|21052021)', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(setdiff(1 : 96, 57 : 72));
        % Last line to remove unreliable YWKD073
    elseif ~isempty(regexp(List_GUI_output_files{file_no, 1}, '(06082021|24092021)', 'once'))
        Outputs_cell{file_no, 1}    = Output_temp.Wells(setdiff(1 : 96, 81 : 88));
        % Remove mislabelled YWKD065b
    else
        Outputs_cell{file_no, 1}    = Output_temp.Wells(1 : end);
    end
    [Outputs_cell{file_no, 1}.Experiment]   = deal(List_GUI_output_files{file_no});
    % Keep track which fields are common to all data files (to allow
    % combining all data later)
    if file_no == 1
        All_fields                  = fieldnames(Outputs_cell{file_no, 1});
    else
        All_fields                  = intersect(All_fields, fieldnames(Outputs_cell{file_no, 1}));
    end
end
% Trim fields that are not always present for some reason
for file_no = 1 : numel(List_GUI_output_files)
    Current_fields              = fieldnames(Outputs_cell{file_no, 1});
    Removed_fields              = setdiff(Current_fields, All_fields);
    Outputs_cell{file_no, 1}    = rmfield(Outputs_cell{file_no, 1}, Removed_fields);
end
% Vectorize
Outputs_vec = [Outputs_cell{:}]';
% Keep relevant variables
clearvars -except Outputs_vec filter_suppressors data_set with_pie rates_normalized

%% Get convenient data format
% Keep only the relevant fields
fields_to_keep  = {'t_doubl', 't_doubl_err', 'Data', 'Genotype', 'Medium', 'base', 'end_OD', 'Size_Fit_Window', 'Strain', 'Experiment'};
Outputs_vec2    = rmfield(Outputs_vec, setdiff(fieldnames(Outputs_vec), fields_to_keep));
Outputs_vec2    = orderfields(Outputs_vec2, fields_to_keep);
Outputs_vec2    = struct2cell(Outputs_vec2)';

% Keep track which field is in which index
t_doubl_col     = find(strcmp('t_doubl', fields_to_keep));
t_err_col       = find(strcmp('t_doubl_err', fields_to_keep));
window_col      = find(strcmp('Size_Fit_Window', fields_to_keep));
data_col        = find(strcmp('Data', fields_to_keep));
base_col        = find(strcmp('base', fields_to_keep));
end_col         = find(strcmp('end_OD', fields_to_keep));
gen_col         = find(strcmp('Genotype', fields_to_keep));
med_col         = find(strcmp('Medium', fields_to_keep));
strain_col      = find(strcmp('Strain', fields_to_keep));
exp_col         = find(strcmp('Experiment', fields_to_keep));

% Convert doubling times back to growth rates
for i = 1 : size(Outputs_vec2, 1)
    Outputs_vec2{i, 1} = 1 / Outputs_vec2{i, t_doubl_col};
    Outputs_vec2{i, 2} = - diff(1 ./ Outputs_vec2{i, t_err_col}) / (2 * tinv(0.975, Outputs_vec2{i, window_col} - 2));
end

% Delete empty wells
Outputs_vec2    = Outputs_vec2( ~strcmp(Outputs_vec2(:, gen_col), 'None'), :);

% Trim medium text
for i = 1 : size(Outputs_vec2, 1)
    Outputs_vec2(i, med_col) = strcat(regexp(Outputs_vec2{i, med_col}, '(?<=Raf \+ )[\w.]{1,5}', 'match'), '% Gal');
end

% Replace inconsistent naming of WT
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'WT'), gen_col} ]= deal('WT CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'WT + CDC42'), gen_col} ]= deal('WT CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'WT + Gal1-sfGFP-CDC42'), gen_col} ]= deal('WT Gal1-sfGFP-CDC42');
% Replace inconsistent naming of dbem1
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + CDC42'), gen_col} ]= deal('dbem1 CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + Gal1-CDC42'), gen_col} ]= deal('dbem1 Gal1-CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + Gal1-sfGFP-CDC42'), gen_col} ]= deal('dbem1 Gal1-sfGFP-CDC42');
% Replace inconsistent naming of dbem1dbem3
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + dbem3 + CDC42'), gen_col} ]= deal('dbem1 dbem3 CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + dbem3 + Gal1-CDC42'), gen_col} ]= deal('dbem1 dbem3 Gal1-CDC42');
[Outputs_vec2{ strcmp(Outputs_vec2(:,gen_col), 'dbem1 + dbem3 + Gal1-sfGFP-CDC42'), gen_col} ]= deal('dbem1 dbem3 Gal1-sfGFP-CDC42');

%% Quick check suspected suppressor growth
% Fit tanh function on data, and extrapolate what the implied initial OD
% should have been. If this is unreasonably low, mark this well as
% suppressor
[Outputs_vec2{:, end + 1}] = deal(false(1));
for i = 1 : size(Outputs_vec2, 1)
    x                       = Outputs_vec2{i, data_col}(1, :) / 60;
    y                       = Outputs_vec2{i, data_col}(2, :)' - Outputs_vec2{i, base_col};
    x2                      = x(~isnan(y))';
    y2                      = y(~isnan(y));
    
    if ~isempty(x2) || numel(x2)>3
        p0_a                = max(y2) / 2;
        p0_b                = 2 * atanh(0.8) / (x2(find(y2 >= min(y2) + 0.9 *(max(y2) - min(y2)), 1)) - ...
                                x2(find(y2 >= min(y2) + 0.1 *(max(y2) - min(y2)), 1)));
        if isinf(p0_b)
            p0_b            = 0.1;
        end
        p0_c                = x2(find(y2 > p0_a, 1));
        
        obj                 = fit(x2, y2, fittype('a * (tanh(b * (x - c)) + 1)', ...
                                'options', fitoptions('Method', 'NonLinearLeastSquares', 'StartPoint', [p0_a; p0_b; p0_c], ...
                                'Lower', [0 0 -max(x2)/p0_b * 10], 'Upper', [10 100 max(x2)/p0_b * 10])));
        
        Outputs_vec2{i, end}= feval(obj, 0) < 1e-3 && ~isnan(Outputs_vec2{i, t_doubl_col});
    end
end

%% Further rejection of growths (if OD rise < 0.01)
% Sometimes growth stochadtically possible, but colony ultimately dies out before a significant OD rise.
% Analysis of the end ODs of wells labelled as growth/non-growth suggests
% this threshold is around 0.01-0.02 OD rise (clear cut-off there)

OD_rise_threshold       = 0.01;
OD_rise_growth          = zeros(size(Outputs_vec2, 1), 2);
for i = 1 : size(Outputs_vec2, 1)
    OD_rise_growth(i, 1 : 2) = [Outputs_vec2{i, end_col} - Outputs_vec2{i, base_col} ~isnan(Outputs_vec2{i, t_doubl_col})];
end
[Outputs_vec2{ OD_rise_growth(:,1) < OD_rise_threshold , 1 }] = deal(NaN);
[Outputs_vec2{ OD_rise_growth(:,1) < OD_rise_threshold , 2 }] = deal(NaN);

if filter_suppressors
    Outputs_vec_full    = Outputs_vec2;
    Outputs_vec2        = Outputs_vec2( ~[Outputs_vec2{:, end}] , 1 : end - 1);
end

%% Combine replicates
Genotypes                       = table2cell(unique(cell2table(Outputs_vec2(:, gen_col))));
Media                           = table2cell(unique(cell2table(Outputs_vec2(:, med_col))));

if data_set == 2
    % Keep genotypes relevant for the plot later (with BEM1 BEM3 background)
    ind_gen_keep                = cellfun(@(x) ~isempty(regexp(x, 'BEM1 BEM3', 'match')), Genotypes);
    Genotypes                   = Genotypes(ind_gen_keep);
end

Summary_cell                    = cell(numel(Genotypes) + 1, numel(Media) + 1 );
Summary_cell(2 : end, 1)        = Genotypes;
Summary_cell(1, 2 : end)        = Media';

for genotype_no = 1 : numel(Genotypes)
    ind_this_genotype           = strcmp(Outputs_vec2(:, gen_col), Genotypes(genotype_no, 1));
    
    for media_no = 1 : numel(Media)
        ind_this_medium         = strcmp(Outputs_vec2(:, med_col), Media(media_no, 1));
        if isempty(Outputs_vec2(ind_this_genotype & ind_this_medium, 1 : 2))
            Summary_cell{genotype_no + 1, media_no + 1} = NaN(1, 5);
        else
            rates                   = cell2mat(Outputs_vec2(ind_this_genotype & ind_this_medium, 1 : 2));
            num_replicates          = size(rates, 1);
            num_growths             = nnz(~isnan(rates(:, 1)));
            weights                 = 1 ./ rates(:, 2) .^ 2;
            rates_comb              = nansum(weights .* rates(:, 1)) / nansum(weights);
            rates_err_comb          = sqrt(1 / nansum(weights));
            if rates_err_comb == Inf
                rates_comb          = 0;
            end
            chisq_red               = 1/(num_growths - 1)*nansum((rates(:, 1)-rates_comb).^2./rates(:, 2).^2);
            Summary_cell{genotype_no + 1, media_no + 1} = [num_growths num_replicates rates_comb rates_err_comb chisq_red];
        end
        Summary_cell{genotype_no + 1, media_no + 1}(6) = numel(unique(Outputs_vec2(ind_this_genotype & ind_this_medium, strain_col)));
        Summary_cell{genotype_no + 1, media_no + 1}(7) = numel(unique(Outputs_vec2(ind_this_genotype & ind_this_medium, exp_col)));
    end
end
chisq_red_mat           = cellfun(@(z) z(5), Summary_cell(2 : end, 2 : end));
num_obs_WT_mat          = cellfun(@(z) z(1), Summary_cell(2, 2 : end));
ind_media_discard       = (num_obs_WT_mat <= 3); % Discard these media as not enough WT control replicates
over_dispersion_corr    = sqrt((num_obs_WT_mat(~ind_media_discard) - 1) / (num_obs_WT_mat(~ind_media_discard) - 3)).* ...
                            sqrt(mean(chisq_red_mat(1, ~ind_media_discard)));
% Overdispersion correction factor (see [1]), correct errors with this
ind_media_discard       = false(size(ind_media_discard));
% only discard strains from media over calculation in overdispersion correction
Summary_cell            = Summary_cell(:, [true(1) ~ind_media_discard]);
Media                   = Media(~ind_media_discard');

for genotype_no = 1 : numel(Genotypes)
    for media_no = 1 : numel(Media)
        Summary_cell{genotype_no + 1, media_no + 1}(4) = over_dispersion_corr * Summary_cell{genotype_no + 1, media_no + 1}(4);
    end
end
for i = 1 : size(Outputs_vec2, 1)
    Outputs_vec2{i, 2} = over_dispersion_corr * Outputs_vec2{i, 2};
end

% Keep also the full data set as a reference, including suppressors
if filter_suppressors
    Genotypes_full                      = table2cell(unique(cell2table(Outputs_vec_full(:, gen_col))));
    Media_full                          = table2cell(unique(cell2table(Outputs_vec_full(:, med_col))));

    Summary_cell_full                   = cell(numel(Genotypes_full) + 1, numel(Media_full) + 1 );
    Summary_cell_full(2 : end, 1)       = Genotypes_full;
    Summary_cell_full(1, 2 : end)       = Media_full';

    for genotype_no = 1 : numel(Genotypes_full)
        ind_this_genotype               = strcmp(Outputs_vec_full(:, gen_col), Genotypes_full(genotype_no, 1));
        for media_no = 1 : numel(Media_full)
            ind_this_medium             = strcmp(Outputs_vec_full(:, med_col), Media_full(media_no, 1));
            if isempty(Outputs_vec_full(ind_this_genotype & ind_this_medium, 1 : 2))
                Summary_cell_full{genotype_no + 1, media_no + 1} = NaN(1, 2);
            else
                rates                   = cell2mat(Outputs_vec_full(ind_this_genotype & ind_this_medium, 1 : 2));
                num_replicates          = size(rates, 1);
                num_growths             = nnz(~isnan(rates(:, 1)));
                Summary_cell_full{genotype_no + 1, media_no + 1} = [num_growths num_replicates];
            end
        end
    end
end

clearvars -except Outputs_vec2 Summary_cell* Genotypes Media *_col filter_suppressors data_set with_pie rates_normalized

%% Bayesian distributions for dispersion corrected growth rates
rng default  % For reproducibility
r_min           = 0;
r_max           = 1/70; % Not faster than 70 min. doubling time seems likely prior to measurement
delta           = Summary_cell{2, 2}(4)/2;
num_samples     = 5e4;
Samples         = cell(numel(Genotypes) + 1, numel(Media) + 1 );
Summary_cell2   = cell(numel(Genotypes) + 1, numel(Media) + 1 );
Summary_cell3   = cell(numel(Genotypes) + 1, numel(Media) + 1 );
n_sim           = 1e4;

for genotype_no = 1 : numel(Genotypes)
    ind_this_genotype           = strcmp(Outputs_vec2(:, gen_col), Genotypes(genotype_no, 1));  
    for media_no = 1 : numel(Media)
        fprintf('Genotype %0.f of %0.f, Media %0.f of %0.f \n', genotype_no , numel(Genotypes), media_no , numel(Media))
        ind_this_medium         = strcmp(Outputs_vec2(:, med_col), Media(media_no, 1));
        Samples{genotype_no + 1, media_no + 1} = NaN(1, 5);
        Summary_cell2{genotype_no + 1, media_no + 1} = NaN(1, 5);
        if isempty(Outputs_vec2(ind_this_genotype & ind_this_medium, [t_doubl_col t_err_col]))
        else
            rates           = cell2mat(Outputs_vec2(ind_this_genotype & ind_this_medium, [t_doubl_col t_err_col]));
            dof             = cell2mat(Outputs_vec2(ind_this_genotype & ind_this_medium, window_col)) - 2;
            dof             = dof(~isnan(rates(:,1)), :);
            rates           = rates(~isnan(rates(:,1)), :);
            if isempty(rates), continue, end
            likelihood      = @(mu) prod(tpdf((mu - rates(:,1)) ./ rates(:,2), dof));
            prior           = @(mu) unifpdf(mu, r_min, r_max);
            lik_prior_dist  = @(mu) likelihood(mu) * prior(mu);
            random_num_func = @(x) x + rand * 2 * delta - delta;
            start_dist      = @(x,y) normpdf(x, Summary_cell{genotype_no + 1, media_no + 1}(3), ...
                                    10 * Summary_cell{genotype_no + 1, media_no + 1}(4));
            x               = mhsample(Summary_cell{genotype_no + 1, media_no + 1}(3), num_samples, ...
                                'pdf', lik_prior_dist, 'proppdf', start_dist, 'proprnd', random_num_func, 'burnin', 1e3);
            Samples{genotype_no + 1, media_no + 1} = x;
            Summary_cell2{genotype_no + 1, media_no + 1} = [quantile(x, 0.025) mean(x) quantile(x, 0.975) ...
                                                            quantile(x, 0.16) mean(x) quantile(x, 0.84)];
            
            if genotype_no > 1 % relative fitness numbers relative to genotype 1
                relative_x = x(randi(num_samples, n_sim, 1)) ./ Samples{2, media_no + 1}(randi(num_samples, n_sim, 1));
                Summary_cell3{genotype_no + 1, media_no + 1}= [quantile(relative_x, 0.025) mean(relative_x) quantile(relative_x, 0.975) ...
                                                                quantile(relative_x, 0.16) mean(relative_x) quantile(relative_x, 0.84) ];
            else % genotype 1 as fitness reference
                Summary_cell3{genotype_no + 1, media_no + 1} = ...
                    Summary_cell2{genotype_no + 1, media_no + 1} / Summary_cell2{genotype_no + 1, media_no + 1}(2);
            end        
        end
    end
end
clearvars -except Genotypes Media Summary_cell* Samples n_MC filter_suppressors data_set with_pie rates_normalized

if filter_suppressors
    if data_set == 1
        save('Main_data_without_suppressors.mat')
    else
        save('Supplementary_data_without_suppressors.mat')
    end
else
    if data_set == 1
        save('Main_data.mat')
    else
        save('Supplementary_data.mat')
    end
end

%% Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General plot properties %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

font_size                       = 9;
font_name                       = 'Arial';
line_width                      = 1; % Line width
if data_set == 1
    target_width                = 5.5; %inch
    target_height               = 3.9; %inch
    ax_position                 = [0.7 0.4 target_width - 0.8 target_height - 0.5];
else
    target_width                = 3.5; %inch
    target_height               = 2.5; %inch
    ax_position                     = [0.5 0.5 target_width - 0.6 target_height - 0.6];
end
target_res                      = 300; %dpi
mark_size                       = 8; % Marker size
alpha_patch                     = 0.05;
colors_tot                      = lines(4);
color_no_growth                 = [1 1 1];
x_shift                         = 0.2;
r_pie                           = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for plotting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of genotypes of interest
switch data_set
    case 1
        Selected_genotypes      = {'WT CDC42'; 'WT Gal1-sfGFP-CDC42'; 'dbem1 Gal1-sfGFP-CDC42'; 'dbem1 dbem3 Gal1-sfGFP-CDC42'};
        le_names                = {'{\it{CDC42}} {\it{BEM1}} {\it{BEM3}}'; ...
                                    '{\it{pGAL1-CDC42-sfGFP^{SW}}} {\it{BEM1}} {\it{BEM3}}'; ...
                                    '{\it{pGAL1-CDC42-sfGFP^{SW}}} \Delta{\it{bem1}} {\it{BEM3}}'; ...
                                    '{\it{pGAL1-CDC42}} \Delta{\it{bem1}} \Delta{\it{bem3}}'};
    case 2
        Selected_genotypes      = {'BEM1 BEM3 CDC42'; 'BEM1 BEM3 GAL1-CDC42'; 'BEM1 BEM3 GAL1-sfGFP-CDC42'};
        le_names                = {'{\it{CDC42}}'; '{\it{pGAL1-CDC42}}'; '{\it{pGAL1-sfGFP-CDC42^{SW}}}'};
end

ind_sel_genotypes               = arrayfun(@(x) find(strcmp(Genotypes, x)), Selected_genotypes);
Summary_cell                    = Summary_cell( 1 + ind_sel_genotypes, 2 : end);
Summary_cell2                   = Summary_cell2(1 + ind_sel_genotypes, 2 : end);
Summary_cell3                   = Summary_cell3(1 + ind_sel_genotypes, 2 : end);
Samples                         = Samples(1 + ind_sel_genotypes, 2 : end);

% Media axis
x_Gal_cell                      = cellfun(@(z) regexp(z, '[\d.]{1,5}', 'match'), Media)';
x_plot_tick                     = ( 1 : numel(Media) )';

% Pie marker data
num_growth                      = cellfun(@(z) z(1), Summary_cell);
num_exp                         = cellfun(@(z) z(2), Summary_cell);

% Fitness axis
errcellfun                      = @(x, varargin) NaN;
if rates_normalized
    rel_fit                     = cellfun(@(z) z(2), Summary_cell3, 'ErrorHandler', errcellfun)';
    rel_fit_lb                  = rel_fit - cellfun(@(z) z(4), Summary_cell3, 'ErrorHandler', errcellfun)';
    rel_fit_ub                  = cellfun(@(z) z(6), Summary_cell3, 'ErrorHandler', errcellfun)' - rel_fit;
else
    rel_fit                     = cellfun(@(z) z(2), Summary_cell2, 'ErrorHandler', errcellfun)';
    rel_fit_lb                  = rel_fit - cellfun(@(z) z(4), Summary_cell2, 'ErrorHandler', errcellfun)';
    rel_fit_ub                  = cellfun(@(z) z(6), Summary_cell2, 'ErrorHandler', errcellfun)' - rel_fit;
end
rel_fit(   num_growth' == 0)    = 0;
rel_fit_lb(num_growth' == 0)    = 0;
rel_fit_ub(num_growth' == 0)    = 0;

%%%%%%%%%%%%%%%%%%%%%%%
% Set figure and axes %
%%%%%%%%%%%%%%%%%%%%%%%

fig_rel_fit                     = figure;
set(fig_rel_fit, 'PaperPositionMode', 'auto', 'Units', 'inches', ...
    'Renderer', 'painters', 'Position', [1 1 target_width target_height]);
ax_fig_rel_fit                  = axes('Parent', fig_rel_fit);
ax_fig_rel_fit_2                = axes('Parent', fig_rel_fit);
set(ax_fig_rel_fit, 'Units', 'inches', 'Position', ax_position, 'FontSize', font_size, 'FontName', font_name, ...
    'XLim', [min(x_plot_tick) - 0.5 max(x_plot_tick + 0.5)], 'XTick', x_plot_tick, 'XTickLabel', x_Gal_cell, ...
    'YLim', [0.0 1.3], 'TickLength', [0 0],  ...
    'NextPlot', 'add', 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
set(ax_fig_rel_fit_2, 'Units', 'inches', 'Position', ax_fig_rel_fit.Position, 'FontSize', font_size, 'FontName', font_name, ...
    'XLim', ax_fig_rel_fit.XLim, 'XTick', [], ...
    'YLim', ax_fig_rel_fit.YLim, 'Color', 'none', ...
    'NextPlot', 'add', 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
linkaxes([ax_fig_rel_fit, ax_fig_rel_fit_2])

xlabel(ax_fig_rel_fit,   '% galactose');
if rates_normalized
    ylabel(ax_fig_rel_fit,      'Fitness rel. to WT');
    ylabel(ax_fig_rel_fit_2,    'Fitness rel. to WT');
else
    ylabel(ax_fig_rel_fit,      'Fitness (1/doubling time) [1/min]');
    ylabel(ax_fig_rel_fit_2,    'Fitness (1/doubling time) [1/min]');
end

if data_set == 1 && ~rates_normalized
    ax_fig_rel_fit.YLim     = [0.0 0.015];
    ax_fig_rel_fit_2.YLim   = [0.0 0.015];
elseif data_set == 1
    ax_fig_rel_fit.YLim     = [0.0 2.5];
    ax_fig_rel_fit_2.YLim   = [0.0 2.5];
end
    
%%%%%%%%%%%%%%
% Make plots %
%%%%%%%%%%%%%%

% Error bar lines
rel_fit_lines                   = cell(size(rel_fit, 2), 2 + numel(x_plot_tick));

for i = 1 : size(rel_fit, 2)
    ind_meas                    = ~isnan(num_exp(i, :));
    x                           = x_plot_tick(ind_meas) - 0.3 + (i - 1) * x_shift;
    y                           = rel_fit(ind_meas, i);
    y_low                       = rel_fit_lb(ind_meas, i);
    y_high                      = rel_fit_ub(ind_meas, i);
    pie_data                    = bsxfun(@rdivide, [num_growth(i, :); num_exp(i, :) - num_growth(i, :)], num_exp(i, :));
    pie_data                    = pie_data(:, ind_meas);
    
    if with_pie
        rel_fit_lines{i, 1}     = errorbar(x, y, y_low, y_high, 'Parent', ax_fig_rel_fit, ...
                                    'DisplayName', le_names{i}, 'LineWidth', line_width, 'LineStyle', 'none', ...
                                    'Color', colors_tot(i, :), 'Marker', 'none', 'MarkerSize', mark_size);
    else
        rel_fit_lines{i, 1}     = errorbar(x, y, y_low, y_high, 'Parent', ax_fig_rel_fit, ...
                                    'DisplayName', le_names{i}, 'LineWidth', line_width, 'LineStyle', 'none', ...
                                    'Color', colors_tot(i, :), 'Marker', 'd', 'MarkerSize', mark_size);
    end
    rel_fit_lines{i, 1}.Annotation.LegendInformation.IconDisplayStyle = 'off';
    rel_fit_lines{i, 2}         = line(x, y, 'Parent', ax_fig_rel_fit, 'DisplayName', le_names{i}, ...
                                    'LineWidth', line_width / 3, 'LineStyle', '-', 'Color', colors_tot(i, :));
    
    color_growth                = colors_tot(i, :);
    
    if with_pie
        for j = 1 : numel(y)        
            rel_fit_lines{i, 2 + j} = make_pie_symbol(pie_data(:, j), x(j), y(j), r_pie, ...
                                        {color_growth; color_no_growth}, ax_fig_rel_fit);
        end
    end
end

if with_pie
    marker_pie{1, 1}            = line([1 1], [1e5 1e5], 'Parent', ax_fig_rel_fit, 'DisplayName', 'Growth fraction', ...
                                    'LineWidth', line_width, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', [1 1 1], ...
                                    'MarkerEdgeColor', [1 1 1], 'Color', colors_tot(i,:));
end

% Patches to split mediums
medium_patches                  = cell(numel(x_plot_tick), 1);
for i = 1 : 2 : numel(x_plot_tick)
    medium_patches{i}           = patch(x_plot_tick(i) + [-0.5 0.5 0.5 -0.5], [-1 -1 3 3], [0 0 0], ...
                                    'Parent', ax_fig_rel_fit, 'EdgeColor', 'none', 'FaceAlpha', alpha_patch);
    medium_patches{i}.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% WT reference line
WT_line                         = line([x_plot_tick(1) - 0.5 x_plot_tick(end) + 0.5], [1 1], 'Parent', ax_fig_rel_fit, ...
                                    'Color', rel_fit_lines{1,1}.Color, 'LineWidth', line_width);
WT_line.Annotation.LegendInformation.IconDisplayStyle = 'off';

%%%%%%%%%%%%%%%%%%%%%%%
% Legend and printing %
%%%%%%%%%%%%%%%%%%%%%%%

[le_fitness, le_fitness_icons]          = legend(ax_fig_rel_fit, 'show');
if data_set == 1
    set(le_fitness, 'FontSize', font_size, 'FontName', font_name, 'Location', 'NorthWest', 'Units', 'inches');
else
    set(le_fitness, 'FontSize', font_size, 'FontName', font_name, 'Location', 'SouthEast', 'Units', 'inches');
end
[le_fitness_icons(1 : 4).FontSize]      = deal(font_size);
[le_fitness_icons(1 : 4).FontName]      = deal(font_name);
le_fitness.Position(3)                  = 1.03 * le_fitness.Position(3);

le_pos                          = le_fitness.Position(1 : 2);       % inches, lower left corner rel. to figure
ax_pos                          = ax_fig_rel_fit.Position(1 : 2);   % inches, origin rel. to figure
ax_w                            = ax_fig_rel_fit.Position(3);       % width, inches
ax_h                            = ax_fig_rel_fit.Position(4);       % height, inches
% conversion position icon to units data in axis
icon_pos                        = [le_fitness_icons(numel(Selected_genotypes) + 2).XData(1) + ...
    diff(le_fitness_icons(numel(Selected_genotypes) + 2).XData)/2 ...
    le_fitness_icons(end - 1).YData(1)]; % normalized units inside legend
icons_pos2                      = icon_pos .* le_fitness.Position(3 : 4); % icon position in inches rel. to lower left corner legend
icons_pos3                      = icons_pos2 + le_pos; % icon position in inches rel. to figure
icons_pos4                      = icons_pos3 - ax_pos; % icon position in inches rel. to axis
icons_pos5                      = [ax_fig_rel_fit.XLim(1) ax_fig_rel_fit.YLim(1)] + icons_pos4 ./ [ax_w ax_h] ...
    .* [diff(ax_fig_rel_fit.XLim) diff(ax_fig_rel_fit.YLim)]; % icon position in units data of axis

if with_pie
    marker_pie{2}               = make_pie_symbol([0.5; 0.5], icons_pos5(1) , icons_pos5(2), ...
        0.1, {[0 0 0]; color_no_growth}, ax_fig_rel_fit_2);
    marker_pie{2}(1).Clipping   = 'off';
    delete(marker_pie{2}([2 4]))
end

switch data_set
    case 1
        if filter_suppressors
            if rates_normalized
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Main data without suppressors')
            else
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Main data without suppressors non-normalized')
            end
        else
            if rates_normalized
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Main data')
            else
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Main data non-normalized')
            end
        end
    case 2
        if filter_suppressors
            if rates_normalized
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Supplementary data without suppressors')
            else
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Supplementary data without suppressors non-normalized')
            end
        else
            if rates_normalized
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Supplementary data')
            else
                print(fig_rel_fit, '-dtiffn', strcat('-r', num2str(target_res)), './Supplementary data non-normalized')
            end
        end
end

function pie_handle = make_pie_symbol(data, x, y, r, c, ax)
    % This goal of this function is to make a pie diagram at a specific
    % location in a plot, as a marker.

    warning('off','MATLAB:pie:NonPositiveData')
    pie_handle  = pie(ax, data);
    ax_pos      = ax.Position;
    ax_xlim     = ax.XLim;
    ax_ylim     = ax.YLim;
    ind_data    = find(data, numel(pie_handle)/2);

    for i = 1 : numel(ind_data)
        pie_handle(2 * i - 1).FaceColor = c{ind_data(i)};
        pie_handle(2 * i - 1).XData     = r * pie_handle(2 * i - 1).XData + x;
        pie_handle(2 * i - 1).YData     = r * ax_ylim / ax_xlim * ax_pos(3) / ax_pos(4) * pie_handle(2 * i - 1).YData + y;
        pie_handle(2 * i - 1).EdgeColor = c{1};
        pie_handle(2 * i - 1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        pie_handle(2 * i).String        = '';
    end

end

end