% This script generates the figures (in .tif format) and tables (in .txt
% format) regarding the growth rate data starting from processed plate
% reader data.

%% Figures
Growth_data_analysis(1, false, true, false)     % Figure 3A analogue
Growth_data_analysis(1, false, false, false)    % Supplementary Figure S5

Growth_data_analysis(2, true, true, true)       % Supplementary Figure S4
Growth_data_analysis(2, true, true, false)      % Figure not relevant, just the .mat file needed later

%% Tables
% Supplementary table S5
load('Main_data.mat')
data_for_table  = cell2mat(cellfun(@(x) x(1 : 2)', Summary_cell(2 : end, 2: end), 'UniformOutput', false))';
Col_headers     = strrep(strrep(Summary_cell(2 : end, 1), ' ', '_'), '-', '_');
Col_headers     = [{'Medium'}; reshape([strrep(Col_headers, '42', '42_TG') strrep(Col_headers, '42', '42_R')]', [], 1)];
table_S5        = cell2table([Summary_cell(1, 2 : end)' num2cell(data_for_table)], 'VariableNames', Col_headers');
writetable(table_S5, 'Supplementary_table_S5.txt', 'Delimiter', '\t');

% Supplementary table S6
load('Supplementary_data_without_suppressors.mat')
data_for_table  = cell2mat(cellfun(@(x) x(1), Summary_cell(2 : end, 2: end), 'UniformOutput', false))';
data_for_table2 = cell2mat(cellfun(@(x) x(2), Summary_cell(2 : end, 2: end), 'UniformOutput', false))';
load('Supplementary_data.mat')
data_for_table3 = cell2mat(cellfun(@(x) x(1), Summary_cell(2 : end, 2: end), 'UniformOutput', false))' - data_for_table;
data_for_table4 = reshape([data_for_table3; data_for_table; data_for_table2], size(data_for_table, 1), []);

Col_headers     = strrep(strrep(Summary_cell(2 : end, 1), ' ', '_'), '-', '_');
Col_headers     = [{'Medium'}; reshape([strrep(Col_headers, '42', '42_D') ...
                    strrep(Col_headers, '42', '42_TG') strrep(Col_headers, '42', '42_R')]', [], 1)];
table_S6        = cell2table([Summary_cell(1, 2 : end)' num2cell(data_for_table4)], 'VariableNames', Col_headers');
writetable(table_S6, 'Supplementary_table_S6.txt', 'Delimiter', '\t');

% Supplementary table S7
load('Supplementary_data.mat')
sfGFP_col               = find(cellfun(@(x) ~isempty(regexp(x, 'GAL1-sfGFP-CDC42', 'match')), Genotypes)) + 1;
min_sfGFP_col           = find(cellfun(@(x) ~isempty(regexp(x, 'GAL1-CDC42', 'match')), Genotypes)) + 1;
sfGFP_comparison        = NaN(2, size(Samples, 2) - 1);

for i = 2 : size(Samples, 2)
    if ~(any(isnan(Samples{min_sfGFP_col, i})) || any(isnan(Samples{sfGFP_col, i})))
        frac_larger         = mean(Samples{min_sfGFP_col, i} > Samples{sfGFP_col, i});
        sfGFP_comparison(1, i - 1) = max(frac_larger / (1 - frac_larger), (1 - frac_larger) / frac_larger);
    end
end

load('Supplementary_data_without_suppressors.mat')
for i = 2 : size(Samples, 2)
    if ~(any(isnan(Samples{min_sfGFP_col, i})) || any(isnan(Samples{sfGFP_col, i})))
        frac_larger         = mean(Samples{min_sfGFP_col, i} > Samples{sfGFP_col, i});
        sfGFP_comparison(2, i - 1) = max(frac_larger / (1 - frac_larger), (1 - frac_larger) / frac_larger);
    end
end

ind_notNaN              = ~any(isnan(sfGFP_comparison));
table_S7                = cell2table([Media(ind_notNaN)'; num2cell(round(sfGFP_comparison(:, ind_notNaN)))], ...
    'RowNames', {'Medium'; 'Bayes_factor_without_filtering'; 'Bayes_factor_with_filtering'});
writetable(table_S7, 'Supplementary_table_S7.txt', 'Delimiter', '\t');