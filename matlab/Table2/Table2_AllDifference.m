clear;
clc;

OAR = readmatrix('TH_OAR_opt_reverse.csv');
OA  = readmatrix('TH_OA_max2_reverse.csv');
OS  = readmatrix('TH_OS_opt.csv');
OAC = readmatrix('TH_cooi_reverse.csv');

all_th_from_excel = [];
all_th_from_oar   = [];

num_files = 5;

for i = 1:num_files
    excel_filename = sprintf('forloop_right_V%d.xlsx', i);
    dataTable = readtable(excel_filename);

    % Read 125 throughput points from Excel (column vector)
    excel_chunk = dataTable.Max_Throughput_Raw;

    % Read the corresponding 125 points from OAR (row vector)
    %oar_chunk_row = OAR(i, :);

    % Read the corresponding 125 points from OA (row vector)
    %oar_chunk_row = OA(i, :);

    % Read the corresponding 125 points from OS (row vector)
    %oar_chunk_row = OS(i, :);

    % Read the corresponding 125 points from OAC (row vector)
    oar_chunk_row = OAC(i, :);

    % Transpose the OAR row vector into a column vector
    oar_chunk_col = oar_chunk_row';

    % Vertically concatenate the current 125 points to the main vectors
    all_th_from_excel = [all_th_from_excel; excel_chunk];
    all_th_from_oar   = [all_th_from_oar;   oar_chunk_col];
end

% Keep only entries where the Excel benchmark is nonzero
keep_indices = (all_th_from_excel ~= 0);

filtered_th_from_excel = all_th_from_excel(keep_indices);
filtered_th_from_oar   = all_th_from_oar(keep_indices);

% Percentage loss relative to the Excel benchmark
Delta = 100 * max(filtered_th_from_excel - filtered_th_from_oar, 0) ...
        ./ filtered_th_from_excel;

% --- Output overall statistics for all combined points ---
fprintf('--- Overall statistics for all combined points ---\n');
fprintf('Mean:    %.4f%%\n', mean(Delta));
fprintf('Median:  %.4f%%\n', median(Delta));
fprintf('Max:     %.4f%%\n', max(Delta));
fprintf('Min:     %.4f%%\n', min(Delta));
fprintf('Q1:      %.4f%%\n', quantile(Delta, 0.25));
fprintf('Q3:      %.4f%%\n', quantile(Delta, 0.75));
fprintf('\n');
