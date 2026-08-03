% SimplePunctaAnalysis_v4.m
% Interactive tool for synaptic puncta engulfment analysis
%
% WORKFLOW:
%   1. Run in 'tune' mode on one representative dataset
%   2. Adjust sliders until segmentation looks right
%   3. Save parameters
%   4. Run in 'batch' mode to process all datasets
%
% PHILOSOPHY: Simple per-image processing with visual verification
%
% v2 changes:
%   - Added binary mask output (mask_green, mask_red, mask_ctb, 
%     mask_red_engulfed, mask_ctb_engulfed)
%
% v3 changes:
%   - Added Maximum_intensity_projections folder per dataset with:
%     * MIPs of all mask stacks
%     * Raw/normalized RGB MIPs (original + adjusted color schemes)
%     * Masked RGB MIPs
%     * Engulfed-only MIPs (glial + engulfed puncta)
%     * All with/without 10µm scalebar, 900 DPI TIFF
%   - Added montage generation per sample/slide in root directories
%
% v4 changes:
%   - Fixed montage generation to group by Sample1A-7A and Sample1B-7B
%   - Added 4 montage types per group:
%     * RGB Original MIP (R=red, G=green, B=CTB)
%     * RGB Adjusted MIP (grey=green, magenta=red, green=CTB)
%     * Green channel MIP (grayscale)
%     * Green channel binary mask MIP (0/255)
%   - Labels above each image with full filename (wrapped) + section count
%   - Montages saved to Sample_A/Montages/ and Sample_B/Montages/
%   - Fixed EPS output to be true vector format (removed transparency which
%     caused rasterization; uses painters renderer)
%   - Added PDF output alongside EPS for better Illustrator compatibility

clear; clc;
fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
fprintf('  SimplePunctaAnalysis v4 - Interactive Tuning + Batch\n');
fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');

%% ===================== USER CONFIG =====================

% --- Bio-Formats path (required for reading voxel sizes from CZI files) ---
BFMATLAB_DIR = 'C:/Program Files/MATLAB/R2025a/bfmatlab/';

% --- MODE ---
% 'tune'  : Interactive parameter tuning on single dataset
% 'batch' : Process all datasets with saved/specified parameters
MODE = 'batch';

% --- PATHS ---
% For TUNE mode: path to one representative dataset
TUNE_DATASET = 'Y:\Colenso\2026\Global_PerSample\Sample_A\Sample_2A\Top_right_example_2';

% Alternative dataset for testing:
% TUNE_DATASET = 'Y:\Colenso\2026\Global_PerSample\Sample_B\Sample_6B\bottom_left_cell_4';

% For BATCH mode: root folder containing Sample_A and Sample_B
BATCH_ROOTS = { 'Y:\Colenso\2026\Global_PerSample' };

% Output folder for batch results
BATCH_OUTPUT = 'Y:\Colenso\2026\Global_PerSample\BatchResults';

% --- PARAMETERS (adjust in tune mode, then copy here for batch) ---
params = struct();

% Stretch percentiles (per-image normalization)
params.stretch_low = 1;      % percentile mapped to 0
params.stretch_high = 99.5;  % percentile mapped to 255

% Gaussian blur sigma (pixels)
params.blur_sigma = 1.5;

% --- THRESHOLD METHOD ---
% 'fixed'      : Use fixed threshold values (requires tuning per-image - NOT recommended)
% 'percentile' : Threshold at Nth percentile of image (adapts to brightness - RECOMMENDED)
% 'otsu'       : Automatic Otsu threshold per image
params.threshold_method = 'percentile';

% --- GREEN CHANNEL (microglia) ---
params.green_threshold = 0.25;          % For 'fixed' method only
params.green_percentile = 97.4;         % For 'percentile': top 2.6% of pixels
params.green_min_size_um3 = 50;         % minimum object volume in ÂµmÂ³

% --- RED CHANNEL (vGlut2 puncta) ---
params.red_threshold = 0.50;            % For 'fixed' method only
params.red_percentile = 99.0;           % For 'percentile': top 1% (puncta are small/sparse)
params.red_min_size_um3 = 0.1;          % minimum puncta volume
params.red_max_size_um3 = 300.0;         % maximum puncta volume

% --- CTB CHANNEL ---
params.ctb_threshold = 0.25;            % For 'fixed' method only
params.ctb_percentile = 99.0;           % For 'percentile': top 3.1%
params.ctb_min_size_um3 = 0.1;
params.ctb_max_size_um3 = 300.0;

% --- PARAMETER FILE (for auto-save/load between tune and batch modes) ---
PARAMS_FILE = fullfile(BATCH_OUTPUT, 'tuned_parameters.mat');

% --- QC OUTPUT OPTIONS ---
% Save QC images to verify segmentation quality
params.qc.enabled = true;               % Master switch for QC outputs
params.qc.save_all = true;              % true = save for all datasets, false = save subset
params.qc.subset_n = 10;                % If save_all=false, save first N from each group
params.qc.save_overlays = true;         % MIP with mask overlay for each channel
params.qc.save_histograms = true;       % Intensity histogram with threshold line
params.qc.save_composite = true;        % Combined RGB overlay of all channels

% --- MASKED IMAGE OUTPUT ---
% Save raw intensity images with mask applied (signal that passed filtering)
params.save_masked_images = true;       % Save masked image stacks for each channel

% --- BINARY MASK OUTPUT (v2) ---
% Save binary masks as 8-bit images (0=background, 255=mask)
params.save_binary_masks = true;        % Save binary mask stacks for each channel

% --- MIP OUTPUT OPTIONS (v3) ---
params.save_mips = true;                % Save maximum intensity projection images
params.mip_dpi = 900;                   % Resolution for MIP TIFF outputs
params.scalebar_um = 10;                % Scalebar length in microns

% --- MONTAGE OUTPUT OPTIONS (v3) ---
params.save_montages = true;            % Save montage images per sample/slide

%% ===================== INITIALIZATION =====================

% Load Bio-Formats (required for reading voxel sizes from CZI files)
if isfolder(BFMATLAB_DIR)
    addpath(genpath(BFMATLAB_DIR));
    try
        bfCheckJavaPath(true);
        fprintf('Bio-Formats loaded successfully.\n\n');
    catch ME
        warning('Bio-Formats could not be initialized: %s', ME.message);
        fprintf('Voxel sizes will be read from TIFF metadata if available.\n\n');
    end
else
    warning('Bio-Formats directory not found: %s', BFMATLAB_DIR);
    fprintf('Voxel sizes will be read from TIFF metadata if available.\n\n');
end

%% ===================== RUN =====================

switch MODE
    case 'tune'
        run_interactive_tuning(TUNE_DATASET, params, PARAMS_FILE);
        
    case 'batch'
        % Try to load tuned parameters if available
        if isfile(PARAMS_FILE)
            fprintf('Loading tuned parameters from: %s\n', PARAMS_FILE);
            loaded = load(PARAMS_FILE);
            
            % Update only the tunable parameters
            params.green_percentile = loaded.params.green_percentile;
            params.red_percentile = loaded.params.red_percentile;
            params.red_max_size_um3 = loaded.params.red_max_size_um3;
            params.ctb_percentile = loaded.params.ctb_percentile;
            params.ctb_max_size_um3 = loaded.params.ctb_max_size_um3;
            
            fprintf('  Green %%ile: %.1f\n', params.green_percentile);
            fprintf('  Red %%ile: %.1f, max size: %.1f ÂµmÂ³\n', params.red_percentile, params.red_max_size_um3);
            fprintf('  CTB %%ile: %.1f, max size: %.1f ÂµmÂ³\n\n', params.ctb_percentile, params.ctb_max_size_um3);
        else
            fprintf('No tuned parameters file found. Using defaults.\n');
            fprintf('Run in tune mode first to optimize parameters.\n\n');
        end
        
        run_batch_processing(BATCH_ROOTS, BATCH_OUTPUT, params);
        
    otherwise
        error('Unknown MODE: %s', MODE);
end

%% =================== INTERACTIVE TUNING ===================

function run_interactive_tuning(dataset_path, params, params_file)
    
    fprintf('Loading dataset: %s\n', dataset_path);
    
    % Load stacks
    green_raw = load_stack(fullfile(dataset_path, 'green'));
    red_raw = load_stack(fullfile(dataset_path, 'red'));
    ctb_raw = load_stack(fullfile(dataset_path, 'blue'));
    
    [H, W, nZ] = size(green_raw);
    fprintf('  Size: %d x %d x %d\n', H, W, nZ);
    
    % Get voxel size
    voxel_um = get_voxel_size(dataset_path);
    voxel_vol = prod(voxel_um);
    fprintf('  Voxel: %.3f x %.3f x %.3f Âµm (vol=%.4f ÂµmÂ³)\n', voxel_um, voxel_vol);
    
    % Normalize (per-image stretch)
    fprintf('  Normalizing...\n');
    green_norm = normalize_stack(green_raw, params.stretch_low, params.stretch_high);
    red_norm = normalize_stack(red_raw, params.stretch_low, params.stretch_high);
    ctb_norm = normalize_stack(ctb_raw, params.stretch_low, params.stretch_high);
    
    % Apply blur
    green_blur = apply_blur_3d(green_norm, params.blur_sigma);
    red_blur = apply_blur_3d(red_norm, params.blur_sigma);
    ctb_blur = apply_blur_3d(ctb_norm, params.blur_sigma);
    
    % Create MIPs for display
    green_mip = max(green_norm, [], 3);
    red_mip = max(red_norm, [], 3);
    ctb_mip = max(ctb_norm, [], 3);
    
    % Create blurred MIPs for fast 2D preview
    fprintf('  Creating 2D preview images...\n');
    green_mip_blur = max(green_blur, [], 3);
    red_mip_blur = max(red_blur, [], 3);
    ctb_mip_blur = max(ctb_blur, [], 3);
    
    % Create figure
    fig = figure('Name', 'Puncta Analysis - Parameter Tuning', ...
        'NumberTitle', 'off', 'Position', [50 50 1600 900]);
    
    % Store data
    data.green_blur = green_blur;
    data.red_blur = red_blur;
    data.ctb_blur = ctb_blur;
    data.green_mip = green_mip;
    data.red_mip = red_mip;
    data.ctb_mip = ctb_mip;
    data.green_mip_blur = green_mip_blur;
    data.red_mip_blur = red_mip_blur;
    data.ctb_mip_blur = ctb_mip_blur;
    data.voxel_vol = voxel_vol;
    data.voxel_um = voxel_um;
    data.params = params;
    data.params_file = params_file;
    data.H = H; data.W = W; data.nZ = nZ;
    
    % Create axes for display
    data.ax_green = subplot(2, 4, 1);
    data.ax_red = subplot(2, 4, 2);
    data.ax_ctb = subplot(2, 4, 3);
    data.ax_overlay = subplot(2, 4, 4);
    data.ax_green_mask = subplot(2, 4, 5);
    data.ax_red_mask = subplot(2, 4, 6);
    data.ax_ctb_mask = subplot(2, 4, 7);
    data.ax_stats = subplot(2, 4, 8);
    
    % Show raw MIPs
    axes(data.ax_green);
    imshow(green_mip, []); title('Green (raw MIP)');
    
    axes(data.ax_red);
    imshow(red_mip, []); title('Red (raw MIP)');
    
    axes(data.ax_ctb);
    imshow(ctb_mip, []); title('CTB (raw MIP)');
    
    % Create slider panel
    panel = uipanel('Title', 'Parameters (Percentile-based thresholding)', 'Position', [0.01 0.01 0.98 0.12]);
    
    % GREEN sliders - percentile
    uicontrol(panel, 'Style', 'text', 'String', 'GREEN %ile:', ...
        'Units', 'normalized', 'Position', [0.01 0.6 0.08 0.3]);
    data.slider_green_pct = uicontrol(panel, 'Style', 'slider', ...
        'Min', 80, 'Max', 99.9, 'Value', params.green_percentile, ...
        'Units', 'normalized', 'Position', [0.09 0.65 0.12 0.25], ...
        'Callback', @(~,~) update_display(fig));
    data.text_green_pct = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.1f', params.green_percentile), ...
        'Units', 'normalized', 'Position', [0.21 0.6 0.04 0.3]);
    
    % RED sliders - percentile and size
    uicontrol(panel, 'Style', 'text', 'String', 'RED %ile:', ...
        'Units', 'normalized', 'Position', [0.26 0.6 0.07 0.3]);
    data.slider_red_pct = uicontrol(panel, 'Style', 'slider', ...
        'Min', 80, 'Max', 99.9, 'Value', params.red_percentile, ...
        'Units', 'normalized', 'Position', [0.33 0.65 0.12 0.25], ...
        'Callback', @(~,~) update_display(fig));
    data.text_red_pct = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.1f', params.red_percentile), ...
        'Units', 'normalized', 'Position', [0.45 0.6 0.04 0.3]);
    
    uicontrol(panel, 'Style', 'text', 'String', 'RED max ÂµmÂ³:', ...
        'Units', 'normalized', 'Position', [0.26 0.1 0.07 0.3]);
    data.slider_red_max = uicontrol(panel, 'Style', 'slider', ...
        'Min', 1, 'Max', 150, 'Value', params.red_max_size_um3, ...
        'Units', 'normalized', 'Position', [0.33 0.15 0.12 0.25], ...
        'Callback', @(~,~) update_display(fig));
    data.text_red_max = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.1f', params.red_max_size_um3), ...
        'Units', 'normalized', 'Position', [0.45 0.1 0.04 0.3]);
    
    % CTB sliders - percentile and size
    uicontrol(panel, 'Style', 'text', 'String', 'CTB %ile:', ...
        'Units', 'normalized', 'Position', [0.50 0.6 0.07 0.3]);
    data.slider_ctb_pct = uicontrol(panel, 'Style', 'slider', ...
        'Min', 80, 'Max', 99.9, 'Value', params.ctb_percentile, ...
        'Units', 'normalized', 'Position', [0.57 0.65 0.12 0.25], ...
        'Callback', @(~,~) update_display(fig));
    data.text_ctb_pct = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.1f', params.ctb_percentile), ...
        'Units', 'normalized', 'Position', [0.69 0.6 0.04 0.3]);
    
    uicontrol(panel, 'Style', 'text', 'String', 'CTB max ÂµmÂ³:', ...
        'Units', 'normalized', 'Position', [0.50 0.1 0.07 0.3]);
    data.slider_ctb_max = uicontrol(panel, 'Style', 'slider', ...
        'Min', 1, 'Max', 200, 'Value', params.ctb_max_size_um3, ...
        'Units', 'normalized', 'Position', [0.57 0.15 0.12 0.25], ...
        'Callback', @(~,~) update_display(fig));
    data.text_ctb_max = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.1f', params.ctb_max_size_um3), ...
        'Units', 'normalized', 'Position', [0.69 0.1 0.04 0.3]);
    
    % Buttons
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Calculate 3D Stats', ...
        'Units', 'normalized', 'Position', [0.75 0.55 0.1 0.35], ...
        'BackgroundColor', [0.3 0.6 0.9], 'FontWeight', 'bold', ...
        'Callback', @(~,~) calculate_3d_stats(fig));
    
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Print Parameters', ...
        'Units', 'normalized', 'Position', [0.75 0.1 0.1 0.35], ...
        'Callback', @(~,~) print_params(fig));
    
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Save & Close', ...
        'Units', 'normalized', 'Position', [0.86 0.3 0.1 0.4], ...
        'BackgroundColor', [0.4 0.8 0.4], ...
        'Callback', @(~,~) save_and_close(fig));
    
    guidata(fig, data);
    
    % Initial display
    update_display(fig);
    
    fprintf('\n  âœ“ Interactive tuning ready.\n');
    fprintf('    1. Adjust sliders until masks look correct (fast 2D preview)\n');
    fprintf('    2. Click "Calculate 3D Stats" to compute actual engulfment\n');
    fprintf('    3. Click "Print Parameters" to get values for batch processing\n');
end

function update_display(fig)
    data = guidata(fig);
    
    % Get current parameter values (percentiles)
    green_pct = data.slider_green_pct.Value;
    red_pct = data.slider_red_pct.Value;
    red_max = data.slider_red_max.Value;
    ctb_pct = data.slider_ctb_pct.Value;
    ctb_max = data.slider_ctb_max.Value;
    
    % Update text displays
    data.text_green_pct.String = sprintf('%.1f', green_pct);
    data.text_red_pct.String = sprintf('%.1f', red_pct);
    data.text_red_max.String = sprintf('%.1f', red_max);
    data.text_ctb_pct.String = sprintf('%.1f', ctb_pct);
    data.text_ctb_max.String = sprintf('%.1f', ctb_max);
    
    % === FAST 2D PROCESSING ON MIPs ===
    % Compute percentile-based thresholds on MIP images
    green_vals = double(data.green_mip_blur(data.green_mip_blur > 0));
    red_vals = double(data.red_mip_blur(data.red_mip_blur > 0));
    ctb_vals = double(data.ctb_mip_blur(data.ctb_mip_blur > 0));
    
    green_thresh = prctile(green_vals, green_pct);
    red_thresh = prctile(red_vals, red_pct);
    ctb_thresh = prctile(ctb_vals, ctb_pct);
    
    % Threshold on MIP images (instant)
    green_mask_mip = data.green_mip_blur > green_thresh;
    red_mask_mip = data.red_mip_blur > red_thresh;
    ctb_mask_mip = data.ctb_mip_blur > ctb_thresh;
    
    % Quick 2D size filter (for visual feedback only)
    green_mask_mip = bwareaopen(green_mask_mip, 500);
    red_mask_mip = bwareaopen(red_mask_mip, 5);
    ctb_mask_mip = bwareaopen(ctb_mask_mip, 10);
    
    % Display mask overlays
    axes(data.ax_green_mask);
    show_mask_overlay(data.green_mip, green_mask_mip, [0 1 0]);
    title(sprintf('Green mask (2D)\n%ile=%.1f â†’ thresh=%.0f', green_pct, green_thresh));
    
    axes(data.ax_red_mask);
    show_mask_overlay(data.red_mip, red_mask_mip, [1 0 0]);
    title(sprintf('Red mask (2D)\n%ile=%.1f â†’ thresh=%.0f', red_pct, red_thresh));
    
    axes(data.ax_ctb_mask);
    show_mask_overlay(data.ctb_mip, ctb_mask_mip, [0 0.5 1]);
    title(sprintf('CTB mask (2D)\n%ile=%.1f â†’ thresh=%.0f', ctb_pct, ctb_thresh));
    
    % Combined overlay
    axes(data.ax_overlay);
    rgb = create_rgb_overlay(data.green_mip, data.red_mip, data.ctb_mip, ...
        green_mask_mip, red_mask_mip, ctb_mask_mip);
    imshow(rgb);
    title('Combined overlay (2D)');
    
    % Stats display - show message to run 3D analysis
    axes(data.ax_stats);
    cla;
    axis off;
    
    stats_text = {
        'â•â•â• 2D PREVIEW MODE â•â•â•', ...
        '', ...
        sprintf('Green: top %.1f%% of pixels', 100-green_pct), ...
        sprintf('Red: top %.1f%% of pixels', 100-red_pct), ...
        sprintf('CTB: top %.1f%% of pixels', 100-ctb_pct), ...
        '', ...
        'Click "Calculate 3D Stats"', ...
        'for actual engulfment numbers.', ...
        '', ...
        'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•'
    };
    
    text(0.1, 0.9, stats_text, 'FontName', 'FixedWidth', 'FontSize', 10, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    guidata(fig, data);
end

function calculate_3d_stats(fig)
    % Run full 3D analysis with current parameters
    data = guidata(fig);
    
    fprintf('Calculating 3D statistics...\n');
    
    % Get current percentile values from sliders
    green_pct = data.slider_green_pct.Value;
    red_pct = data.slider_red_pct.Value;
    red_max = data.slider_red_max.Value;
    ctb_pct = data.slider_ctb_pct.Value;
    ctb_max = data.slider_ctb_max.Value;
    
    % Convert size limits to voxels
    voxel_vol = data.voxel_vol;
    green_min_vox = data.params.green_min_size_um3 / voxel_vol;
    red_min_vox = data.params.red_min_size_um3 / voxel_vol;
    red_max_vox = red_max / voxel_vol;
    ctb_min_vox = data.params.ctb_min_size_um3 / voxel_vol;
    ctb_max_vox = ctb_max / voxel_vol;
    
    % Compute percentile thresholds from 3D data
    fprintf('  Computing thresholds...\n');
    green_vals = double(data.green_blur(data.green_blur > 0));
    red_vals = double(data.red_blur(data.red_blur > 0));
    ctb_vals = double(data.ctb_blur(data.ctb_blur > 0));
    
    green_thresh = prctile(green_vals, green_pct);
    red_thresh = prctile(red_vals, red_pct);
    ctb_thresh = prctile(ctb_vals, ctb_pct);
    
    fprintf('    Green: %ile=%.1f â†’ thresh=%.1f\n', green_pct, green_thresh);
    fprintf('    Red:   %ile=%.1f â†’ thresh=%.1f\n', red_pct, red_thresh);
    fprintf('    CTB:   %ile=%.1f â†’ thresh=%.1f\n', ctb_pct, ctb_thresh);
    
    % Full 3D segmentation
    fprintf('  Segmenting green...\n');
    green_mask = data.green_blur > green_thresh;
    green_mask = filter_by_size(green_mask, green_min_vox, inf);
    
    fprintf('  Segmenting red...\n');
    red_mask = data.red_blur > red_thresh;
    red_mask = filter_by_size(red_mask, red_min_vox, red_max_vox);
    
    fprintf('  Segmenting CTB...\n');
    ctb_mask = data.ctb_blur > ctb_thresh;
    ctb_mask = filter_by_size(ctb_mask, ctb_min_vox, ctb_max_vox);
    
    % Compute statistics
    fprintf('  Computing statistics...\n');
    n_green = nnz(green_mask);
    n_red_in_green = nnz(red_mask & green_mask);
    n_ctb_in_green = nnz(ctb_mask & green_mask);
    
    if n_green > 0
        pct_red = 100 * n_red_in_green / n_green;
        pct_ctb = 100 * n_ctb_in_green / n_green;
    else
        pct_red = 0;
        pct_ctb = 0;
    end
    
    % Count objects
    CC_red = bwconncomp(red_mask & green_mask, 26);
    CC_ctb = bwconncomp(ctb_mask & green_mask, 26);
    CC_red_total = bwconncomp(red_mask, 26);
    n_red_obj = CC_red.NumObjects;
    n_ctb_obj = CC_ctb.NumObjects;
    n_red_total = CC_red_total.NumObjects;
    
    fprintf('  Done!\n\n');
    
    % Update MIP displays with 3D mask projections
    green_mask_mip = max(green_mask, [], 3);
    red_mask_mip = max(red_mask, [], 3);
    ctb_mask_mip = max(ctb_mask, [], 3);
    
    axes(data.ax_green_mask);
    show_mask_overlay(data.green_mip, green_mask_mip, [0 1 0]);
    title(sprintf('Green mask (3D)\n%d voxels', n_green));
    
    axes(data.ax_red_mask);
    show_mask_overlay(data.red_mip, red_mask_mip, [1 0 0]);
    title(sprintf('Red mask (3D)\n%d total objects', n_red_total));
    
    axes(data.ax_ctb_mask);
    show_mask_overlay(data.ctb_mip, ctb_mask_mip, [0 0.5 1]);
    title(sprintf('CTB mask (3D)\n%d objects', bwconncomp(ctb_mask, 26).NumObjects));
    
    % Combined overlay with 3D masks
    axes(data.ax_overlay);
    rgb = create_rgb_overlay(data.green_mip, data.red_mip, data.ctb_mip, ...
        green_mask_mip, red_mask_mip, ctb_mask_mip);
    imshow(rgb);
    title('Combined overlay (3D masks)');
    
    % Stats display
    axes(data.ax_stats);
    cla;
    axis off;
    
    stats_text = {
        'â•â•â• 3D ENGULFMENT RESULTS â•â•â•', ...
        '', ...
        sprintf('Green volume: %.0f voxels', n_green), ...
        sprintf('             (%.2f ÂµmÂ³)', n_green * voxel_vol), ...
        '', ...
        'â”€â”€â”€ RED (vGlut2) â”€â”€â”€', ...
        sprintf('  Engulfment: %.4f %%', pct_red), ...
        sprintf('  Objects in green: %d', n_red_obj), ...
        sprintf('  Total objects: %d', n_red_total), ...
        '', ...
        'â”€â”€â”€ CTB â”€â”€â”€', ...
        sprintf('  Engulfment: %.4f %%', pct_ctb), ...
        sprintf('  Objects in green: %d', n_ctb_obj), ...
        '', ...
        'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•'
    };
    
    text(0.1, 0.95, stats_text, 'FontName', 'FixedWidth', 'FontSize', 10, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % Color-code warnings
    if pct_red > 10
        text(0.1, 0.05, 'âš  Red engulfment >10% - check threshold!', ...
            'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 10);
    elseif pct_red < 0.1
        text(0.1, 0.05, 'âš  Red engulfment very low - threshold too high?', ...
            'Color', [0.8 0.5 0], 'FontWeight', 'bold', 'FontSize', 10);
    else
        text(0.1, 0.05, 'âœ“ Red engulfment in reasonable range', ...
            'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 10);
    end
    
    % Print to console too
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('  3D ENGULFMENT RESULTS\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('  Green volume: %.0f voxels (%.2f ÂµmÂ³)\n', n_green, n_green * voxel_vol);
    fprintf('  Red engulfment:  %.4f %% (%d objects in green, %d total)\n', pct_red, n_red_obj, n_red_total);
    fprintf('  CTB engulfment:  %.4f %% (%d objects in green)\n', pct_ctb, n_ctb_obj);
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
end

function print_params(fig)
    data = guidata(fig);
    
    fprintf('\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('  COPY THESE PARAMETERS TO BATCH PROCESSING:\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('params.threshold_method = ''percentile'';\n');
    fprintf('params.green_percentile = %.1f;\n', data.slider_green_pct.Value);
    fprintf('params.red_percentile = %.1f;\n', data.slider_red_pct.Value);
    fprintf('params.red_max_size_um3 = %.1f;\n', data.slider_red_max.Value);
    fprintf('params.ctb_percentile = %.1f;\n', data.slider_ctb_pct.Value);
    fprintf('params.ctb_max_size_um3 = %.1f;\n', data.slider_ctb_max.Value);
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
end

function save_and_close(fig)
    data = guidata(fig);
    
    % Create params struct with current values
    params = data.params;
    params.green_percentile = data.slider_green_pct.Value;
    params.red_percentile = data.slider_red_pct.Value;
    params.red_max_size_um3 = data.slider_red_max.Value;
    params.ctb_percentile = data.slider_ctb_pct.Value;
    params.ctb_max_size_um3 = data.slider_ctb_max.Value;
    
    % Save to file
    params_file = data.params_file;
    
    % Ensure directory exists
    [params_dir, ~, ~] = fileparts(params_file);
    if ~isfolder(params_dir)
        mkdir(params_dir);
    end
    
    save(params_file, 'params');
    
    fprintf('\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('  PARAMETERS SAVED!\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('  File: %s\n', params_file);
    fprintf('\n');
    fprintf('  These parameters will be automatically loaded in batch mode.\n');
    fprintf('  Just set MODE = ''batch'' and run the script.\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    
    print_params(fig);
    close(fig);
end

%% =================== BATCH PROCESSING ===================

function run_batch_processing(roots, output_dir, params)
    
    fprintf('Scanning for datasets...\n');
    
    % Find all datasets
    datasets = find_all_datasets(roots);
    n_total = numel(datasets);
    
    fprintf('Found %d datasets\n\n', n_total);
    
    if n_total == 0
        error('No datasets found!');
    end
    
    % === TEST CONNECTIVITY ===
    % Try to read one file from first dataset to verify access
    fprintf('Testing file access...\n');
    test_folder = fullfile(datasets(1).path, 'green');
    test_files = dir(fullfile(test_folder, '*.tif'));
    if isempty(test_files)
        test_files = dir(fullfile(test_folder, '*.TIF'));
    end
    
    if ~isempty(test_files)
        test_file = fullfile(test_folder, test_files(1).name);
        fprintf('  Testing: %s\n', test_file);
        try
            test_img = imread(test_file);
            fprintf('  âœ“ File access OK (%dx%d image)\n\n', size(test_img, 1), size(test_img, 2));
            clear test_img;
        catch ME
            fprintf('\n');
            fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
            fprintf('  ERROR: Cannot read files from network drive!\n');
            fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
            fprintf('  Test file: %s\n', test_file);
            fprintf('  Error: %s\n', ME.message);
            fprintf('\n');
            fprintf('  Possible causes:\n');
            fprintf('    1. Network drive disconnected - check Y:\\ in Explorer\n');
            fprintf('    2. Files locked by another program (ImageJ, MATLAB, etc.)\n');
            fprintf('    3. Antivirus blocking access\n');
            fprintf('    4. Permission issue\n');
            fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
            error('Cannot access files. Fix network/permission issue and retry.');
        end
    end
    
    % Create output directory
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    
    % Create QC directory if enabled
    qc_dir = fullfile(output_dir, 'QC');
    if params.qc.enabled && ~isfolder(qc_dir)
        mkdir(qc_dir);
    end
    
    % Track counts for QC subset
    qc_count_A = 0;
    qc_count_B = 0;
    
    % Process each dataset
    % Collect all stats in cell array first
    all_stats = {};
    
    % Collect montage data for each dataset (v4 - updated fields)
    montage_data = struct('group', {}, 'replicate', {}, 'name', {}, ...
        'nZ', {}, 'mip_rgb_original', {}, 'mip_rgb_adjusted', {}, ...
        'mip_green_grey', {}, 'mip_green_mask', {});
    
    for i = 1:n_total
        ds = datasets(i);
        fprintf('[%d/%d] %s... ', i, n_total, ds.name);
        
        % Determine if we should save QC for this dataset
        save_qc = false;
        if params.qc.enabled
            if params.qc.save_all
                save_qc = true;
            else
                % Save for first N of each group
                if strcmp(ds.group, 'A') && qc_count_A < params.qc.subset_n
                    save_qc = true;
                    qc_count_A = qc_count_A + 1;
                elseif strcmp(ds.group, 'B') && qc_count_B < params.qc.subset_n
                    save_qc = true;
                    qc_count_B = qc_count_B + 1;
                end
            end
        end
        
        try
            % Set up QC info
            qc_info = struct();
            qc_info.enabled = save_qc;
            qc_info.output_dir = fullfile(qc_dir, sprintf('%s_%s', ds.group, ds.name));
            qc_info.dataset_name = ds.name;
            qc_info.params = params.qc;
            
            [stats, mip_info] = process_single_dataset(ds.path, params, qc_info);
            
            % Add metadata
            stats.dataset = ds.name;
            stats.group = ds.group;
            stats.replicate = ds.replicate;
            stats.path = ds.path;
            
            all_stats{end+1} = stats; %#ok<AGROW>
            
            % Store montage data if MIPs were saved (v4 - updated fields)
            if params.save_montages && ~isempty(mip_info)
                md_idx = numel(montage_data) + 1;
                montage_data(md_idx).group = ds.group;
                montage_data(md_idx).replicate = ds.replicate;
                montage_data(md_idx).name = ds.name;
                montage_data(md_idx).nZ = mip_info.nZ;
                montage_data(md_idx).mip_rgb_original = mip_info.mip_rgb_original;
                montage_data(md_idx).mip_rgb_adjusted = mip_info.mip_rgb_adjusted;
                montage_data(md_idx).mip_green_grey = mip_info.mip_green_grey;
                montage_data(md_idx).mip_green_mask = mip_info.mip_green_mask;
            end
            
            fprintf('Red: %.3f%% (thresh=%.1f), CTB: %.3f%%\n', ...
                stats.red_engulfment_pct, stats.red_threshold_used, stats.ctb_engulfment_pct);
            
        catch ME
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    
    % Check if we have results
    if isempty(all_stats)
        fprintf('\nERROR: No datasets were successfully processed!\n');
        return;
    end
    
    % Convert to table
    fprintf('\nConverting results to table...\n');
    
    % Get field names from first result
    fnames = fieldnames(all_stats{1});
    n_results = numel(all_stats);
    
    % Build table manually for robustness
    results = table();
    for f = 1:numel(fnames)
        fn = fnames{f};
        vals = cell(n_results, 1);
        for r = 1:n_results
            vals{r} = all_stats{r}.(fn);
        end
        
        % Check if numeric or string
        if isnumeric(vals{1})
            results.(fn) = cell2mat(vals);
        else
            results.(fn) = vals;
        end
    end
    
    % Save results (with fallback for locked files)
    % Rename channel columns to color+marker convention for on-disk clarity
    % (red=VGluT2, green=microglia, blue=CTB). Internal processing above uses the
    % original short names; this only affects the exported CSV headers.
    results_out = rename_channel_columns(results);
    results_file = fullfile(output_dir, 'engulfment_results.csv');
    try
        writetable(results_out, results_file);
        fprintf('Results saved to: %s\n', results_file);
    catch
        % File might be open in Excel - try timestamped name
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        results_file = fullfile(output_dir, sprintf('engulfment_results_%s.csv', timestamp));
        writetable(results_out, results_file);
        fprintf('WARNING: Original CSV was locked. Results saved to: %s\n', results_file);
    end
    
    % Generate summary
    generate_summary(results, output_dir, params);
    
    % === GENERATE MONTAGES (v3) ===
    if params.save_montages && ~isempty(montage_data)
        fprintf('\nGenerating montage images...\n');
        generate_montages(montage_data, roots{1}, params);
    end
    
    fprintf('\nâœ“ Batch processing complete!\n');
    fprintf('  Successful: %d datasets\n', height(results));
    fprintf('  Failed: %d datasets (empty folders or missing data)\n', n_total - height(results));
end

function [stats, mip_info] = process_single_dataset(dataset_path, params, qc_info)
    
    mip_info = [];  % Will be populated if MIPs are saved
    
    % Load stacks
    green_raw = load_stack(fullfile(dataset_path, 'green'));
    red_raw = load_stack(fullfile(dataset_path, 'red'));
    ctb_raw = load_stack(fullfile(dataset_path, 'blue'));
    
    [~, ~, nZ] = size(green_raw);
    
    % Get voxel size
    voxel_um = get_voxel_size(dataset_path);
    voxel_vol = prod(voxel_um);
    
    % Normalize
    green_norm = normalize_stack(green_raw, params.stretch_low, params.stretch_high);
    red_norm = normalize_stack(red_raw, params.stretch_low, params.stretch_high);
    ctb_norm = normalize_stack(ctb_raw, params.stretch_low, params.stretch_high);
    
    % Blur
    green_blur = apply_blur_3d(green_norm, params.blur_sigma);
    red_blur = apply_blur_3d(red_norm, params.blur_sigma);
    ctb_blur = apply_blur_3d(ctb_norm, params.blur_sigma);
    
    % Convert size limits to voxels
    green_min_vox = params.green_min_size_um3 / voxel_vol;
    red_min_vox = params.red_min_size_um3 / voxel_vol;
    red_max_vox = params.red_max_size_um3 / voxel_vol;
    ctb_min_vox = params.ctb_min_size_um3 / voxel_vol;
    ctb_max_vox = params.ctb_max_size_um3 / voxel_vol;
    
    % Determine thresholds based on method
    switch params.threshold_method
        case 'fixed'
            % Fixed thresholds are in 0-1 range, convert to 0-255
            green_thresh = params.green_threshold * 255;
            red_thresh = params.red_threshold * 255;
            ctb_thresh = params.ctb_threshold * 255;
            
        case 'percentile'
            % Threshold at Nth percentile of each image
            % This adapts to brightness - brighter images get higher thresholds
            green_vals = double(green_blur(green_blur > 0));
            red_vals = double(red_blur(red_blur > 0));
            ctb_vals = double(ctb_blur(ctb_blur > 0));
            
            green_thresh = prctile(green_vals, params.green_percentile);
            red_thresh = prctile(red_vals, params.red_percentile);
            ctb_thresh = prctile(ctb_vals, params.ctb_percentile);
            
        case 'otsu'
            % Compute Otsu threshold for each channel (returns 0-1, convert to 0-255)
            green_thresh = graythresh(green_blur) * 255;
            red_thresh = graythresh(red_blur) * 255;
            ctb_thresh = graythresh(ctb_blur) * 255;
            
        otherwise
            error('Unknown threshold_method: %s', params.threshold_method);
    end
    
    % Segment (thresholds are now in 0-255 range)
    green_mask = green_blur > green_thresh;
    green_mask = filter_by_size(green_mask, green_min_vox, inf);
    
    red_mask = red_blur > red_thresh;
    red_mask = filter_by_size(red_mask, red_min_vox, red_max_vox);
    
    ctb_mask = ctb_blur > ctb_thresh;
    ctb_mask = filter_by_size(ctb_mask, ctb_min_vox, ctb_max_vox);
    
    % Compute statistics
    n_green = nnz(green_mask);
    n_red_total = nnz(red_mask);
    n_ctb_total = nnz(ctb_mask);
    n_red_in_green = nnz(red_mask & green_mask);
    n_ctb_in_green = nnz(ctb_mask & green_mask);
    
    % Object counts
    CC_red_eng = bwconncomp(red_mask & green_mask, 26);
    CC_ctb_eng = bwconncomp(ctb_mask & green_mask, 26);
    CC_red_total = bwconncomp(red_mask, 26);
    
    % === SAVE QC OUTPUTS ===
    if qc_info.enabled
        save_qc_outputs(qc_info, green_norm, red_norm, ctb_norm, ...
            green_blur, red_blur, ctb_blur, ...
            green_mask, red_mask, ctb_mask, ...
            green_thresh, red_thresh, ctb_thresh, ...
            n_green, n_red_in_green, n_ctb_in_green, params);
    end
    
    % === SAVE MASKED IMAGE STACKS ===
    % These show only the signal that passed the filtering criteria
    if params.save_masked_images
        save_masked_images(dataset_path, ...
            green_raw, red_raw, ctb_raw, ...       % Raw intensities
            green_norm, red_norm, ctb_norm, ...   % Normalized intensities
            green_mask, red_mask, ctb_mask);
    end
    
    % === SAVE BINARY MASKS (v2) ===
    % These are the binary masks as 8-bit images (0=background, 255=mask)
    if params.save_binary_masks
        save_binary_masks(dataset_path, green_mask, red_mask, ctb_mask);
    end
    
    % === SAVE MIP OUTPUTS (v3) ===
    if params.save_mips
        mip_info = save_mip_outputs(dataset_path, ...
            green_raw, red_raw, ctb_raw, ...
            green_norm, red_norm, ctb_norm, ...
            green_mask, red_mask, ctb_mask, ...
            voxel_um, params);
        mip_info.nZ = nZ;
    end
    
    % Populate stats
    stats = struct();
    stats.voxel_vol_um3 = voxel_vol;
    
    % Store thresholds used
    stats.green_threshold_used = green_thresh;
    stats.red_threshold_used = red_thresh;
    stats.ctb_threshold_used = ctb_thresh;
    
    stats.green_volume_um3 = n_green * voxel_vol;
    stats.green_voxels = n_green;
    
    % Image volume for density calculations
    [H, W, nZ] = size(green_mask);
    image_volume_um3 = H * W * nZ * voxel_vol;
    stats.image_volume_um3 = image_volume_um3;
    
    % Red channel stats
    stats.red_total_objects = CC_red_total.NumObjects;
    stats.red_total_volume_um3 = n_red_total * voxel_vol;
    stats.red_engulfed_objects = CC_red_eng.NumObjects;
    stats.red_engulfed_volume_um3 = n_red_in_green * voxel_vol;
    stats.red_engulfment_pct = 100 * n_red_in_green / max(n_green, 1);
    
    % Mean volume per engulfed red object
    if CC_red_eng.NumObjects > 0
        stats.red_mean_object_volume_um3 = stats.red_engulfed_volume_um3 / CC_red_eng.NumObjects;
    else
        stats.red_mean_object_volume_um3 = 0;
    end
    
    % Red object density (objects per 1000 ÂµmÂ³ of image volume)
    stats.red_density_per_1000um3 = CC_red_total.NumObjects / (image_volume_um3 / 1000);
    
    % CTB channel stats
    stats.ctb_total_volume_um3 = n_ctb_total * voxel_vol;
    stats.ctb_engulfed_objects = CC_ctb_eng.NumObjects;
    stats.ctb_engulfed_volume_um3 = n_ctb_in_green * voxel_vol;
    stats.ctb_engulfment_pct = 100 * n_ctb_in_green / max(n_green, 1);
    
    % Mean volume per engulfed CTB object
    if CC_ctb_eng.NumObjects > 0
        stats.ctb_mean_object_volume_um3 = stats.ctb_engulfed_volume_um3 / CC_ctb_eng.NumObjects;
    else
        stats.ctb_mean_object_volume_um3 = 0;
    end
end

function save_qc_outputs(qc_info, green_norm, red_norm, ctb_norm, ...
    green_blur, red_blur, ctb_blur, ...
    green_mask, red_mask, ctb_mask, ...
    green_thresh, red_thresh, ctb_thresh, ...
    n_green, n_red_in_green, n_ctb_in_green, params)
    
    % Create QC output directory for this dataset
    qc_dir = qc_info.output_dir;
    if ~isfolder(qc_dir)
        mkdir(qc_dir);
    end
    
    % Create MIPs
    green_mip = max(green_norm, [], 3);
    red_mip = max(red_norm, [], 3);
    ctb_mip = max(ctb_norm, [], 3);
    
    green_mask_mip = max(green_mask, [], 3);
    red_mask_mip = max(red_mask, [], 3);
    ctb_mask_mip = max(ctb_mask, [], 3);
    
    % === SAVE OVERLAYS ===
    if qc_info.params.save_overlays
        % Green overlay
        save_mask_overlay_image(green_mip, green_mask_mip, [0 1 0], ...
            fullfile(qc_dir, 'green_overlay.png'), ...
            sprintf('Green: thresh=%.1f, %d voxels', green_thresh, nnz(green_mask)));
        
        % Red overlay
        save_mask_overlay_image(red_mip, red_mask_mip, [1 0 0], ...
            fullfile(qc_dir, 'red_overlay.png'), ...
            sprintf('Red: thresh=%.1f, engulf=%.3f%%', red_thresh, 100*n_red_in_green/max(n_green,1)));
        
        % CTB overlay
        save_mask_overlay_image(ctb_mip, ctb_mask_mip, [0 0.5 1], ...
            fullfile(qc_dir, 'ctb_overlay.png'), ...
            sprintf('CTB: thresh=%.1f, engulf=%.3f%%', ctb_thresh, 100*n_ctb_in_green/max(n_green,1)));
    end
    
    % === SAVE COMPOSITE ===
    if qc_info.params.save_composite
        fig = figure('Visible', 'off', 'Position', [100 100 1200 400]);
        
        % Raw composite
        subplot(1, 3, 1);
        rgb_raw = zeros([size(green_mip), 3], 'uint8');
        rgb_raw(:,:,1) = red_mip;
        rgb_raw(:,:,2) = green_mip;
        rgb_raw(:,:,3) = ctb_mip;
        imshow(rgb_raw);
        title('Raw MIP (R=Red, G=Green, B=CTB)', 'FontSize', 10);
        
        % Mask composite
        subplot(1, 3, 2);
        rgb_mask = zeros([size(green_mip), 3], 'uint8');
        rgb_mask(:,:,1) = uint8(red_mask_mip) * 255;
        rgb_mask(:,:,2) = uint8(green_mask_mip) * 255;
        rgb_mask(:,:,3) = uint8(ctb_mask_mip) * 255;
        imshow(rgb_mask);
        title('Masks (R=Red, G=Green, B=CTB)', 'FontSize', 10);
        
        % Overlay composite
        subplot(1, 3, 3);
        gray_bg = double(max(cat(3, green_mip, red_mip, ctb_mip), [], 3)) / 255;
        R = gray_bg; G = gray_bg; B = gray_bg;
        
        % Color the masks
        R(green_mask_mip) = 0.3; G(green_mask_mip) = 1; B(green_mask_mip) = 0.3;
        R(red_mask_mip) = 1; G(red_mask_mip) = G(red_mask_mip) * 0.3; B(red_mask_mip) = B(red_mask_mip) * 0.3;
        R(ctb_mask_mip) = R(ctb_mask_mip) * 0.3; G(ctb_mask_mip) = G(ctb_mask_mip) * 0.5; B(ctb_mask_mip) = 1;
        
        imshow(cat(3, R, G, B));
        title('Overlay', 'FontSize', 10);
        
        sgtitle(sprintf('%s\nRed eng: %.3f%%, CTB eng: %.3f%%', ...
            qc_info.dataset_name, 100*n_red_in_green/max(n_green,1), 100*n_ctb_in_green/max(n_green,1)), ...
            'FontSize', 12, 'FontWeight', 'bold');
        
        saveas(fig, fullfile(qc_dir, 'composite.png'));
        close(fig);
    end
    
    % === SAVE HISTOGRAMS ===
    if qc_info.params.save_histograms
        fig = figure('Visible', 'off', 'Position', [100 100 1200 350]);
        
        % Green histogram
        subplot(1, 3, 1);
        green_vals = double(green_blur(:));
        green_vals = green_vals(green_vals > 0);
        histogram(green_vals, 100, 'FaceColor', [0 0.7 0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        xline(green_thresh, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('Intensity');
        ylabel('Count');
        title(sprintf('Green: thresh=%.1f (%.1f%%ile)', green_thresh, params.green_percentile));
        xlim([0 255]);
        
        % Red histogram
        subplot(1, 3, 2);
        red_vals = double(red_blur(:));
        red_vals = red_vals(red_vals > 0);
        histogram(red_vals, 100, 'FaceColor', [0.8 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        xline(red_thresh, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('Intensity');
        ylabel('Count');
        title(sprintf('Red: thresh=%.1f (%.1f%%ile)', red_thresh, params.red_percentile));
        xlim([0 255]);
        
        % CTB histogram
        subplot(1, 3, 3);
        ctb_vals = double(ctb_blur(:));
        ctb_vals = ctb_vals(ctb_vals > 0);
        histogram(ctb_vals, 100, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        xline(ctb_thresh, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('Intensity');
        ylabel('Count');
        title(sprintf('CTB: thresh=%.1f (%.1f%%ile)', ctb_thresh, params.ctb_percentile));
        xlim([0 255]);
        
        sgtitle(sprintf('%s - Intensity Histograms', qc_info.dataset_name), 'FontSize', 12);
        
        saveas(fig, fullfile(qc_dir, 'histograms.png'));
        close(fig);
    end
    
    % === SAVE INFO TEXT FILE ===
    fid = fopen(fullfile(qc_dir, 'qc_info.txt'), 'w');
    fprintf(fid, 'QC Info for: %s\n', qc_info.dataset_name);
    fprintf(fid, '=====================================\n\n');
    fprintf(fid, 'Threshold method: %s\n\n', params.threshold_method);
    fprintf(fid, 'GREEN:\n');
    fprintf(fid, '  Percentile: %.1f\n', params.green_percentile);
    fprintf(fid, '  Threshold (0-255): %.2f\n', green_thresh);
    fprintf(fid, '  Mask voxels: %d\n\n', n_green);
    fprintf(fid, 'RED:\n');
    fprintf(fid, '  Percentile: %.1f\n', params.red_percentile);
    fprintf(fid, '  Threshold (0-255): %.2f\n', red_thresh);
    fprintf(fid, '  Size filter: %.2f - %.2f ÂµmÂ³\n', params.red_min_size_um3, params.red_max_size_um3);
    fprintf(fid, '  Voxels in green: %d\n', n_red_in_green);
    fprintf(fid, '  Engulfment: %.4f%%\n\n', 100*n_red_in_green/max(n_green,1));
    fprintf(fid, 'CTB:\n');
    fprintf(fid, '  Percentile: %.1f\n', params.ctb_percentile);
    fprintf(fid, '  Threshold (0-255): %.2f\n', ctb_thresh);
    fprintf(fid, '  Size filter: %.2f - %.2f ÂµmÂ³\n', params.ctb_min_size_um3, params.ctb_max_size_um3);
    fprintf(fid, '  Voxels in green: %d\n', n_ctb_in_green);
    fprintf(fid, '  Engulfment: %.4f%%\n', 100*n_ctb_in_green/max(n_green,1));
    fclose(fid);
end

function save_masked_images(dataset_path, ...
    green_raw, red_raw, ctb_raw, ...
    green_norm, red_norm, ctb_norm, ...
    green_mask, red_mask, ctb_mask)
    % Save both raw and normalized images with masks applied
    % Output shows only the signal that passed filtering criteria
    
    % Create output directories for RAW masked images
    green_raw_dir = fullfile(dataset_path, 'masked_green_raw');
    red_raw_dir = fullfile(dataset_path, 'masked_red_raw');
    ctb_raw_dir = fullfile(dataset_path, 'masked_ctb_raw');
    
    % Create output directories for NORMALIZED masked images
    green_norm_dir = fullfile(dataset_path, 'masked_green_norm');
    red_norm_dir = fullfile(dataset_path, 'masked_red_norm');
    ctb_norm_dir = fullfile(dataset_path, 'masked_ctb_norm');
    
    % Make all directories
    dirs = {green_raw_dir, red_raw_dir, ctb_raw_dir, ...
            green_norm_dir, red_norm_dir, ctb_norm_dir};
    for i = 1:numel(dirs)
        if ~isfolder(dirs{i}), mkdir(dirs{i}); end
    end
    
    % Get stack size
    [~, ~, nZ] = size(green_raw);
    
    % === SAVE RAW MASKED IMAGES ===
    % Green channel - raw masked
    green_masked_raw = green_raw .* uint8(green_mask);
    for z = 1:nZ
        imwrite(green_masked_raw(:,:,z), fullfile(green_raw_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % Red channel - raw masked
    red_masked_raw = red_raw .* uint8(red_mask);
    for z = 1:nZ
        imwrite(red_masked_raw(:,:,z), fullfile(red_raw_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % CTB channel - raw masked
    ctb_masked_raw = ctb_raw .* uint8(ctb_mask);
    for z = 1:nZ
        imwrite(ctb_masked_raw(:,:,z), fullfile(ctb_raw_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % === SAVE NORMALIZED MASKED IMAGES ===
    % Green channel - normalized masked
    green_masked_norm = green_norm .* uint8(green_mask);
    for z = 1:nZ
        imwrite(green_masked_norm(:,:,z), fullfile(green_norm_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % Red channel - normalized masked
    red_masked_norm = red_norm .* uint8(red_mask);
    for z = 1:nZ
        imwrite(red_masked_norm(:,:,z), fullfile(red_norm_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % CTB channel - normalized masked
    ctb_masked_norm = ctb_norm .* uint8(ctb_mask);
    for z = 1:nZ
        imwrite(ctb_masked_norm(:,:,z), fullfile(ctb_norm_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
end

function save_mask_overlay_image(mip, mask_mip, color, filepath, title_str)
    % Save a MIP with mask outline overlay
    fig = figure('Visible', 'off', 'Position', [100 100 600 600]);
    
    % Create RGB from grayscale MIP
    mip_norm = double(mip) / 255;
    R = mip_norm; G = mip_norm; B = mip_norm;
    
    % Create mask outline (edge of mask)
    outline = bwperim(mask_mip);
    
    % Dilate outline for visibility
    outline = imdilate(outline, strel('disk', 1));
    
    % Color the outline
    R(outline) = color(1);
    G(outline) = color(2);
    B(outline) = color(3);
    
    imshow(cat(3, R, G, B));
    title(title_str, 'FontSize', 11);
    
    saveas(fig, filepath);
    close(fig);
end

function generate_summary(results, output_dir, params)
    
    % Separate by group
    grp_A = results(strcmp(results.group, 'A'), :);
    grp_B = results(strcmp(results.group, 'B'), :);
    
    if isempty(grp_A) || isempty(grp_B)
        fprintf('Warning: Need both groups for comparison\n');
        return;
    end
    
    % === COMPUTE REPLICATE-LEVEL AVERAGES ===
    % Each replicate (biological sample) may have multiple imaging fields
    % Statistics should be computed on replicate means, not individual fields
    
    fprintf('\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('           COMPUTING REPLICATE-LEVEL AVERAGES\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    
    % Get unique replicates for each group
    reps_A = unique(grp_A.replicate);
    reps_B = unique(grp_B.replicate);
    
    fprintf('Group A: %d replicates (biological samples)\n', numel(reps_A));
    for i = 1:numel(reps_A)
        n_fields = sum(strcmp(grp_A.replicate, reps_A{i}));
        fprintf('  %s: %d imaging fields\n', reps_A{i}, n_fields);
    end
    
    fprintf('Group B: %d replicates (biological samples)\n', numel(reps_B));
    for i = 1:numel(reps_B)
        n_fields = sum(strcmp(grp_B.replicate, reps_B{i}));
        fprintf('  %s: %d imaging fields\n', reps_B{i}, n_fields);
    end
    
    % Compute replicate averages
    rep_means_A = compute_replicate_means(grp_A);
    rep_means_B = compute_replicate_means(grp_B);
    
    % Save replicate averages to file
    rep_results = [rep_means_A; rep_means_B];
    writetable(rep_results, fullfile(output_dir, 'replicate_averages.csv'));
    fprintf('\nReplicate averages saved to: replicate_averages.csv\n');
    
    % Open summary file
    fid = fopen(fullfile(output_dir, 'summary_report.txt'), 'w');
    
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(fid, '              ENGULFMENT ANALYSIS SUMMARY\n');
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    
    fprintf(fid, 'Parameters used:\n');
    fprintf(fid, '  Threshold method: %s\n', params.threshold_method);
    if strcmp(params.threshold_method, 'percentile')
        fprintf(fid, '  Green percentile: %.1f\n', params.green_percentile);
        fprintf(fid, '  Red percentile: %.1f\n', params.red_percentile);
        fprintf(fid, '  CTB percentile: %.1f\n', params.ctb_percentile);
    end
    fprintf(fid, '  Red size filter: %.2f - %.2f ÂµmÂ³\n', params.red_min_size_um3, params.red_max_size_um3);
    fprintf(fid, '  CTB size filter: %.2f - %.2f ÂµmÂ³\n\n', params.ctb_min_size_um3, params.ctb_max_size_um3);
    
    fprintf(fid, 'Sample structure:\n');
    fprintf(fid, '  Group A: %d biological replicates, %d imaging fields total\n', numel(reps_A), height(grp_A));
    fprintf(fid, '  Group B: %d biological replicates, %d imaging fields total\n\n', numel(reps_B), height(grp_B));
    
    fprintf(fid, 'IMPORTANT: Statistics below are computed on REPLICATE MEANS\n');
    fprintf(fid, '(each biological sample averaged across its imaging fields)\n\n');
    
    % Compare metrics using replicate averages
    metrics = {
        'red_engulfment_pct', 'ctb_engulfment_pct', ...
        'red_engulfed_objects', 'ctb_engulfed_objects', ...
        'red_engulfed_volume_um3', 'ctb_engulfed_volume_um3', ...
        'red_mean_object_volume_um3', 'ctb_mean_object_volume_um3', ...
        'red_density_per_1000um3'
    };
    metric_names = {
        'Red Engulfment (%)', 'CTB Engulfment (%)', ...
        'Red Engulfed Objects', 'CTB Engulfed Objects', ...
        'Red Engulfed Volume (ÂµmÂ³)', 'CTB Engulfed Volume (ÂµmÂ³)', ...
        'Red Mean Object Volume (ÂµmÂ³)', 'CTB Mean Object Volume (ÂµmÂ³)', ...
        'Red Object Density (per 1000 ÂµmÂ³)'
    };
    
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(fid, '          GROUP COMPARISONS (REPLICATE-LEVEL)\n');
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    
    for m = 1:numel(metrics)
        mn = metrics{m};
        
        vals_A = rep_means_A.(mn);
        vals_B = rep_means_B.(mn);
        
        n_A = numel(vals_A);
        n_B = numel(vals_B);
        mean_A = mean(vals_A);
        std_A = std(vals_A);
        mean_B = mean(vals_B);
        std_B = std(vals_B);
        
        % Two-sample t-test on replicate means
        [~, p, ~, stats] = ttest2(vals_A, vals_B);
        t_stat = stats.tstat;
        df = stats.df;
        
        % Cohen's d
        s_pooled = sqrt(((n_A-1)*std_A^2 + (n_B-1)*std_B^2) / (n_A+n_B-2));
        d = (mean_A - mean_B) / s_pooled;
        
        % % Difference (B relative to A)
        pct_diff = (mean_B - mean_A) / mean_A * 100;
        
        fprintf(fid, '%s:\n', metric_names{m});
        fprintf(fid, '  Group A (n=%d replicates): %.4f Â± %.4f\n', n_A, mean_A, std_A);
        fprintf(fid, '  Group B (n=%d replicates): %.4f Â± %.4f\n', n_B, mean_B, std_B);
        fprintf(fid, '  t(%g) = %.3f, p = %.4f\n', df, t_stat, p);
        fprintf(fid, '  Cohen''s d = %.3f\n', d);
        fprintf(fid, '  %% Difference (B vs A): %+.1f%%\n\n', pct_diff);
    end
    
    fclose(fid);
    
    % Create comparison plots (pass both individual data and replicate means)
    create_comparison_plots(grp_A, grp_B, rep_means_A, rep_means_B, output_dir, params);
    
    % === OUTLIER DETECTION ===
    detect_and_report_outliers(results, output_dir);
    
    fprintf('Summary saved to: %s\n', fullfile(output_dir, 'summary_report.txt'));
end

function rep_means = compute_replicate_means(grp_data)
    % Compute mean values for each biological replicate
    
    replicates = unique(grp_data.replicate);
    n_reps = numel(replicates);
    
    % Initialize output table
    rep_means = table();
    rep_means.replicate = replicates;
    rep_means.group = repmat(grp_data.group(1), n_reps, 1);
    rep_means.n_fields = zeros(n_reps, 1);
    
    % Engulfment percentages
    rep_means.red_engulfment_pct = zeros(n_reps, 1);
    rep_means.ctb_engulfment_pct = zeros(n_reps, 1);
    
    % Object counts
    rep_means.red_engulfed_objects = zeros(n_reps, 1);
    rep_means.ctb_engulfed_objects = zeros(n_reps, 1);
    
    % Volumes
    rep_means.green_volume_um3 = zeros(n_reps, 1);
    rep_means.red_engulfed_volume_um3 = zeros(n_reps, 1);
    rep_means.ctb_engulfed_volume_um3 = zeros(n_reps, 1);
    
    % Mean object volumes
    rep_means.red_mean_object_volume_um3 = zeros(n_reps, 1);
    rep_means.ctb_mean_object_volume_um3 = zeros(n_reps, 1);
    
    % Density
    rep_means.red_density_per_1000um3 = zeros(n_reps, 1);
    
    for i = 1:n_reps
        mask = strcmp(grp_data.replicate, replicates{i});
        rep_data = grp_data(mask, :);
        
        rep_means.n_fields(i) = height(rep_data);
        
        % Engulfment percentages
        rep_means.red_engulfment_pct(i) = mean(rep_data.red_engulfment_pct);
        rep_means.ctb_engulfment_pct(i) = mean(rep_data.ctb_engulfment_pct);
        
        % Object counts
        rep_means.red_engulfed_objects(i) = mean(rep_data.red_engulfed_objects);
        rep_means.ctb_engulfed_objects(i) = mean(rep_data.ctb_engulfed_objects);
        
        % Volumes
        rep_means.green_volume_um3(i) = mean(rep_data.green_volume_um3);
        rep_means.red_engulfed_volume_um3(i) = mean(rep_data.red_engulfed_volume_um3);
        rep_means.ctb_engulfed_volume_um3(i) = mean(rep_data.ctb_engulfed_volume_um3);
        
        % Mean object volumes
        rep_means.red_mean_object_volume_um3(i) = mean(rep_data.red_mean_object_volume_um3);
        rep_means.ctb_mean_object_volume_um3(i) = mean(rep_data.ctb_mean_object_volume_um3);
        
        % Density
        rep_means.red_density_per_1000um3(i) = mean(rep_data.red_density_per_1000um3);
    end
end

function detect_and_report_outliers(results, output_dir)
    % Detect datasets with unusually high engulfment (> mean + 2*SD)
    
    fprintf('\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('                    OUTLIER DETECTION\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    
    % Metrics to check for outliers
    metrics = {
        'red_engulfment_pct', 'ctb_engulfment_pct', ...
        'red_engulfed_objects', 'ctb_engulfed_objects', ...
        'red_engulfed_volume_um3', 'ctb_engulfed_volume_um3', ...
        'red_mean_object_volume_um3', 'ctb_mean_object_volume_um3', ...
        'red_density_per_1000um3'
    };
    metric_names = {
        'Red Engulfment (%)', 'CTB Engulfment (%)', ...
        'Red Engulfed Objects', 'CTB Engulfed Objects', ...
        'Red Engulfed Volume (ÂµmÂ³)', 'CTB Engulfed Volume (ÂµmÂ³)', ...
        'Red Mean Object Volume (ÂµmÂ³)', 'CTB Mean Object Volume (ÂµmÂ³)', ...
        'Red Object Density (per 1000 ÂµmÂ³)'
    };
    
    outlier_threshold = 2;  % Number of SDs above mean to flag as outlier
    
    all_outliers = {};
    
    % Open outlier report file
    fid = fopen(fullfile(output_dir, 'outliers_report.txt'), 'w');
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(fid, '              OUTLIER DETECTION REPORT\n');
    fprintf(fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(fid, 'Outliers defined as: value > mean + %.1f Ã— SD\n\n', outlier_threshold);
    
    for m = 1:numel(metrics)
        mn = metrics{m};
        vals = results.(mn);
        
        overall_mean = mean(vals);
        overall_std = std(vals);
        threshold = overall_mean + outlier_threshold * overall_std;
        
        fprintf('\n%s:\n', metric_names{m});
        fprintf('  Overall: %.4f Â± %.4f (threshold for outlier: > %.4f)\n', ...
            overall_mean, overall_std, threshold);
        
        fprintf(fid, '\n%s:\n', metric_names{m});
        fprintf(fid, '  Overall mean: %.4f\n', overall_mean);
        fprintf(fid, '  Overall SD: %.4f\n', overall_std);
        fprintf(fid, '  Outlier threshold (mean + %.1fÃ—SD): %.4f\n\n', outlier_threshold, threshold);
        
        % Find outliers
        outlier_idx = find(vals > threshold);
        
        if isempty(outlier_idx)
            fprintf('  No outliers detected.\n');
            fprintf(fid, '  No outliers detected.\n');
        else
            fprintf('  OUTLIERS DETECTED (%d):\n', numel(outlier_idx));
            fprintf(fid, '  OUTLIERS (%d datasets):\n', numel(outlier_idx));
            
            for i = 1:numel(outlier_idx)
                idx = outlier_idx(i);
                ds_name = results.dataset{idx};
                ds_group = results.group{idx};
                ds_val = vals(idx);
                ds_zscore = (ds_val - overall_mean) / overall_std;
                
                fprintf('    â€¢ %s (Group %s): %.4f (z=%.2f)\n', ds_name, ds_group, ds_val, ds_zscore);
                fprintf(fid, '    %s (Group %s)\n', ds_name, ds_group);
                fprintf(fid, '      Value: %.4f\n', ds_val);
                fprintf(fid, '      Z-score: %.2f\n', ds_zscore);
                fprintf(fid, '      Path: %s\n\n', results.path{idx});
                
                all_outliers{end+1} = struct('dataset', ds_name, 'group', ds_group, ...
                    'metric', metric_names{m}, 'value', ds_val, 'zscore', ds_zscore, ...
                    'path', results.path{idx}); %#ok<AGROW>
            end
        end
    end
    
    fclose(fid);
    
    % Summary
    fprintf('\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    if isempty(all_outliers)
        fprintf('No outliers detected in any metric.\n');
    else
        fprintf('Total outliers found: %d\n', numel(all_outliers));
        fprintf('See: %s\n', fullfile(output_dir, 'outliers_report.txt'));
        
        % Also save as CSV for easy filtering
        if ~isempty(all_outliers)
            outlier_table = struct2table([all_outliers{:}]);
            writetable(outlier_table, fullfile(output_dir, 'outliers.csv'));
            fprintf('Outlier list saved to: %s\n', fullfile(output_dir, 'outliers.csv'));
        end
    end
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
end

function create_comparison_plots(grp_A, grp_B, rep_means_A, rep_means_B, output_dir, params)
    
    figures_dir = fullfile(output_dir, 'figures');
    if ~isfolder(figures_dir), mkdir(figures_dir); end
    
    % Combine data for LMM
    all_data = [grp_A; grp_B];
    
    % Define all metrics to plot
    metrics = {
        'red_engulfment_pct', ...
        'ctb_engulfment_pct', ...
        'red_engulfed_objects', ...
        'ctb_engulfed_objects', ...
        'red_engulfed_volume_um3', ...
        'ctb_engulfed_volume_um3', ...
        'red_mean_object_volume_um3', ...
        'ctb_mean_object_volume_um3', ...
        'red_density_per_1000um3'
    };
    
    titles = {
        'Red (vGlut2) Engulfment', ...
        'CTB Engulfment', ...
        'Red Engulfed Objects', ...
        'CTB Engulfed Objects', ...
        'Red Engulfed Volume', ...
        'CTB Engulfed Volume', ...
        'Red Mean Object Volume', ...
        'CTB Mean Object Volume', ...
        'Red Object Density (whole image)'
    };
    
    ylabels = {
        '% of green volume', ...
        '% of green volume', ...
        'Number of objects', ...
        'Number of objects', ...
        'Volume (ÂµmÂ³)', ...
        'Volume (ÂµmÂ³)', ...
        'Volume per object (ÂµmÂ³)', ...
        'Volume per object (ÂµmÂ³)', ...
        'Objects per 1000 ÂµmÂ³'
    };
    
    % Print summary statistics to console and file
    fprintf('\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('      FINAL STATISTICS (LINEAR MIXED MODEL)\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    
    stats_fid = fopen(fullfile(output_dir, 'final_statistics.txt'), 'w');
    fprintf(stats_fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(stats_fid, '      FINAL STATISTICS (LINEAR MIXED MODEL)\n');
    fprintf(stats_fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    fprintf(stats_fid, 'Model: response ~ Group + (1|Animal)\n');
    fprintf(stats_fid, '  - Fixed effect: Group (A vs B)\n');
    fprintf(stats_fid, '  - Random effect: Animal (random intercept)\n');
    fprintf(stats_fid, '  - Observations: individual imaging fields\n\n');
    
    % Store results for methods section
    lmm_results = struct();
    
    for m = 1:numel(metrics)
        metric_name = metrics{m};
        
        % Individual imaging field values (for plotting)
        vals_A_fields = grp_A.(metric_name);
        vals_B_fields = grp_B.(metric_name);
        
        % Replicate mean values (for plotting means Â± SD)
        vals_A_rep = rep_means_A.(metric_name);
        vals_B_rep = rep_means_B.(metric_name);
        
        n_A_rep = numel(vals_A_rep);
        n_B_rep = numel(vals_B_rep);
        mean_A = mean(vals_A_rep);
        mean_B = mean(vals_B_rep);
        std_A = std(vals_A_rep);
        std_B = std(vals_B_rep);
        
        % === FIT LINEAR MIXED MODEL ===
        % Create table for fitlme
        lmm_tbl = table();
        lmm_tbl.response = all_data.(metric_name);
        lmm_tbl.group = categorical(all_data.group);
        lmm_tbl.animal = categorical(all_data.replicate);
        
        % Fit the model: response ~ group + (1|animal)
        try
            lme = fitlme(lmm_tbl, 'response ~ group + (1|animal)');
            
            % Extract fixed effects using anova or coefTest
            % Method that works across MATLAB versions
            [~, ~, FEstats] = fixedEffects(lme);
            
            % Group effect (B - A)
            % In MATLAB, categorical sorts alphabetically, so 'A' is reference
            % Second row is the group effect
            group_effect = FEstats.Estimate(2);  % Coefficient for group_B
            group_se = FEstats.SE(2);
            group_tstat = FEstats.tStat(2);
            group_pval = FEstats.pValue(2);
            group_df = FEstats.DF(2);
            
            % Variance components
            % covarianceParameters returns SD estimates
            [psi, ~, stats_re] = covarianceParameters(lme);
            
            % Random effect variance (animal)
            % psi{1} contains the covariance matrix for random effects
            % For random intercept only, this is a scalar variance
            if iscell(psi) && ~isempty(psi)
                var_animal = psi{1};  % This is already variance for (1|animal)
            else
                var_animal = 0;
            end
            var_residual = lme.MSE;  % Residual variance
            
            % Effect size: standardized coefficient
            % d = coefficient / sqrt(var_animal + var_residual)
            total_sd = sqrt(var_animal + var_residual);
            if total_sd > 0
                cohens_d = group_effect / total_sd;
            else
                cohens_d = 0;
            end
            
            % % Difference (B relative to A)
            if mean_A ~= 0
                pct_diff = (mean_B - mean_A) / mean_A * 100;
            else
                pct_diff = 0;
            end
            
            % Intraclass correlation (ICC) - proportion of variance due to animal
            total_var = var_animal + var_residual;
            if total_var > 0
                icc = var_animal / total_var;
            else
                icc = 0;
            end
            icc = max(0, min(1, icc));  % Ensure ICC is in [0, 1]
            
            % Post-hoc power approximation for LMM
            % Using design effect approach: effective n = n / (1 + (k-1)*ICC)
            % where k = average cluster size
            n_total = height(all_data);
            n_animals = numel(unique(all_data.replicate));
            avg_cluster_size = n_total / n_animals;
            design_effect = 1 + (avg_cluster_size - 1) * icc;
            design_effect = max(design_effect, 1);  % Ensure design effect >= 1
            effective_n_per_group = (n_total / 2) / design_effect;
            power = compute_ttest_power(abs(cohens_d), effective_n_per_group, effective_n_per_group, 0.05);
            
            lmm_success = true;
            
        catch ME
            % Fallback if LMM fails
            fprintf('  Warning: LMM failed for %s (%s). Using t-test.\n', metric_name, ME.message);
            [~, group_pval, ~, stats_t] = ttest2(vals_A_rep, vals_B_rep);
            group_tstat = stats_t.tstat;
            group_df = stats_t.df;
            group_effect = mean_B - mean_A;
            group_se = sqrt(std_A^2/n_A_rep + std_B^2/n_B_rep);
            s_pooled = sqrt(((n_A_rep-1)*std_A^2 + (n_B_rep-1)*std_B^2) / (n_A_rep+n_B_rep-2));
            cohens_d = group_effect / s_pooled;
            pct_diff = group_effect / mean_A * 100;
            icc = NaN;
            var_animal = NaN;
            var_residual = NaN;
            power = compute_ttest_power(abs(cohens_d), n_A_rep, n_B_rep, 0.05);
            lmm_success = false;
        end
        
        % Store results
        lmm_results(m).metric = metric_name;
        lmm_results(m).effect = group_effect;
        lmm_results(m).se = group_se;
        lmm_results(m).pval = group_pval;
        lmm_results(m).d = cohens_d;
        lmm_results(m).power = power;
        
        % Print to console
        fprintf('%s:\n', titles{m});
        fprintf('  Group A: n=%d replicates (%d fields), Mean Â± SD: %.4f Â± %.4f\n', ...
            n_A_rep, numel(vals_A_fields), mean_A, std_A);
        fprintf('  Group B: n=%d replicates (%d fields), Mean Â± SD: %.4f Â± %.4f\n', ...
            n_B_rep, numel(vals_B_fields), mean_B, std_B);
        if lmm_success
            fprintf('  LMM: Î² = %.4f Â± %.4f, t(%.1f) = %.3f, p = %.4f\n', ...
                group_effect, group_se, group_df, group_tstat, group_pval);
            fprintf('  ICC (animal): %.3f\n', icc);
        else
            fprintf('  t-test: t(%.1f) = %.3f, p = %.4f\n', group_df, group_tstat, group_pval);
        end
        fprintf('  Cohen''s d = %.3f, Power = %.1f%%\n', cohens_d, power * 100);
        fprintf('  %% Difference (B vs A): %+.1f%%\n', pct_diff);
        if power < 0.8
            fprintf('  âš  UNDERPOWERED (power < 80%%)\n');
        end
        fprintf('\n');
        
        % Print to file
        fprintf(stats_fid, '%s:\n', titles{m});
        fprintf(stats_fid, '  Group A: n=%d replicates (%d imaging fields)\n', n_A_rep, numel(vals_A_fields));
        fprintf(stats_fid, '           Mean Â± SD: %.4f Â± %.4f\n', mean_A, std_A);
        fprintf(stats_fid, '  Group B: n=%d replicates (%d imaging fields)\n', n_B_rep, numel(vals_B_fields));
        fprintf(stats_fid, '           Mean Â± SD: %.4f Â± %.4f\n', mean_B, std_B);
        if lmm_success
            fprintf(stats_fid, '  LMM Results:\n');
            fprintf(stats_fid, '    Fixed effect (B - A): Î² = %.4f Â± %.4f (SE)\n', group_effect, group_se);
            fprintf(stats_fid, '    t(%.1f) = %.3f, p = %.4f\n', group_df, group_tstat, group_pval);
            fprintf(stats_fid, '    Variance (animal): %.4f\n', var_animal);
            fprintf(stats_fid, '    Variance (residual): %.4f\n', var_residual);
            fprintf(stats_fid, '    ICC (animal): %.3f\n', icc);
        else
            fprintf(stats_fid, '  t-test (fallback): t(%.1f) = %.3f, p = %.4f\n', group_df, group_tstat, group_pval);
        end
        fprintf(stats_fid, '  Cohen''s d = %.3f\n', cohens_d);
        fprintf(stats_fid, '  Post-hoc power (Î±=0.05): %.1f%%\n', power * 100);
        fprintf(stats_fid, '  %% Difference (B vs A): %+.1f%%\n', pct_diff);
        if power < 0.8
            fprintf(stats_fid, '  âš  Study may be UNDERPOWERED (power < 80%%)\n');
        end
        fprintf(stats_fid, '\n');
        
        % === CREATE FIGURE ===
        fig = figure('Visible', 'off', 'Position', [100 100 500 650]);
        hold on;
        
        % Colors
        color_A = [0.5 0.5 0.5];         % Gray for individual points
        color_B = [1 0.5 0.5];           % Light red for individual points
        color_A_rep = [0 0 0];           % Black for replicate means
        color_B_rep = [0.8 0 0];         % Red for replicate means
        
        positions = [1 2];
        
        % Draw individual imaging field points (smaller)
        % NOTE: Using lighter colors instead of alpha for EPS vector compatibility
        swarm_x_A = beeswarm(vals_A_fields, positions(1), 0.25);
        swarm_x_B = beeswarm(vals_B_fields, positions(2), 0.25);
        
        % Simulate transparency with lighter colors (for vector EPS output)
        color_A_light = color_A + (1 - color_A) * 0.65;  % Blend with white
        color_B_light = color_B + (1 - color_B) * 0.65;
        
        scatter(swarm_x_A, vals_A_fields, 20, color_A_light, 'filled', ...
            'MarkerEdgeColor', 'none');
        scatter(swarm_x_B, vals_B_fields, 20, color_B_light, 'filled', ...
            'MarkerEdgeColor', 'none');
        
        % Draw replicate means as larger diamonds
        swarm_x_A_rep = beeswarm(vals_A_rep, positions(1), 0.12);
        swarm_x_B_rep = beeswarm(vals_B_rep, positions(2), 0.12);
        
        h_rep_A = scatter(swarm_x_A_rep, vals_A_rep, 100, color_A_rep, 'diamond', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        h_rep_B = scatter(swarm_x_B_rep, vals_B_rep, 100, color_B_rep, 'diamond', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        
        % Mean Â± SD error bars
        h_bar_A = errorbar(positions(1) + 0.35, mean_A, std_A, 'o', ...
            'Color', color_A_rep, 'MarkerFaceColor', color_A_rep, 'MarkerSize', 10, ...
            'LineWidth', 2, 'CapSize', 12);
        h_bar_B = errorbar(positions(2) + 0.35, mean_B, std_B, 'o', ...
            'Color', color_B_rep, 'MarkerFaceColor', color_B_rep, 'MarkerSize', 10, ...
            'LineWidth', 2, 'CapSize', 12);
        
        % Formatting
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Group A', 'Group B'});
        xlim([0.4 2.8]);
        ylabel(ylabels{m}, 'FontSize', 12);
        
        % Significance indicator
        if group_pval < 0.001
            sig_str = '***';
        elseif group_pval < 0.01
            sig_str = '**';
        elseif group_pval < 0.05
            sig_str = '*';
        else
            sig_str = 'ns';
        end
        
        title_str = sprintf('%s\nLMM: p = %.4f (%s), d = %.2f, Power = %.0f%%', ...
            titles{m}, group_pval, sig_str, cohens_d, power*100);
        title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
        
        % Add legend
        h_field = scatter(nan, nan, 20, color_A_light, 'filled');
        legend([h_field, h_rep_A, h_bar_A], ...
            {'Imaging fields', 'Animal means', 'Mean Â± SD'}, ...
            'Location', 'northwest', 'FontSize', 9);
        
        % Add significance bracket if significant
        all_vals = [vals_A_fields(:); vals_B_fields(:)];
        y_max = max(all_vals);
        y_min = min(all_vals);
        y_range = y_max - y_min;
        if y_range == 0, y_range = y_max * 0.1; end
        if y_range == 0, y_range = 1; end
        
        if group_pval < 0.05
            bracket_y = y_max + 0.12 * y_range;
            
            plot([1 1 2 2], [bracket_y bracket_y+0.03*y_range bracket_y+0.03*y_range bracket_y], ...
                'k-', 'LineWidth', 1.5);
            text(1.5, bracket_y + 0.06*y_range, sig_str, 'HorizontalAlignment', 'center', ...
                'FontSize', 14, 'FontWeight', 'bold');
            
            ylim([y_min - 0.05*y_range, bracket_y + 0.2*y_range]);
        else
            ylim([y_min - 0.05*y_range, y_max + 0.15*y_range]);
        end
        
        set(gca, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
        
        hold off;
        
        % Save figure
        saveas(fig, fullfile(figures_dir, sprintf('%s.png', metric_name)));
        
        % Save vector formats using painters renderer for true vector output
        % EPS (vector) - use print with painters renderer
        print(fig, fullfile(figures_dir, sprintf('%s.eps', metric_name)), '-depsc', '-painters');
        
        % PDF (vector) - often works better with Illustrator than EPS
        print(fig, fullfile(figures_dir, sprintf('%s.pdf', metric_name)), '-dpdf', '-painters');
        
        close(fig);
    end
    
    % === POWER ANALYSIS NOTES ===
    fprintf(stats_fid, '\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf(stats_fid, '                    POWER ANALYSIS NOTES\n');
    fprintf(stats_fid, 'â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');
    fprintf(stats_fid, 'Power approximation for LMM uses design effect approach:\n');
    fprintf(stats_fid, '  effective_n = n / (1 + (cluster_size - 1) Ã— ICC)\n');
    fprintf(stats_fid, 'Then standard two-sample power calculation on effective n.\n\n');
    fprintf(stats_fid, 'Interpretation:\n');
    fprintf(stats_fid, '  - Power â‰¥ 80%% is conventionally considered adequate\n');
    fprintf(stats_fid, '  - Power < 80%% suggests the study may be underpowered\n');
    fprintf(stats_fid, '  - ICC indicates proportion of variance due to animal\n');
    fprintf(stats_fid, '    (higher ICC = more clustering = less effective sample size)\n');
    
    fclose(stats_fid);
    fprintf('Statistics saved to: %s\n', fullfile(output_dir, 'final_statistics.txt'));
    fprintf('Figures saved to: %s\n', figures_dir);
    
    % === GENERATE METHODS SECTION ===
    n_A_reps = numel(unique(grp_A.replicate));
    n_B_reps = numel(unique(grp_B.replicate));
    n_A_fields = height(grp_A);
    n_B_fields = height(grp_B);
    save_methods_section(output_dir, params, n_A_reps, n_B_reps, n_A_fields, n_B_fields);
end

function save_methods_section(output_dir, params, n_A_reps, n_B_reps, n_A_fields, n_B_fields)
    % Save a methods section for publication
    
    methods_file = fullfile(output_dir, 'methods_section.txt');
    fid = fopen(methods_file, 'w');
    
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                           METHODS SECTION\n');
    fprintf(fid, '================================================================================\n\n');
    fprintf(fid, 'This text can be adapted for use in publications.\n\n');
    
    fprintf(fid, '--------------------------------------------------------------------------------\n');
    fprintf(fid, 'IMAGE PROCESSING\n');
    fprintf(fid, '--------------------------------------------------------------------------------\n\n');
    
    fprintf(fid, 'Confocal z-stack images were processed using custom MATLAB scripts\n');
    fprintf(fid, '(MathWorks, R2025a). For each imaging field, three channels were analyzed:\n');
    fprintf(fid, 'green (microglia/YFP), red (vGlut2), and blue (CTB).\n\n');
    
    fprintf(fid, 'Intensity normalization: Each channel was normalized independently using\n');
    fprintf(fid, 'linear contrast stretching. Pixel intensities were mapped from the %.1fth\n', params.stretch_low);
    fprintf(fid, 'to %.1fth percentile of non-zero values to the 0-255 range.\n\n', params.stretch_high);
    
    fprintf(fid, 'Spatial filtering: A 3D Gaussian blur (sigma = %.1f pixels) was applied to\n', params.blur_sigma);
    fprintf(fid, 'reduce noise while preserving object boundaries.\n\n');
    
    fprintf(fid, 'Thresholding: Adaptive percentile-based thresholding was used to account\n');
    fprintf(fid, 'for brightness variation across images. Pixels above the following\n');
    fprintf(fid, 'percentiles were considered foreground:\n');
    fprintf(fid, '  - Green (microglia): %.1fth percentile\n', params.green_percentile);
    fprintf(fid, '  - Red (vGlut2): %.1fth percentile\n', params.red_percentile);
    fprintf(fid, '  - CTB: %.1fth percentile\n\n', params.ctb_percentile);
    
    fprintf(fid, 'Object size filtering: Connected components were identified using 26-\n');
    fprintf(fid, 'connectivity in 3D. Objects were retained based on volume criteria:\n');
    fprintf(fid, '  - Green (microglia): minimum %.1f um^3\n', params.green_min_size_um3);
    fprintf(fid, '  - Red (vGlut2): %.2f - %.1f um^3\n', params.red_min_size_um3, params.red_max_size_um3);
    fprintf(fid, '  - CTB: %.2f - %.1f um^3\n\n', params.ctb_min_size_um3, params.ctb_max_size_um3);
    
    fprintf(fid, 'Engulfment quantification: Engulfed puncta were defined as red or CTB\n');
    fprintf(fid, 'voxels that overlapped spatially with the microglial (green) mask.\n');
    fprintf(fid, 'Engulfment was calculated as the percentage of microglial volume occupied\n');
    fprintf(fid, 'by engulfed puncta. Additional metrics included: number of engulfed\n');
    fprintf(fid, 'objects, total engulfed volume, mean volume per engulfed object, and\n');
    fprintf(fid, 'red channel object density across the entire imaging field.\n\n');
    
    fprintf(fid, '--------------------------------------------------------------------------------\n');
    fprintf(fid, 'STATISTICAL ANALYSIS\n');
    fprintf(fid, '--------------------------------------------------------------------------------\n\n');
    
    fprintf(fid, 'Statistical analyses were performed using MATLAB (MathWorks, R2025a).\n');
    fprintf(fid, 'Data were collected from %d animals in Group A (%d imaging fields) and\n', n_A_reps, n_A_fields);
    fprintf(fid, '%d animals in Group B (%d imaging fields).\n\n', n_B_reps, n_B_fields);
    
    fprintf(fid, 'To account for the hierarchical structure of the data (multiple imaging\n');
    fprintf(fid, 'fields nested within animals), linear mixed-effects models (LMM) were\n');
    fprintf(fid, 'used for all comparisons. Models were fit using restricted maximum\n');
    fprintf(fid, 'likelihood (REML) estimation with the following structure:\n\n');
    
    fprintf(fid, '    Response ~ Group + (1|Animal)\n\n');
    
    fprintf(fid, 'where Group (A vs B) was included as a fixed effect and Animal was\n');
    fprintf(fid, 'included as a random intercept to account for within-animal correlation\n');
    fprintf(fid, 'of observations. Degrees of freedom were estimated using the\n');
    fprintf(fid, 'Satterthwaite approximation.\n\n');
    
    fprintf(fid, 'Effect sizes were calculated as Cohen''s d, computed as the fixed effect\n');
    fprintf(fid, 'estimate (difference between groups) divided by the total standard\n');
    fprintf(fid, 'deviation (square root of the sum of between-animal and residual\n');
    fprintf(fid, 'variance components).\n\n');
    
    fprintf(fid, 'The intraclass correlation coefficient (ICC) was computed as:\n');
    fprintf(fid, '    ICC = variance_animal / (variance_animal + variance_residual)\n');
    fprintf(fid, 'representing the proportion of total variance attributable to\n');
    fprintf(fid, 'between-animal differences.\n\n');
    
    fprintf(fid, 'Post-hoc power was approximated using the design effect approach,\n');
    fprintf(fid, 'which accounts for within-animal clustering:\n');
    fprintf(fid, '    n_effective = n_observations / [1 + (cluster_size - 1) x ICC]\n');
    fprintf(fid, 'followed by standard two-sample power calculation.\n\n');
    
    fprintf(fid, 'Statistical significance was set at alpha = 0.05 (two-tailed).\n');
    fprintf(fid, 'Descriptive statistics are reported as mean +/- SD of animal-level\n');
    fprintf(fid, 'averages. LMM results include the fixed effect estimate (beta),\n');
    fprintf(fid, 'standard error, t-statistic, degrees of freedom, and p-value.\n\n');
    
    fprintf(fid, '--------------------------------------------------------------------------------\n');
    fprintf(fid, 'SOFTWARE\n');
    fprintf(fid, '--------------------------------------------------------------------------------\n\n');
    
    fprintf(fid, 'Image processing and statistical analysis were performed using custom\n');
    fprintf(fid, 'MATLAB scripts (version R2025a, MathWorks). Linear mixed-effects models\n');
    fprintf(fid, 'were fit using the fitlme function from the Statistics and Machine\n');
    fprintf(fid, 'Learning Toolbox. Code is available at [repository URL].\n\n');
    
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '================================================================================\n');
    
    fclose(fid);
    fprintf('Methods section saved to: %s\n', methods_file);
end

function power = compute_ttest_power(d, n1, n2, alpha)
    % Compute power for two-sample t-test
    % d = Cohen's d (absolute value)
    % n1, n2 = sample sizes
    % alpha = significance level (two-tailed)
    
    % Degrees of freedom
    df = n1 + n2 - 2;
    
    % Non-centrality parameter
    ncp = d * sqrt(n1 * n2 / (n1 + n2));
    
    % Critical t-value for alpha (two-tailed)
    t_crit = tinv(1 - alpha/2, df);
    
    % Power = P(|T| > t_crit | H1)
    % Under H1, T follows non-central t-distribution with ncp
    % Power = P(T > t_crit) + P(T < -t_crit) under non-central t
    
    % Using non-central t CDF
    power = 1 - nctcdf(t_crit, df, ncp) + nctcdf(-t_crit, df, ncp);
    
    % Ensure power is between 0 and 1
    power = max(0, min(1, power));
end

function n_needed = compute_sample_size_for_power(d, target_power, alpha)
    % Compute sample size per group needed to achieve target power
    % Uses iterative search
    
    if abs(d) < 0.001
        n_needed = Inf;  % Can't detect zero effect
        return;
    end
    
    % Search for n that gives target power
    for n = 2:500
        power = compute_ttest_power(abs(d), n, n, alpha);
        if power >= target_power
            n_needed = n;
            return;
        end
    end
    
    n_needed = Inf;  % Need more than 500 per group
end

function draw_violin(data, pos, width, color, alpha)
    % Draw a violin (kernel density estimate) at position pos
    if numel(data) < 4
        return;  % Need enough points for KDE
    end
    
    % Kernel density estimate
    [f, xi] = ksdensity(data, 'NumPoints', 100);
    
    % Normalize width
    f = f / max(f) * width;
    
    % Draw filled violin
    fill([pos - f, fliplr(pos + f)], [xi, fliplr(xi)], color, ...
        'FaceAlpha', alpha, 'EdgeColor', color, 'LineWidth', 1);
end

function draw_boxplot(data, pos, width, color)
    % Draw a simple box and whiskers
    q1 = prctile(data, 25);
    q2 = prctile(data, 50);  % median
    q3 = prctile(data, 75);
    iqr = q3 - q1;
    
    whisker_lo = max(min(data), q1 - 1.5*iqr);
    whisker_hi = min(max(data), q3 + 1.5*iqr);
    
    % Box
    rectangle('Position', [pos - width/2, q1, width, q3-q1], ...
        'EdgeColor', color, 'LineWidth', 1.5);
    
    % Median line
    plot([pos - width/2, pos + width/2], [q2, q2], '-', ...
        'Color', color, 'LineWidth', 2);
    
    % Whiskers
    plot([pos, pos], [whisker_lo, q1], '-', 'Color', color, 'LineWidth', 1.5);
    plot([pos, pos], [q3, whisker_hi], '-', 'Color', color, 'LineWidth', 1.5);
    
    % Whisker caps
    cap_width = width * 0.5;
    plot([pos - cap_width/2, pos + cap_width/2], [whisker_lo, whisker_lo], '-', 'Color', color, 'LineWidth', 1.5);
    plot([pos - cap_width/2, pos + cap_width/2], [whisker_hi, whisker_hi], '-', 'Color', color, 'LineWidth', 1.5);
end

function x_positions = beeswarm(data, center, width)
    % Create beeswarm (non-overlapping) x positions for points
    % Points are spread horizontally based on local density
    
    n = numel(data);
    x_positions = zeros(n, 1) + center;
    
    if n <= 1
        return;
    end
    
    % Sort by value
    [sorted_data, sort_idx] = sort(data);
    
    % Estimate point size in data units
    data_range = max(data) - min(data);
    if data_range == 0
        data_range = 1;
    end
    point_radius = data_range * 0.02;  % Approximate point size
    
    % Place points one by one
    placed_x = zeros(n, 1);
    placed_y = zeros(n, 1);
    
    for i = 1:n
        y = sorted_data(i);
        
        if i == 1
            placed_x(i) = center;
            placed_y(i) = y;
            continue;
        end
        
        % Find nearby already-placed points
        nearby_mask = abs(placed_y(1:i-1) - y) < point_radius * 2;
        nearby_x = placed_x(nearby_mask);
        
        if isempty(nearby_x)
            placed_x(i) = center;
        else
            % Find x position that doesn't overlap
            found = false;
            for offset = 0:0.02:width
                % Try right
                test_x = center + offset;
                if ~any(abs(nearby_x - test_x) < point_radius * 1.5)
                    placed_x(i) = test_x;
                    found = true;
                    break;
                end
                
                % Try left
                test_x = center - offset;
                if ~any(abs(nearby_x - test_x) < point_radius * 1.5)
                    placed_x(i) = test_x;
                    found = true;
                    break;
                end
            end
            
            if ~found
                placed_x(i) = center + (rand - 0.5) * width;
            end
        end
        
        placed_y(i) = y;
    end
    
    % Unsort to original order
    x_positions(sort_idx) = placed_x;
end

%% =================== HELPER FUNCTIONS ===================

function stack = load_stack(folder)
    files = dir(fullfile(folder, '*.tif'));
    if isempty(files)
        files = dir(fullfile(folder, '*.TIF'));
    end
    if isempty(files)
        error('No TIFF files in: %s', folder);
    end
    
    % Sort files
    names = {files.name};
    try
        names = natsortfiles(names);
    catch
        names = sort(names);
    end
    
    % Read first to get dimensions
    info = imfinfo(fullfile(folder, names{1}));
    H = info.Height; W = info.Width; nZ = numel(names);
    
    stack = zeros(H, W, nZ, 'uint8');
    for k = 1:nZ
        I = imread(fullfile(folder, names{k}));
        if ~isa(I, 'uint8'), I = im2uint8(I); end
        if ndims(I) == 3, I = I(:,:,1); end
        stack(:,:,k) = I;
    end
end

function out = normalize_stack(stack, lo_pct, hi_pct)
    % Per-image linear stretch
    vals = double(stack(:));
    vals_nz = vals(vals > 0);
    
    if isempty(vals_nz)
        out = stack;
        return;
    end
    
    lo = prctile(vals_nz, lo_pct);
    hi = prctile(vals_nz, hi_pct);
    
    if hi <= lo
        out = stack;
        return;
    end
    
    out = double(stack);
    out = (out - lo) / (hi - lo) * 255;
    out = uint8(min(max(out, 0), 255));
end

function out = apply_blur_3d(stack, sigma)
    if sigma <= 0
        out = stack;
        return;
    end
    
    ksz = max(3, 2*ceil(2*sigma)+1);
    G = fspecial('gaussian', ksz, sigma);
    
    out = zeros(size(stack), 'uint8');
    for k = 1:size(stack, 3)
        out(:,:,k) = imfilter(stack(:,:,k), G, 'replicate');
    end
end

function mask = filter_by_size(mask, min_vox, max_vox)
    CC = bwconncomp(mask, 26);
    if CC.NumObjects == 0
        return;
    end
    
    sizes = cellfun(@numel, CC.PixelIdxList);
    keep = (sizes >= min_vox) & (sizes <= max_vox);
    
    mask(:) = false;
    for i = find(keep)
        mask(CC.PixelIdxList{i}) = true;
    end
end

function voxel_um = get_voxel_size(dataset_path)
    % Extract voxel size from CZI file or TIFF metadata
    % Returns voxel size in MICRONS [x, y, z]
    % Errors if voxel size cannot be determined
    
    [parent_dir, dataset_name] = fileparts(dataset_path);
    [grandparent_dir, ~] = fileparts(parent_dir);
    
    % Try to find CZI file in multiple locations
    czi_locations = {
        fullfile(parent_dir, [dataset_name, '.czi']),           % Same folder as dataset
        fullfile(dataset_path, [dataset_name, '.czi']),         % Inside dataset folder
        fullfile(grandparent_dir, [dataset_name, '.czi'])       % Two levels up
    };
    
    for i = 1:numel(czi_locations)
        czi_file = czi_locations{i};
        if isfile(czi_file)
            try
                r = bfGetReader(czi_file);
                md = r.getMetadataStore();
                px = md.getPixelsPhysicalSizeX(0);
                py = md.getPixelsPhysicalSizeY(0);
                pz = md.getPixelsPhysicalSizeZ(0);
                if ~isempty(px) && ~isempty(py) && ~isempty(pz)
                    voxel_um = [px.value().doubleValue(), py.value().doubleValue(), pz.value().doubleValue()];
                    r.close();
                    fprintf('    Voxel size from CZI: [%.4f, %.4f, %.4f] Âµm\n', voxel_um);
                    return;
                end
                r.close();
            catch ME
                fprintf('    Warning: Could not read CZI metadata: %s\n', ME.message);
            end
        end
    end
    
    % Try to read from TIFF metadata
    green_dir = fullfile(dataset_path, 'green');
    if isfolder(green_dir)
        tif_files = dir(fullfile(green_dir, '*.tif'));
        if isempty(tif_files)
            tif_files = dir(fullfile(green_dir, '*.TIF'));
        end
        
        if ~isempty(tif_files)
            tif_path = fullfile(green_dir, tif_files(1).name);
            info = imfinfo(tif_path);
            
            voxel_um = [NaN, NaN, NaN];
            
            % Try XResolution/YResolution
            if isfield(info, 'XResolution') && info.XResolution > 0
                ru = '';
                if isfield(info, 'ResolutionUnit')
                    ru = lower(string(info.ResolutionUnit));
                end
                
                if contains(ru, 'centimeter')
                    voxel_um(1) = 1e4 / info.XResolution;  % cm to Âµm
                    voxel_um(2) = 1e4 / info.YResolution;
                elseif contains(ru, 'inch')
                    voxel_um(1) = 25400 / info.XResolution;  % inch to Âµm
                    voxel_um(2) = 25400 / info.YResolution;
                end
            end
            
            % Try ImageDescription for Z spacing
            if isfield(info, 'ImageDescription')
                desc = info.ImageDescription;
                % Try "spacing=X" format (ImageJ)
                m = regexp(desc, 'spacing=([0-9\.]+)', 'tokens', 'once');
                if ~isempty(m)
                    voxel_um(3) = str2double(m{1});
                end
                % Try other common formats
                if isnan(voxel_um(3))
                    m = regexp(desc, 'Z[:\s]*([0-9\.]+)\s*[uÂµ]m', 'tokens', 'once');
                    if ~isempty(m)
                        voxel_um(3) = str2double(m{1});
                    end
                end
            end
            
            if all(~isnan(voxel_um))
                fprintf('    Voxel size from TIFF: [%.4f, %.4f, %.4f] Âµm\n', voxel_um);
                return;
            end
        end
    end
    
    % Could not determine voxel size
    error('Could not determine voxel size for dataset: %s\nPlease ensure CZI file exists or TIFF metadata contains resolution info.', dataset_path);
end

function show_mask_overlay(img, mask, color)
    rgb = repmat(double(img)/255, [1 1 3]);
    
    for c = 1:3
        ch = rgb(:,:,c);
        ch(mask) = ch(mask) * 0.5 + color(c) * 0.5;
        rgb(:,:,c) = ch;
    end
    
    % Add edge
    edge_mask = edge(mask, 'Canny');
    for c = 1:3
        ch = rgb(:,:,c);
        ch(edge_mask) = color(c);
        rgb(:,:,c) = ch;
    end
    
    imshow(rgb);
end

function rgb = create_rgb_overlay(green_img, red_img, ctb_img, green_mask, red_mask, ctb_mask)
    % Create composite showing masks on original images
    
    R = double(red_img) / 255;
    G = double(green_img) / 255;
    B = double(ctb_img) / 255;
    
    % Highlight engulfed regions
    red_in_green = red_mask & green_mask;
    ctb_in_green = ctb_mask & green_mask;
    
    % Brighten engulfed puncta
    R(red_in_green) = min(R(red_in_green) + 0.5, 1);
    B(ctb_in_green) = min(B(ctb_in_green) + 0.5, 1);
    
    % Green outline
    green_edge = edge(green_mask, 'Canny');
    G(green_edge) = 1;
    R(green_edge) = 0;
    B(green_edge) = 0;
    
    rgb = cat(3, R, G, B);
end

function datasets = find_all_datasets(roots)
    datasets = struct('path', {}, 'name', {}, 'group', {}, 'replicate', {});
    
    for r = 1:numel(roots)
        root = roots{r};
        if ~isfolder(root), continue; end
        
        % Find Sample_A and Sample_B folders
        sample_dirs = dir(fullfile(root, 'Sample_*'));
        
        for s = 1:numel(sample_dirs)
            if ~sample_dirs(s).isdir, continue; end
            
            sample_path = fullfile(root, sample_dirs(s).name);
            
            % Determine group from sample folder name
            if contains(sample_dirs(s).name, 'Sample_A') || contains(sample_dirs(s).name, '_A')
                group = 'A';
            elseif contains(sample_dirs(s).name, 'Sample_B') || contains(sample_dirs(s).name, '_B')
                group = 'B';
            else
                continue;
            end
            
            % Recursively find all datasets under this sample folder
            found = find_datasets_recursive(sample_path, group, '');
            
            for i = 1:numel(found)
                datasets(end+1) = found(i); %#ok<AGROW>
            end
        end
    end
end

function datasets = find_datasets_recursive(folder, group, parent_name)
    % Recursively find all folders that contain red/green/blue subfolders WITH TIFFs
    % parent_name is used to track the replicate (biological sample) ID
    datasets = struct('path', {}, 'name', {}, 'group', {}, 'replicate', {});
    
    % Check if this folder is a dataset
    has_red = isfolder(fullfile(folder, 'red'));
    has_green = isfolder(fullfile(folder, 'green'));
    has_blue = isfolder(fullfile(folder, 'blue'));
    
    if has_red && has_green && has_blue
        % Also check that there are actual TIFF files
        green_tifs = dir(fullfile(folder, 'green', '*.tif'));
        if isempty(green_tifs)
            green_tifs = dir(fullfile(folder, 'green', '*.TIF'));
        end
        
        if ~isempty(green_tifs)
            [~, name] = fileparts(folder);
            datasets(1).path = folder;
            datasets(1).name = name;
            datasets(1).group = group;
            
            % Replicate is the parent folder (e.g., "Sample_1A" from path like Sample_A/Sample_1A/field_1)
            if isempty(parent_name)
                % If no parent tracked, extract from path
                path_parts = strsplit(folder, filesep);
                % Find the folder right after Sample_A or Sample_B
                replicate = name;  % Default to dataset name
                for p = 1:numel(path_parts)-1
                    if contains(path_parts{p}, 'Sample_A') || contains(path_parts{p}, 'Sample_B')
                        if p < numel(path_parts)
                            replicate = path_parts{p+1};
                            break;
                        end
                    end
                end
                datasets(1).replicate = replicate;
            else
                datasets(1).replicate = parent_name;
            end
        end
        return;  % Don't recurse into datasets (even if empty)
    end
    
    % Otherwise, recurse into subfolders
    contents = dir(folder);
    contents = contents([contents.isdir]);
    contents = contents(~ismember({contents.name}, {'.', '..'}));
    
    [~, current_name] = fileparts(folder);
    
    for i = 1:numel(contents)
        subfolder = fullfile(folder, contents(i).name);
        
        % Track replicate: the first level of subfolders under Sample_A/Sample_B
        % is the replicate (animal) ID
        if isempty(parent_name) && (startsWith(contents(i).name, 'Sample_') || ...
                contains(contents(i).name, '_A') || contains(contents(i).name, '_B'))
            % This is likely the replicate folder
            sub_datasets = find_datasets_recursive(subfolder, group, contents(i).name);
        else
            sub_datasets = find_datasets_recursive(subfolder, group, parent_name);
        end
        
        for j = 1:numel(sub_datasets)
            datasets(end+1) = sub_datasets(j); %#ok<AGROW>
        end
    end
end

%% =================== V2 FUNCTIONS - BINARY MASKS ===================

function save_binary_masks(dataset_path, green_mask, red_mask, ctb_mask)
    % Save binary masks as 8-bit images (0 = background, 255 = mask)
    % Also saves engulfed masks (puncta inside microglia)
    
    % Create output directories for individual channel masks
    green_mask_dir = fullfile(dataset_path, 'mask_green');
    red_mask_dir = fullfile(dataset_path, 'mask_red');
    ctb_mask_dir = fullfile(dataset_path, 'mask_ctb');
    
    % Create output directories for engulfed masks
    red_eng_mask_dir = fullfile(dataset_path, 'mask_red_engulfed');
    ctb_eng_mask_dir = fullfile(dataset_path, 'mask_ctb_engulfed');
    
    % Make all directories
    dirs = {green_mask_dir, red_mask_dir, ctb_mask_dir, ...
            red_eng_mask_dir, ctb_eng_mask_dir};
    for i = 1:numel(dirs)
        if ~isfolder(dirs{i}), mkdir(dirs{i}); end
    end
    
    % Compute engulfed masks (puncta inside microglia)
    red_engulfed_mask = red_mask & green_mask;
    ctb_engulfed_mask = ctb_mask & green_mask;
    
    % Get stack size
    [~, ~, nZ] = size(green_mask);
    
    % === SAVE GREEN MASK ===
    for z = 1:nZ
        imwrite(uint8(green_mask(:,:,z)) * 255, ...
            fullfile(green_mask_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % === SAVE RED MASK ===
    for z = 1:nZ
        imwrite(uint8(red_mask(:,:,z)) * 255, ...
            fullfile(red_mask_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % === SAVE CTB MASK ===
    for z = 1:nZ
        imwrite(uint8(ctb_mask(:,:,z)) * 255, ...
            fullfile(ctb_mask_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % === SAVE RED ENGULFED MASK ===
    for z = 1:nZ
        imwrite(uint8(red_engulfed_mask(:,:,z)) * 255, ...
            fullfile(red_eng_mask_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
    
    % === SAVE CTB ENGULFED MASK ===
    for z = 1:nZ
        imwrite(uint8(ctb_engulfed_mask(:,:,z)) * 255, ...
            fullfile(ctb_eng_mask_dir, sprintf('%04d.tif', z)), ...
            'Compression', 'none');
    end
end

%% =================== V3 FUNCTIONS - MIP OUTPUTS ===================

function mip_info = save_mip_outputs(dataset_path, ...
    green_raw, red_raw, ctb_raw, ...
    green_norm, red_norm, ctb_norm, ...
    green_mask, red_mask, ctb_mask, ...
    voxel_um, params)
    % Save maximum intensity projection images with various color schemes
    % Returns mip_info struct for montage generation
    
    % Create MIP output directory
    mip_dir = fullfile(dataset_path, 'Maximum_intensity_projections');
    if ~isfolder(mip_dir)
        mkdir(mip_dir);
    end
    
    % Calculate scalebar length in pixels
    scalebar_pixels = round(params.scalebar_um / voxel_um(1));
    
    % Get image dimensions
    [H, W, ~] = size(green_raw);
    
    % === CREATE MIPs ===
    % Raw MIPs
    green_mip_raw = max(green_raw, [], 3);
    red_mip_raw = max(red_raw, [], 3);
    ctb_mip_raw = max(ctb_raw, [], 3);
    
    % Normalized MIPs
    green_mip_norm = max(green_norm, [], 3);
    red_mip_norm = max(red_norm, [], 3);
    ctb_mip_norm = max(ctb_norm, [], 3);
    
    % Mask MIPs
    green_mask_mip = max(green_mask, [], 3);
    red_mask_mip = max(red_mask, [], 3);
    ctb_mask_mip = max(ctb_mask, [], 3);
    red_eng_mask_mip = max(red_mask & green_mask, [], 3);
    ctb_eng_mask_mip = max(ctb_mask & green_mask, [], 3);
    
    % Masked data MIPs (intensity where mask is true)
    green_masked_mip = max(green_norm .* uint8(green_mask), [], 3);
    red_masked_mip = max(red_norm .* uint8(red_mask), [], 3);
    ctb_masked_mip = max(ctb_norm .* uint8(ctb_mask), [], 3);
    
    % Engulfed-only MIPs (only show signal inside green mask)
    red_eng_mip = max(red_norm .* uint8(red_mask & green_mask), [], 3);
    ctb_eng_mip = max(ctb_norm .* uint8(ctb_mask & green_mask), [], 3);
    
    % === 1) SAVE MASK MIPs ===
    save_tiff_highres(uint8(green_mask_mip) * 255, fullfile(mip_dir, 'mask_green_MIP.tif'), params.mip_dpi);
    save_tiff_highres(uint8(red_mask_mip) * 255, fullfile(mip_dir, 'mask_red_MIP.tif'), params.mip_dpi);
    save_tiff_highres(uint8(ctb_mask_mip) * 255, fullfile(mip_dir, 'mask_ctb_MIP.tif'), params.mip_dpi);
    save_tiff_highres(uint8(red_eng_mask_mip) * 255, fullfile(mip_dir, 'mask_red_engulfed_MIP.tif'), params.mip_dpi);
    save_tiff_highres(uint8(ctb_eng_mask_mip) * 255, fullfile(mip_dir, 'mask_ctb_engulfed_MIP.tif'), params.mip_dpi);
    
    % === 2) SAVE RAW AND NORMALIZED RGB MIPs ===
    % Create RGB images - Original: R=red, G=green, B=ctb
    raw_rgb_original = create_rgb_original(green_mip_raw, red_mip_raw, ctb_mip_raw);
    norm_rgb_original = create_rgb_original(green_mip_norm, red_mip_norm, ctb_mip_norm);
    
    % Create RGB images - Adjusted: grey=green, magenta=red, green=ctb
    raw_rgb_adjusted = create_rgb_adjusted(green_mip_raw, red_mip_raw, ctb_mip_raw);
    norm_rgb_adjusted = create_rgb_adjusted(green_mip_norm, red_mip_norm, ctb_mip_norm);
    
    % Save without scalebar
    save_tiff_highres(raw_rgb_original, fullfile(mip_dir, 'raw_RGB_original.tif'), params.mip_dpi);
    save_tiff_highres(raw_rgb_adjusted, fullfile(mip_dir, 'raw_RGB_adjusted.tif'), params.mip_dpi);
    save_tiff_highres(norm_rgb_original, fullfile(mip_dir, 'norm_RGB_original.tif'), params.mip_dpi);
    save_tiff_highres(norm_rgb_adjusted, fullfile(mip_dir, 'norm_RGB_adjusted.tif'), params.mip_dpi);
    
    % Save with scalebar
    save_tiff_highres(add_scalebar(raw_rgb_original, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'raw_RGB_original_scalebar.tif'), params.mip_dpi);
    save_tiff_highres(add_scalebar(raw_rgb_adjusted, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'raw_RGB_adjusted_scalebar.tif'), params.mip_dpi);
    save_tiff_highres(add_scalebar(norm_rgb_original, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'norm_RGB_original_scalebar.tif'), params.mip_dpi);
    save_tiff_highres(add_scalebar(norm_rgb_adjusted, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'norm_RGB_adjusted_scalebar.tif'), params.mip_dpi);
    
    % === 3) SAVE MASKED RGB MIPs ===
    masked_rgb_original = create_rgb_original(green_masked_mip, red_masked_mip, ctb_masked_mip);
    masked_rgb_adjusted = create_rgb_adjusted(green_masked_mip, red_masked_mip, ctb_masked_mip);
    
    % Save without scalebar
    save_tiff_highres(masked_rgb_original, fullfile(mip_dir, 'masked_RGB_original.tif'), params.mip_dpi);
    save_tiff_highres(masked_rgb_adjusted, fullfile(mip_dir, 'masked_RGB_adjusted.tif'), params.mip_dpi);
    
    % Save with scalebar
    save_tiff_highres(add_scalebar(masked_rgb_original, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'masked_RGB_original_scalebar.tif'), params.mip_dpi);
    save_tiff_highres(add_scalebar(masked_rgb_adjusted, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'masked_RGB_adjusted_scalebar.tif'), params.mip_dpi);
    
    % === 4) SAVE ENGULFED-ONLY RGB MIPs ===
    % Show glial mask + engulfed puncta only
    engulfed_rgb_original = create_rgb_original(green_masked_mip, red_eng_mip, ctb_eng_mip);
    engulfed_rgb_adjusted = create_rgb_adjusted(green_masked_mip, red_eng_mip, ctb_eng_mip);
    
    % Save without scalebar
    save_tiff_highres(engulfed_rgb_original, fullfile(mip_dir, 'engulfed_RGB_original.tif'), params.mip_dpi);
    save_tiff_highres(engulfed_rgb_adjusted, fullfile(mip_dir, 'engulfed_RGB_adjusted.tif'), params.mip_dpi);
    
    % Save with scalebar
    save_tiff_highres(add_scalebar(engulfed_rgb_original, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'engulfed_RGB_original_scalebar.tif'), params.mip_dpi);
    save_tiff_highres(add_scalebar(engulfed_rgb_adjusted, scalebar_pixels, H, W), ...
        fullfile(mip_dir, 'engulfed_RGB_adjusted_scalebar.tif'), params.mip_dpi);
    
    % === Return data for montages (v4 - additional fields) ===
    mip_info = struct();
    mip_info.mip_rgb_original = norm_rgb_original;   % Original color RGB for montage
    mip_info.mip_rgb_adjusted = norm_rgb_adjusted;   % Color-adjusted RGB for montage
    mip_info.mip_green_grey = green_mip_norm;        % Grey-scale green for montage
    mip_info.mip_green_mask = uint8(green_mask_mip) * 255;  % Binary mask for montage
end

function rgb = create_rgb_original(green_mip, red_mip, ctb_mip)
    % Original color scheme: R=red, G=green, B=ctb
    rgb = zeros([size(green_mip), 3], 'uint8');
    rgb(:,:,1) = red_mip;
    rgb(:,:,2) = green_mip;
    rgb(:,:,3) = ctb_mip;
end

function rgb = create_rgb_adjusted(green_mip, red_mip, ctb_mip)
    % Adjusted color scheme: grey=green (microglia), magenta=red (vGlut2), green=ctb
    % Grey = equal R, G, B
    % Magenta = R + B
    % Green = G only
    
    green_d = double(green_mip);
    red_d = double(red_mip);
    ctb_d = double(ctb_mip);
    
    % Combine: grey base from green, magenta from red, green from ctb
    R = green_d + red_d;           % Grey + Magenta red component
    G = green_d + ctb_d;           % Grey + CTB green component
    B = green_d + red_d;           % Grey + Magenta blue component
    
    % Normalize to prevent overflow
    max_val = max([R(:); G(:); B(:)]);
    if max_val > 255
        R = R / max_val * 255;
        G = G / max_val * 255;
        B = B / max_val * 255;
    end
    
    rgb = zeros([size(green_mip), 3], 'uint8');
    rgb(:,:,1) = uint8(R);
    rgb(:,:,2) = uint8(G);
    rgb(:,:,3) = uint8(B);
end

function img_out = add_scalebar(img, scalebar_pixels, H, W)
    % Add a white scalebar to bottom-right corner of image
    % Scalebar is 3 pixels thick with 10 pixel margin from edges
    
    img_out = img;
    
    margin = 10;
    thickness = 3;
    
    % Calculate scalebar position
    x_end = W - margin;
    x_start = x_end - scalebar_pixels + 1;
    y_start = H - margin - thickness + 1;
    y_end = H - margin;
    
    % Ensure within bounds
    x_start = max(1, x_start);
    y_start = max(1, y_start);
    
    % Draw white scalebar
    if ndims(img_out) == 3
        img_out(y_start:y_end, x_start:x_end, :) = 255;
    else
        img_out(y_start:y_end, x_start:x_end) = 255;
    end
end

function save_tiff_highres(img, filepath, dpi)
    % Save TIFF with specified resolution
    % dpi = dots per inch
    
    imwrite(img, filepath, 'tif', 'Compression', 'none', ...
        'Resolution', [dpi dpi]);
end

function slide_info = extract_slide_info(dataset_name)
    % Extract slide/section info from dataset name
    % Looks for patterns like "Section_1", "Slide_2", etc.
    
    % Try to find Section_X pattern
    match = regexp(dataset_name, '[Ss]ection[_\s]*(\d+)', 'tokens', 'once');
    if ~isempty(match)
        slide_info = sprintf('Section_%s', match{1});
        return;
    end
    
    % Try to find Slide_X pattern
    match = regexp(dataset_name, '[Ss]lide[_\s]*(\d+)', 'tokens', 'once');
    if ~isempty(match)
        slide_info = sprintf('Slide_%s', match{1});
        return;
    end
    
    % Try to find _X_ pattern at start (like "1_side_1")
    match = regexp(dataset_name, '^(\d+)[_\s]', 'tokens', 'once');
    if ~isempty(match)
        slide_info = sprintf('Section_%s', match{1});
        return;
    end
    
    % Default: use "All" to group everything together
    slide_info = 'All';
end

%% =================== V4 FUNCTIONS - MONTAGE GENERATION ===================

function generate_montages(montage_data, root_path, params)
    % Generate montage images for each sample group (Sample1A-7A, Sample1B-7B)
    % montage_data: struct array with fields: group, replicate, name, nZ, 
    %               mip_rgb_original, mip_rgb_adjusted, mip_green_grey, mip_green_mask
    
    % Get unique groups (A and B)
    groups = unique({montage_data.group});
    
    for g = 1:numel(groups)
        grp = groups{g};
        grp_mask = strcmp({montage_data.group}, grp);
        grp_data = montage_data(grp_mask);
        
        if isempty(grp_data)
            continue;
        end
        
        % Find output directory for this group (Sample_A or Sample_B)
        sample_dirs = dir(fullfile(root_path, sprintf('Sample_%s*', grp)));
        sample_dirs = sample_dirs([sample_dirs.isdir]);
        % Filter to just Sample_A or Sample_B (not Sample_1A etc)
        main_sample_idx = [];
        for sd = 1:numel(sample_dirs)
            if strcmp(sample_dirs(sd).name, sprintf('Sample_%s', grp))
                main_sample_idx = sd;
                break;
            end
        end
        
        if ~isempty(main_sample_idx)
            montage_output_dir = fullfile(root_path, sample_dirs(main_sample_idx).name, 'Montages');
        else
            montage_output_dir = fullfile(root_path, sprintf('Sample_%s', grp), 'Montages');
        end
        
        if ~isfolder(montage_output_dir)
            mkdir(montage_output_dir);
        end
        
        % Get unique replicates within this group (Sample1A, Sample2A, etc.)
        replicates = unique({grp_data.replicate});
        
        for r = 1:numel(replicates)
            rep = replicates{r};
            rep_mask = strcmp({grp_data.replicate}, rep);
            rep_data = grp_data(rep_mask);
            
            if isempty(rep_data)
                continue;
            end
            
            fprintf('  Creating montages for %s (%d cells)...\n', rep, numel(rep_data));
            
            % Create montages for this replicate (Sample1A, etc.)
            create_replicate_montages(rep_data, montage_output_dir, rep, params);
        end
    end
end

function create_replicate_montages(rep_data, output_dir, replicate_name, params)
    % Create 4 montage images for a single replicate (e.g., Sample1A)
    % 1. RGB Original MIP
    % 2. RGB Adjusted MIP
    % 3. Green channel MIP (grayscale)
    % 4. Green channel binary mask MIP
    
    n_cells = numel(rep_data);
    
    if n_cells == 0
        return;
    end
    
    % Determine grid size (aim for roughly square, max ~6 columns for readability)
    max_cols = 6;
    n_cols = min(n_cells, max_cols);
    n_rows = ceil(n_cells / n_cols);
    
    % Get image size from first cell
    [cell_h, cell_w, ~] = size(rep_data(1).mip_rgb_adjusted);
    
    % Calculate label area height - need space for filename + section count
    % Estimate ~20 pixels per line, allow up to 3 lines for wrapped text
    label_height = 60;
    cell_h_with_label = cell_h + label_height;
    
    % Create montage canvases
    montage_h = n_rows * cell_h_with_label;
    montage_w = n_cols * cell_w;
    
    rgb_original_montage = zeros(montage_h, montage_w, 3, 'uint8');
    rgb_adjusted_montage = zeros(montage_h, montage_w, 3, 'uint8');
    green_grey_montage = zeros(montage_h, montage_w, 'uint8');
    green_mask_montage = zeros(montage_h, montage_w, 'uint8');
    
    % Place images and collect label info
    label_info = struct('x', {}, 'y', {}, 'text', {});
    
    for i = 1:n_cells
        % Calculate position in grid
        row = ceil(i / n_cols);
        col = mod(i - 1, n_cols) + 1;
        
        % Label area is at the TOP, image below
        y_label_start = (row - 1) * cell_h_with_label + 1;
        y_label_end = y_label_start + label_height - 1;
        y_img_start = y_label_end + 1;
        y_img_end = y_img_start + cell_h - 1;
        x_start = (col - 1) * cell_w + 1;
        x_end = x_start + cell_w - 1;
        
        % Get cell data
        cell_rgb_orig = rep_data(i).mip_rgb_original;
        cell_rgb_adj = rep_data(i).mip_rgb_adjusted;
        cell_green_grey = rep_data(i).mip_green_grey;
        cell_green_mask = rep_data(i).mip_green_mask;
        cell_nZ = rep_data(i).nZ;
        cell_name = rep_data(i).name;
        
        % Resize if necessary (in case of size mismatch)
        if size(cell_rgb_orig, 1) ~= cell_h || size(cell_rgb_orig, 2) ~= cell_w
            cell_rgb_orig = imresize(cell_rgb_orig, [cell_h, cell_w]);
        end
        if size(cell_rgb_adj, 1) ~= cell_h || size(cell_rgb_adj, 2) ~= cell_w
            cell_rgb_adj = imresize(cell_rgb_adj, [cell_h, cell_w]);
        end
        if size(cell_green_grey, 1) ~= cell_h || size(cell_green_grey, 2) ~= cell_w
            cell_green_grey = imresize(cell_green_grey, [cell_h, cell_w]);
        end
        if size(cell_green_mask, 1) ~= cell_h || size(cell_green_mask, 2) ~= cell_w
            cell_green_mask = imresize(cell_green_mask, [cell_h, cell_w], 'nearest');
        end
        
        % Set label area to dark grey background
        rgb_original_montage(y_label_start:y_label_end, x_start:x_end, :) = 40;
        rgb_adjusted_montage(y_label_start:y_label_end, x_start:x_end, :) = 40;
        green_grey_montage(y_label_start:y_label_end, x_start:x_end) = 40;
        green_mask_montage(y_label_start:y_label_end, x_start:x_end) = 40;
        
        % Place images
        rgb_original_montage(y_img_start:y_img_end, x_start:x_end, :) = cell_rgb_orig;
        rgb_adjusted_montage(y_img_start:y_img_end, x_start:x_end, :) = cell_rgb_adj;
        green_grey_montage(y_img_start:y_img_end, x_start:x_end) = cell_green_grey;
        green_mask_montage(y_img_start:y_img_end, x_start:x_end) = cell_green_mask;
        
        % Store label info for later text rendering
        label_info(i).x = x_start + 5;
        label_info(i).y = y_label_start + 5;
        label_info(i).text = sprintf('%s\nn=%d sections', cell_name, cell_nZ);
        label_info(i).cell_w = cell_w;
    end
    
    % Add text labels using insertText if available
    try
        % Check if insertText is available (Computer Vision Toolbox)
        rgb_original_montage = add_wrapped_labels(rgb_original_montage, label_info, cell_w);
        rgb_adjusted_montage = add_wrapped_labels(rgb_adjusted_montage, label_info, cell_w);
        
        % For grayscale images, convert to RGB, add labels, convert back
        green_grey_rgb = repmat(green_grey_montage, [1 1 3]);
        green_grey_rgb = add_wrapped_labels(green_grey_rgb, label_info, cell_w);
        green_grey_montage = green_grey_rgb(:,:,1);
        
        green_mask_rgb = repmat(green_mask_montage, [1 1 3]);
        green_mask_rgb = add_wrapped_labels(green_mask_rgb, label_info, cell_w);
        green_mask_montage = green_mask_rgb(:,:,1);
    catch ME
        % insertText not available - labels will just be dark bars
        fprintf('    Warning: Could not add text labels (%s)\n', ME.message);
    end
    
    % Safe replicate name for filenames
    safe_rep = regexprep(replicate_name, '[^\w]', '_');
    
    % Save all 4 montages
    rgb_orig_filename = fullfile(output_dir, sprintf('montage_%s_RGB_original.tif', safe_rep));
    rgb_adj_filename = fullfile(output_dir, sprintf('montage_%s_RGB_adjusted.tif', safe_rep));
    green_grey_filename = fullfile(output_dir, sprintf('montage_%s_green_MIP.tif', safe_rep));
    green_mask_filename = fullfile(output_dir, sprintf('montage_%s_green_mask.tif', safe_rep));
    
    save_tiff_highres(rgb_original_montage, rgb_orig_filename, params.mip_dpi);
    save_tiff_highres(rgb_adjusted_montage, rgb_adj_filename, params.mip_dpi);
    save_tiff_highres(green_grey_montage, green_grey_filename, params.mip_dpi);
    save_tiff_highres(green_mask_montage, green_mask_filename, params.mip_dpi);
    
    fprintf('    Saved: %s\n', rgb_orig_filename);
    fprintf('    Saved: %s\n', rgb_adj_filename);
    fprintf('    Saved: %s\n', green_grey_filename);
    fprintf('    Saved: %s\n', green_mask_filename);
end

function img_out = add_wrapped_labels(img, label_info, max_width)
    % Add wrapped text labels to image using insertText
    % max_width is the available width for text in pixels
    
    img_out = img;
    
    % Approximate characters per line based on font size and cell width
    font_size = 12;
    chars_per_pixel = 0.15;  % Rough estimate for font size 12
    max_chars = floor(max_width * chars_per_pixel);
    
    for i = 1:numel(label_info)
        % Parse the label text (filename + section count)
        full_text = label_info(i).text;
        lines = strsplit(full_text, '\n');
        
        wrapped_lines = {};
        for ln = 1:numel(lines)
            line = lines{ln};
            % Wrap this line if too long
            while length(line) > max_chars
                % Find a good break point
                break_idx = max_chars;
                % Try to break at underscore or space
                space_idx = find(line(1:max_chars) == ' ' | line(1:max_chars) == '_', 1, 'last');
                if ~isempty(space_idx) && space_idx > max_chars/2
                    break_idx = space_idx;
                end
                wrapped_lines{end+1} = line(1:break_idx); %#ok<AGROW>
                line = line(break_idx+1:end);
            end
            if ~isempty(line)
                wrapped_lines{end+1} = line; %#ok<AGROW>
            end
        end
        
        % Add each line
        y_pos = label_info(i).y;
        x_pos = label_info(i).x;
        line_height = font_size + 4;
        
        for ln = 1:numel(wrapped_lines)
            try
                img_out = insertText(img_out, [x_pos, y_pos], wrapped_lines{ln}, ...
                    'FontSize', font_size, 'TextColor', 'white', 'BoxOpacity', 0);
            catch
                % If insertText fails, just continue
            end
            y_pos = y_pos + line_height;
        end
    end
end

function T = rename_channel_columns(T)
    % Rename channel result columns to a color+marker convention for clarity in
    % the exported CSV: red=VGluT2, green=microglia, blue=CTB. Only columns that
    % are present are renamed, so this is safe if the schema changes.
    old_new = {
        'green_threshold_used',        'green_microglia_threshold_used'
        'green_volume_um3',            'green_microglia_volume_um3'
        'green_voxels',                'green_microglia_voxels'
        'red_threshold_used',          'red_VGluT2_threshold_used'
        'red_total_objects',           'red_VGluT2_total_objects'
        'red_total_volume_um3',        'red_VGluT2_total_volume_um3'
        'red_engulfed_objects',        'red_VGluT2_engulfed_objects'
        'red_engulfed_volume_um3',     'red_VGluT2_engulfed_volume_um3'
        'red_engulfment_pct',          'red_VGluT2_engulfment_pct'
        'red_mean_object_volume_um3',  'red_VGluT2_mean_object_volume_um3'
        'red_density_per_1000um3',     'red_VGluT2_density_per_1000um3'
        'ctb_threshold_used',          'blue_CTB_threshold_used'
        'ctb_total_volume_um3',        'blue_CTB_total_volume_um3'
        'ctb_engulfed_objects',        'blue_CTB_engulfed_objects'
        'ctb_engulfed_volume_um3',     'blue_CTB_engulfed_volume_um3'
        'ctb_engulfment_pct',          'blue_CTB_engulfment_pct'
        'ctb_mean_object_volume_um3',  'blue_CTB_mean_object_volume_um3'
    };
    vn = T.Properties.VariableNames;
    for i = 1:size(old_new, 1)
        idx = strcmp(vn, old_new{i, 1});
        if any(idx)
            vn{idx} = old_new{i, 2};
        end
    end
    T.Properties.VariableNames = vn;
end
