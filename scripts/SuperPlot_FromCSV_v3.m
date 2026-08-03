% SuperPlot_FromCSV_v3.m
% Standalone script to generate publication-ready SuperPlot figures
% from pre-computed batch results CSV
%
% Generates side-by-side: SuperPlot | Box plot (animal means) | Mean±SD
% Two color scheme options included
%
% v3 changes:
%   - Much wider feathering with force-directed positioning
%   - Fixed connecting lines to use same positions as scatter points
%   - Pre-compute all positions before drawing
%
% Output: Vector EPS and PDF files suitable for Illustrator

clear; clc;

%% ===================== USER CONFIGURATION =====================

% Path to your batch results CSV
CSV_FILE = 'Y:\Colenso\2026\Global_PerSample\BatchResults\engulfment_results.csv';

% Output directory for figures
OUTPUT_DIR = 'Y:\Colenso\2026\Global_PerSample\BatchResults\figures_superplot';

% Metrics to plot (column names from CSV)
METRICS = {'red_VGluT2_engulfment_pct', 'blue_CTB_engulfment_pct'};

% Display titles and Y-axis labels for each metric
TITLES = {'vGlut2 Engulfment', 'CTB Engulfment'};
YLABELS = {'% of microglial volume', '% of microglial volume'};

% Color scheme: 'unique' = 14 distinct colors, 'family' = gray/red families
COLOR_SCHEME = 'unique';  % Change to 'unique' for distinct animal colors

% Figure sizing (points) - 250 pt total width
PANEL_WIDTH_PT = 350;
PANEL_HEIGHT_PT = 200;

% Marker sizes (in points)
CELL_MARKER_SIZE_PT = 3.5;    % Individual cells
ANIMAL_MARKER_SIZE_PT = 7;    % Animal means

% Connecting line width (animal mean to cells)
CONNECTING_LINE_WIDTH = 0.5;
CONNECTING_LINE_ALPHA = 0.25;  % Transparency for connecting lines

% Font settings
FONT_NAME = 'Arial';
FONT_SIZE = 8;

% Line settings
LINE_WIDTH = 1;  % 1 pt lines for axes and box plots

%% ===================== LOAD DATA =====================

fprintf('Loading data from: %s\n', CSV_FILE);
data = readtable(CSV_FILE);

% Extract groups
grp_A = data(strcmp(data.group, 'A'), :);
grp_B = data(strcmp(data.group, 'B'), :);

fprintf('  Group A: %d imaging fields from %d animals\n', ...
    height(grp_A), numel(unique(grp_A.replicate)));
fprintf('  Group B: %d imaging fields from %d animals\n', ...
    height(grp_B), numel(unique(grp_B.replicate)));

% Create output directory
if ~isfolder(OUTPUT_DIR)
    mkdir(OUTPUT_DIR);
    fprintf('Created output directory: %s\n', OUTPUT_DIR);
end

%% ===================== DEFINE COLOR SCHEMES =====================

% Get unique animals per group
animals_A = unique(grp_A.replicate, 'stable');
animals_B = unique(grp_B.replicate, 'stable');
n_animals_A = numel(animals_A);
n_animals_B = numel(animals_B);

if strcmp(COLOR_SCHEME, 'unique')
    % 14 distinct, colorblind-friendly colors
    all_colors = [
        0.230, 0.299, 0.754;  % Blue
        0.706, 0.016, 0.150;  % Red
        0.336, 0.706, 0.184;  % Green
        0.800, 0.475, 0.655;  % Pink
        1.000, 0.501, 0.000;  % Orange
        0.400, 0.200, 0.600;  % Purple
        0.650, 0.650, 0.650;  % Gray
        0.900, 0.745, 0.000;  % Gold
        0.200, 0.630, 0.792;  % Cyan
        0.940, 0.340, 0.360;  % Coral
        0.416, 0.540, 0.200;  % Olive
        0.600, 0.400, 0.200;  % Brown
        0.690, 0.612, 0.851;  % Lavender
        0.200, 0.400, 0.400;  % Teal
    ];
    colors_A = all_colors(1:n_animals_A, :);
    colors_B = all_colors(n_animals_A+1:n_animals_A+n_animals_B, :);
    
else  % 'family' scheme
    % Group A: Shades of gray (dark to medium)
    gray_base = linspace(0.25, 0.65, n_animals_A)';
    colors_A = [gray_base, gray_base, gray_base];
    
    % Group B: Shades of red/coral
    colors_B = [
        0.850, 0.200, 0.200;  % Deep red
        0.950, 0.350, 0.300;  % Coral red
        0.800, 0.100, 0.300;  % Crimson
        0.900, 0.450, 0.350;  % Salmon
        0.700, 0.150, 0.150;  % Dark red
        0.950, 0.250, 0.450;  % Rose
        0.750, 0.300, 0.250;  % Terra cotta
    ];
    colors_B = colors_B(1:n_animals_B, :);
end

%% ===================== GENERATE FIGURES =====================

% Convert points to inches for MATLAB figure sizing (72 pt = 1 inch)
fig_width_in = PANEL_WIDTH_PT / 72;
fig_height_in = PANEL_HEIGHT_PT / 72;

% Convert marker sizes from points to MATLAB scatter units
cell_marker_area = CELL_MARKER_SIZE_PT^2 * (pi/4);
animal_marker_area = ANIMAL_MARKER_SIZE_PT^2 * (pi/4);

for m = 1:numel(METRICS)
    metric = METRICS{m};
    fprintf('\nGenerating SuperPlot for: %s\n', metric);
    
    % Extract data for this metric
    vals_A = grp_A.(metric);
    vals_B = grp_B.(metric);
    reps_A = grp_A.replicate;
    reps_B = grp_B.replicate;
    
    % Calculate animal means
    animal_means_A = zeros(n_animals_A, 1);
    animal_cells_A = cell(n_animals_A, 1);  % Store indices for each animal
    for i = 1:n_animals_A
        idx = find(strcmp(reps_A, animals_A{i}));
        animal_cells_A{i} = idx;
        animal_means_A(i) = mean(vals_A(idx));
    end
    
    animal_means_B = zeros(n_animals_B, 1);
    animal_cells_B = cell(n_animals_B, 1);
    for i = 1:n_animals_B
        idx = find(strcmp(reps_B, animals_B{i}));
        animal_cells_B{i} = idx;
        animal_means_B(i) = mean(vals_B(idx));
    end
    
    % Group statistics (of animal means)
    mean_A = mean(animal_means_A);
    mean_B = mean(animal_means_B);
    sd_A = std(animal_means_A);
    sd_B = std(animal_means_B);
    
    % ===== PRE-COMPUTE ALL POSITIONS =====
    
    % X positions for the three display types
    pos_super = [1.0, 3.0];
    pos_box = [5.5, 6.5];
    pos_meansd = [8.0, 9.0];
    
    % Offset for animal means (to the right of cells)
    animal_mean_offset = 1.0;
    
    % Wide feathering for cell swarm
    swarm_width = 0.6;
    
    % Spread for animal means (vertical jitter width)
    animal_mean_spread = 0.9;
    
    % Get y-range for scaling the feathering algorithm
    all_vals = [vals_A; vals_B];
    y_min = min(all_vals);
    y_max = max(all_vals);
    y_range = y_max - y_min;
    if y_range == 0, y_range = 1; end
    
    % Compute swarm positions for ALL cells at once (ensures good global spread)
    cell_x_A = feathered_beeswarm(vals_A, pos_super(1), swarm_width, y_range);
    cell_x_B = feathered_beeswarm(vals_B, pos_super(2), swarm_width, y_range);
    
    % Compute animal mean positions (spread out to avoid overlap)
    mean_x_A = pos_super(1) + animal_mean_offset + spread_positions(animal_means_A, animal_mean_spread, y_range);
    mean_x_B = pos_super(2) + animal_mean_offset + spread_positions(animal_means_B, animal_mean_spread, y_range);
    
    % Create figure
    fig = figure('Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in], ...
        'PaperUnits', 'inches', 'PaperSize', [fig_width_in fig_height_in], ...
        'PaperPosition', [0 0 fig_width_in fig_height_in], ...
        'Color', 'w', 'Visible', 'off');
    
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % ===== SUPERPLOT: DRAW CONNECTING LINES FIRST (behind everything) =====
    
    % Group A connecting lines
    for i = 1:n_animals_A
        col = colors_A(i, :);
        cell_indices = animal_cells_A{i};
        animal_mean_y = animal_means_A(i);
        animal_mean_x = mean_x_A(i);
        
        for j = 1:numel(cell_indices)
            idx = cell_indices(j);
            plot(ax, [cell_x_A(idx), animal_mean_x], [vals_A(idx), animal_mean_y], ...
                '-', 'Color', [col, CONNECTING_LINE_ALPHA], 'LineWidth', CONNECTING_LINE_WIDTH);
        end
    end
    
    % Group B connecting lines
    for i = 1:n_animals_B
        col = colors_B(i, :);
        cell_indices = animal_cells_B{i};
        animal_mean_y = animal_means_B(i);
        animal_mean_x = mean_x_B(i);
        
        for j = 1:numel(cell_indices)
            idx = cell_indices(j);
            plot(ax, [cell_x_B(idx), animal_mean_x], [vals_B(idx), animal_mean_y], ...
                '-', 'Color', [col, CONNECTING_LINE_ALPHA], 'LineWidth', CONNECTING_LINE_WIDTH);
        end
    end
    
    % ===== SUPERPLOT: DRAW CELL POINTS =====
    
    % Group A cells
    for i = 1:n_animals_A
        col = colors_A(i, :);
        cell_indices = animal_cells_A{i};
        scatter(ax, cell_x_A(cell_indices), vals_A(cell_indices), cell_marker_area, col, 'filled', ...
            'MarkerEdgeColor', 'none');
    end
    
    % Group B cells
    for i = 1:n_animals_B
        col = colors_B(i, :);
        cell_indices = animal_cells_B{i};
        scatter(ax, cell_x_B(cell_indices), vals_B(cell_indices), cell_marker_area, col, 'filled', ...
            'MarkerEdgeColor', 'none');
    end
    
    % ===== SUPERPLOT: DRAW ANIMAL MEANS =====
    
    % Group A animal means
    for i = 1:n_animals_A
        col = colors_A(i, :);
        scatter(ax, mean_x_A(i), animal_means_A(i), animal_marker_area, col, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', LINE_WIDTH);
    end
    
    % Group B animal means
    for i = 1:n_animals_B
        col = colors_B(i, :);
        scatter(ax, mean_x_B(i), animal_means_B(i), animal_marker_area, col, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', LINE_WIDTH);
    end
    
    % ===== BOX PLOT (middle) =====
    
    draw_boxplot(ax, animal_means_A, pos_box(1), 0.5, [0.3 0.3 0.3], LINE_WIDTH);
    draw_boxplot(ax, animal_means_B, pos_box(2), 0.5, [0.7 0.1 0.1], LINE_WIDTH);
    
    % ===== MEAN ± SD (right) =====
    
    errorbar(ax, pos_meansd(1), mean_A, sd_A, 'o', ...
        'Color', [0 0 0], 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', ANIMAL_MARKER_SIZE_PT * 0.8, ...
        'LineWidth', LINE_WIDTH, 'CapSize', 8);
    errorbar(ax, pos_meansd(2), mean_B, sd_B, 'o', ...
        'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.1 0.1], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', ANIMAL_MARKER_SIZE_PT * 0.8, ...
        'LineWidth', LINE_WIDTH, 'CapSize', 8);
    
    % ===== FORMATTING =====
    
    ylim(ax, [y_min - 0.08*y_range, y_max + 0.12*y_range]);
    xlim(ax, [0, 10.0]);
    
    % X-axis labels
    super_center = mean([pos_super(1), pos_super(2) + animal_mean_offset]);
    box_center = mean(pos_box);
    meansd_center = mean(pos_meansd);
    
    xticks(ax, [super_center, box_center, meansd_center]);
    xticklabels(ax, {'SuperPlot', 'Box (N=7)', 'Mean±SD'});
    
    % Add group labels below each pair
    y_label_pos = y_min - 0.15*y_range;
    
    text(ax, pos_super(1), y_label_pos, 'A', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    text(ax, pos_super(2), y_label_pos, 'B', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    text(ax, pos_box(1), y_label_pos, 'A', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    text(ax, pos_box(2), y_label_pos, 'B', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    text(ax, pos_meansd(1), y_label_pos, 'A', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    text(ax, pos_meansd(2), y_label_pos, 'B', ...
        'HorizontalAlignment', 'center', 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    
    ylabel(ax, YLABELS{m}, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE);
    title(ax, TITLES{m}, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, 'FontWeight', 'bold');
    
    set(ax, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE, ...
        'Box', 'off', 'TickDir', 'out', 'LineWidth', LINE_WIDTH, ...
        'XColor', 'k', 'YColor', 'k');
    
    % Vertical separator lines
    y_lims = ylim(ax);
    sep1 = pos_super(2) + animal_mean_offset + 0.4;
    sep2 = pos_box(2) + 0.5;
    plot(ax, [sep1 sep1], y_lims, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    plot(ax, [sep2 sep2], y_lims, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    
    hold(ax, 'off');
    
    % ===== SAVE FIGURE =====
    
    base_name = sprintf('%s_superplot_%s', metric, COLOR_SCHEME);
    
    print(fig, fullfile(OUTPUT_DIR, [base_name '.png']), '-dpng', '-r300');
    print(fig, fullfile(OUTPUT_DIR, [base_name '.eps']), '-depsc', '-painters');
    print(fig, fullfile(OUTPUT_DIR, [base_name '.pdf']), '-dpdf', '-painters');
    
    fprintf('  Saved: %s.eps, .pdf, .png\n', base_name);
    
    close(fig);
end

fprintf('\n✓ All figures saved to: %s\n', OUTPUT_DIR);
fprintf('\nTo switch color schemes, change COLOR_SCHEME to ''unique'' or ''family'' and re-run.\n');

%% ===================== HELPER FUNCTIONS =====================

function x_positions = feathered_beeswarm(values, center, width, y_range)
    % Generate well-feathered beeswarm positions using force-directed simulation
    % Ensures minimum spacing between nearby points
    
    n = numel(values);
    x_positions = zeros(n, 1);
    
    if n == 0
        return;
    end
    
    if n == 1
        x_positions = center;
        return;
    end
    
    % Normalize y values to [0, 1] for distance calculations
    y_norm = (values - min(values)) / y_range;
    
    % Target minimum distance between points (in normalized coordinates)
    % This determines how aggressively points are spread
    min_dist = 0.025;  % Smaller = more spread for close points
    
    % Initialize x positions at center
    x_pos = zeros(n, 1);
    
    % Sort by y value for sequential placement
    [~, sort_idx] = sort(values);
    
    % Place points sequentially, pushing away from nearby points
    for i = 1:n
        idx = sort_idx(i);
        y_val = y_norm(idx);
        
        if i == 1
            x_pos(idx) = 0;
            continue;
        end
        
        % Find all previously placed points
        placed_idx = sort_idx(1:i-1);
        
        % Calculate distances in y to all placed points
        y_dists = abs(y_norm(placed_idx) - y_val);
        
        % Find points that are "close" in y (potential overlaps)
        close_mask = y_dists < min_dist * 4;
        
        if ~any(close_mask)
            % No nearby points, place at center
            x_pos(idx) = 0;
        else
            % Find best x position to avoid overlaps
            close_points_x = x_pos(placed_idx(close_mask));
            close_points_y_dist = y_dists(close_mask);
            
            % Weight by y-distance (closer in y = need more x separation)
            weights = 1 - close_points_y_dist / (min_dist * 4);
            
            % Try positions across the width
            test_x = linspace(-width/2, width/2, 101);
            best_x = 0;
            best_score = -inf;
            
            for test_pos = test_x
                % Calculate minimum weighted distance to close points
                x_dists = abs(close_points_x - test_pos);
                
                % Combined distance (x and y)
                combined_dists = sqrt(x_dists.^2 + close_points_y_dist.^2);
                
                % Score: minimum distance to any close point, penalize edges
                min_combined = min(combined_dists);
                edge_penalty = abs(test_pos) * 0.1;  % Slight preference for center
                score = min_combined - edge_penalty;
                
                if score > best_score
                    best_score = score;
                    best_x = test_pos;
                end
            end
            
            x_pos(idx) = best_x;
        end
    end
    
    % Apply force-directed refinement (multiple iterations)
    n_iterations = 50;
    
    for iter = 1:n_iterations
        forces = zeros(n, 1);
        
        for i = 1:n
            for j = 1:n
                if i == j, continue; end
                
                % Distance in y (normalized)
                dy = abs(y_norm(i) - y_norm(j));
                
                % Only consider points close in y
                if dy > min_dist * 5
                    continue;
                end
                
                % Distance in x
                dx = x_pos(i) - x_pos(j);
                
                % Combined distance
                dist = sqrt(dx^2 + dy^2);
                
                if dist < min_dist && dist > 0
                    % Repulsive force (stronger when closer)
                    force_mag = (min_dist - dist) / min_dist;
                    
                    % Direction: push away in x
                    if dx == 0
                        dx = (rand - 0.5) * 0.01;  % Random nudge if exactly aligned
                    end
                    forces(i) = forces(i) + force_mag * sign(dx) * 0.02;
                end
            end
            
            % Centering force (weak pull toward center)
            forces(i) = forces(i) - x_pos(i) * 0.01;
        end
        
        % Apply forces
        x_pos = x_pos + forces;
        
        % Clamp to width
        x_pos = max(-width/2, min(width/2, x_pos));
    end
    
    % Return positions centered on 'center'
    x_positions = center + x_pos;
end

function x_spread = spread_positions(values, width, y_range)
    % Spread points horizontally using force-directed approach
    % All points are spread out, with more separation for points close in y
    
    n = numel(values);
    x_spread = zeros(n, 1);
    
    if n <= 1
        return;
    end
    
    if n == 2
        % Simple case: spread two points apart
        x_spread = [-width/3; width/3];
        [~, sort_idx] = sort(values);
        x_spread(sort_idx) = x_spread;
        return;
    end
    
    % Normalize values for distance calculations
    y_norm = (values - min(values));
    if max(y_norm) > 0
        y_norm = y_norm / max(y_norm);
    end
    
    % Initialize with slight spread based on rank
    [~, sort_idx] = sort(values);
    init_spread = linspace(-width/2, width/2, n);
    x_spread(sort_idx) = init_spread';
    
    % Force-directed refinement
    min_dist = 0.15;  % Minimum desired distance between points
    n_iterations = 100;
    
    for iter = 1:n_iterations
        forces = zeros(n, 1);
        
        for i = 1:n
            for j = 1:n
                if i == j, continue; end
                
                % Distance in y (normalized)
                dy = abs(y_norm(i) - y_norm(j));
                
                % Distance in x
                dx = x_spread(i) - x_spread(j);
                
                % Combined distance
                dist = sqrt(dx^2 + dy^2);
                
                if dist < min_dist && dist > 0
                    % Repulsive force
                    force_mag = (min_dist - dist) / min_dist;
                    
                    if abs(dx) < 0.001
                        dx = (rand - 0.5) * 0.01;
                    end
                    forces(i) = forces(i) + force_mag * sign(dx) * 0.05;
                end
            end
            
            % Weak centering force
            forces(i) = forces(i) - x_spread(i) * 0.005;
        end
        
        % Apply forces
        x_spread = x_spread + forces;
        
        % Clamp to width
        x_spread = max(-width/2, min(width/2, x_spread));
    end
end

function draw_boxplot(ax, data, x_center, width, color, line_width)
    % Draw a custom box-and-whisker plot
    
    q1 = prctile(data, 25);
    q2 = median(data);
    q3 = prctile(data, 75);
    iqr_val = q3 - q1;
    
    whisker_low = max(min(data), q1 - 1.5*iqr_val);
    whisker_high = min(max(data), q3 + 1.5*iqr_val);
    
    half_width = width / 2;
    
    % Box
    rectangle(ax, 'Position', [x_center - half_width, q1, width, q3 - q1], ...
        'EdgeColor', color, 'LineWidth', line_width, 'FaceColor', 'none');
    
    % Median line
    plot(ax, [x_center - half_width, x_center + half_width], [q2, q2], ...
        '-', 'Color', color, 'LineWidth', line_width * 1.5);
    
    % Whiskers
    plot(ax, [x_center, x_center], [q1, whisker_low], '-', 'Color', color, 'LineWidth', line_width);
    plot(ax, [x_center, x_center], [q3, whisker_high], '-', 'Color', color, 'LineWidth', line_width);
    
    % Caps
    cap_width = width * 0.4;
    plot(ax, [x_center - cap_width/2, x_center + cap_width/2], [whisker_low, whisker_low], ...
        '-', 'Color', color, 'LineWidth', line_width);
    plot(ax, [x_center - cap_width/2, x_center + cap_width/2], [whisker_high, whisker_high], ...
        '-', 'Color', color, 'LineWidth', line_width);
end
