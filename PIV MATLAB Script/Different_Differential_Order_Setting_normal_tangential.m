% Analyze and plot using saved preprocessed data
% Each file's Velocity Field with Contours plot is displayed separately and saved as .fig format
% Use grayscale images as background for all 2D plots except velocity field
% Allow epsilon to be tuned as a parameter

clear; clc; close all;

%% Load preprocessed data
fprintf('Loading preprocessed data...\n');
if ~exist('preprocessed_data.mat', 'file')
    error('preprocessed_data.mat not found. Please run the main analysis script first.');
end

load('preprocessed_data.mat', 'preprocessed_data');

num_files = length(preprocessed_data);
fprintf('Loaded preprocessed data for %d files\n', num_files);

%% User-defined input for new parameters
M = input('Enter the bin number M: ');
diff_order = input(['Enter the differential order n for' ...
    ' central difference (e.g., 1, 2, 3): ']);
epsilon = input('Enter the epsilon value: ');

fprintf('Using M = %d for reanalysis\n', M);
fprintf('Using differential order n = %d for central difference\n', diff_order);
fprintf('Using epsilon = %.1f μm for reanalysis\n', epsilon);

%% Create main output folder
main_output_folder = sprintf('Reanalysis_Plots_M%d_n%d_eps%.1f', M, diff_order, epsilon);
if ~exist(main_output_folder, 'dir')
    mkdir(main_output_folder);
    fprintf('Created main output folder: %s\n', main_output_folder);
end

%% Create categorized subfolders for plots without contours
subfolders = {
    'VelocityField_NoContours', 
    'Divergence_NoContours', 
    'H_rr_NoContours',
    'H_theta_theta_NoContours'
};

for i = 1:length(subfolders)
    folder_path = fullfile(main_output_folder, subfolders{i});
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
        fprintf('Created subfolder: %s\n', folder_path);
    end
end

%% Create subfolders for plots with grayscale background (excluding velocity field)
grayscale_subfolders = {
    'Divergence_WithGrayBackground', 
    'H_rr_WithGrayBackground',
    'H_theta_theta_WithGrayBackground'
};

for i = 1:length(grayscale_subfolders)
    folder_path = fullfile(main_output_folder, grayscale_subfolders{i});
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
        fprintf('Created subfolder: %s\n', folder_path);
    end
end

%% Get tif files for grayscale background
tifFiles = dir('*.tif');
if isempty(tifFiles)
    tifFiles = dir('*.tiff');
end

% Sort files in natural order
tifFiles = sort_natural(tifFiles);

if length(tifFiles) ~= num_files
    fprintf('Warning: Number of tif files (%d) does not match number of data files (%d)\n', ...
        length(tifFiles), num_files);
    fprintf('Will attempt to match files by order\n');
end

% Conversion factor for unit alignment
pixels_per_um = 0.8532;
um_per_m = 1e6;
pixels_per_m = pixels_per_um * um_per_m;

fprintf('Using conversion factor: %.4f pixels/micrometer\n', pixels_per_um);

%% Compute regions and mean values with new M and epsilon
fprintf('\nRecomputing regions and mean values with new M = %d and epsilon = %.1f...\n', M, epsilon);

% Initialize arrays to store recomputed results
all_rMean_new = cell(num_files, 1);
all_radial_strain_rate_Mean_new = cell(num_files, 1);
all_tangential_strain_rate_Mean_new = cell(num_files, 1);
all_divergence_Mean_new = cell(num_files, 1);

% Initialize arrays to store gradient data for each file
gradient_data = struct();

% Define colors and labels
colors = lines(num_files);
legendLabels = {'4hr', '8hr', '12hr', '18hr', '24hr'};

for file_idx = 1:num_files
    fprintf('Reanalyzing file %d/%d: %s\n', file_idx, num_files, preprocessed_data(file_idx).filename);
    
    % Extract variables from preprocessed data
    X = preprocessed_data(file_idx).X;
    Y = preprocessed_data(file_idx).Y;
    R = preprocessed_data(file_idx).R;
    U_grid = preprocessed_data(file_idx).U_grid;
    V_grid = preprocessed_data(file_idx).V_grid;
    SDF = preprocessed_data(file_idx).SDF;
    
    % Get unique x and y values
    unique_x = unique(X(1,:));
    unique_y = unique(Y(:,1));

    % Calculate grid spacing
    if length(unique_x) > 1 && length(unique_y) > 1
        dx = unique_x(2) - unique_x(1);
        dy = unique_y(2) - unique_y(1);
    else
        fprintf('Warning: Insufficient grid points, unable to calculate grid spacing\n');
        continue;
    end
    h_x = dx;
    h_y = dy;

    %% Custom Numerical Differentiation with First Contour SDF-aware boundary handling
    fprintf('  Calculating gradients with First Contour SDF-aware numerical differentiation (n=%d, epsilon=%.1f)...\n', diff_order, epsilon);

    % Initialize gradient matrices
    [ny, nx] = size(U_grid);
    du_dx = zeros(ny, nx);
    du_dy = zeros(ny, nx);
    dv_dx = zeros(ny, nx);
    dv_dy = zeros(ny, nx);

    % Calculate gradients with First Contour SDF-aware method using new epsilon
    for i = 1:ny
        for j = 1:nx
            % Current point SDF value
            current_sdf = SDF(i,j);
            
            %% Calculate du/dx
            if current_sdf < epsilon  % Use new epsilon
                % Point inside first contour - set gradient to NaN
                du_dx(i,j) = NaN;
            else
                % Point outside first contour - check neighbors
                left_j = j - diff_order;
                right_j = j + diff_order;
                
                % Check if neighbors are within grid boundaries AND in SDF>=epsilon_new region
                left_valid = (left_j >= 1) && (left_j <= nx) && (SDF(i, left_j) >= epsilon);
                right_valid = (right_j >= 1) && (right_j <= nx) && (SDF(i, right_j) >= epsilon);
                
                if left_valid && right_valid
                    % Both neighbors valid - use central difference
                    du_dx(i,j) = (U_grid(i, right_j) - U_grid(i, left_j)) / (2 * diff_order * h_x);
                elseif right_valid
                    % Only right neighbor valid - use forward difference
                    du_dx(i,j) = (U_grid(i, right_j) - U_grid(i, j)) / (diff_order * h_x);
                elseif left_valid
                    % Only left neighbor valid - use backward difference
                    du_dx(i,j) = (U_grid(i, j) - U_grid(i, left_j)) / (diff_order * h_x);
                else
                    % No valid neighbors - cannot compute gradient
                    du_dx(i,j) = NaN;
                end
            end
            
            %% Calculate du/dy
            if current_sdf < epsilon  % Use new epsilon
                % Point inside first contour - set gradient to NaN
                du_dy(i,j) = NaN;
            else
                % Point outside first contour - check neighbors
                bottom_i = i - diff_order;
                top_i = i + diff_order;
                
                % Check if neighbors are within grid boundaries AND in SDF>=epsilon_new region
                bottom_valid = (bottom_i >= 1) && (bottom_i <= ny) && (SDF(bottom_i, j) >= epsilon);
                top_valid = (top_i >= 1) && (top_i <= ny) && (SDF(top_i, j) >= epsilon);
                
                if bottom_valid && top_valid
                    % Both neighbors valid - use central difference
                    du_dy(i,j) = (U_grid(top_i, j) - U_grid(bottom_i, j)) / (2 * diff_order * h_y);
                elseif top_valid
                    % Only top neighbor valid - use forward difference
                    du_dy(i,j) = (U_grid(top_i, j) - U_grid(i, j)) / (diff_order * h_y);
                elseif bottom_valid
                    % Only bottom neighbor valid - use backward difference
                    du_dy(i,j) = (U_grid(i, j) - U_grid(bottom_i, j)) / (diff_order * h_y);
                else
                    % No valid neighbors - cannot compute gradient
                    du_dy(i,j) = NaN;
                end
            end
            
            %% Calculate dv/dx
            if current_sdf < epsilon  % Use new epsilon
                % Point inside first contour - set gradient to NaN
                dv_dx(i,j) = NaN;
            else
                % Point outside first contour - check neighbors
                left_j = j - diff_order;
                right_j = j + diff_order;
                
                % Check if neighbors are within grid boundaries AND in SDF>=epsilon_new region
                left_valid = (left_j >= 1) && (left_j <= nx) && (SDF(i, left_j) >= epsilon);
                right_valid = (right_j >= 1) && (right_j <= nx) && (SDF(i, right_j) >= epsilon);
                
                if left_valid && right_valid
                    % Both neighbors valid - use central difference
                    dv_dx(i,j) = (V_grid(i, right_j) - V_grid(i, left_j)) / (2 * diff_order * h_x);
                elseif right_valid
                    % Only right neighbor valid - use forward difference
                    dv_dx(i,j) = (V_grid(i, right_j) - V_grid(i, j)) / (diff_order * h_x);
                elseif left_valid
                    % Only left neighbor valid - use backward difference
                    dv_dx(i,j) = (V_grid(i, j) - V_grid(i, left_j)) / (diff_order * h_x);
                else
                    % No valid neighbors - cannot compute gradient
                    dv_dx(i,j) = NaN;
                end
            end
            
            %% Calculate dv/dy
            if current_sdf < epsilon  % Use new epsilon
                % Point inside first contour - set gradient to NaN
                dv_dy(i,j) = NaN;
            else
                % Point outside first contour - check neighbors
                bottom_i = i - diff_order;
                top_i = i + diff_order;
                
                % Check if neighbors are within grid boundaries AND in SDF>=epsilon_new region
                bottom_valid = (bottom_i >= 1) && (bottom_i <= ny) && (SDF(bottom_i, j) >= epsilon);
                top_valid = (top_i >= 1) && (top_i <= ny) && (SDF(top_i, j) >= epsilon);
                
                if bottom_valid && top_valid
                    % Both neighbors valid - use central difference
                    dv_dy(i,j) = (V_grid(top_i, j) - V_grid(bottom_i, j)) / (2 * diff_order * h_y);
                elseif top_valid
                    % Only top neighbor valid - use forward difference
                    dv_dy(i,j) = (V_grid(top_i, j) - V_grid(i, j)) / (diff_order * h_y);
                elseif bottom_valid
                    % Only bottom neighbor valid - use backward difference
                    dv_dy(i,j) = (V_grid(i, j) - V_grid(bottom_i, j)) / (diff_order * h_y);
                else
                    % No valid neighbors - cannot compute gradient
                    dv_dy(i,j) = NaN;
                end
            end
        end
    end

    % Store gradient data for this file
    gradient_data(file_idx).filename = preprocessed_data(file_idx).filename;
    gradient_data(file_idx).du_dx = du_dx;
    gradient_data(file_idx).du_dy = du_dy;
    gradient_data(file_idx).dv_dx = dv_dx;
    gradient_data(file_idx).dv_dy = dv_dy;
    gradient_data(file_idx).diff_order = diff_order;
    gradient_data(file_idx).M = M;
    gradient_data(file_idx).epsilon = epsilon;  % Store epsilon

    % Calculate angle theta for each point
    [dSDF_dx, dSDF_dy] = gradient(SDF);
    theta = atan2(dSDF_dy, dSDF_dx);  % atan2 gives angle in radians from -pi to pi
    
    % Calculate cos_theta and sin_theta
    cos_theta = cos(theta);
    sin_theta = sin(theta);

    div_velocity = du_dx + dv_dy;
    H_rr = du_dx.* (cos_theta).^2 + (du_dy + dv_dx).* cos_theta.* sin_theta + dv_dy.* (sin_theta).^2;
    H_theta_theta = du_dx.* (sin_theta).^2 - (du_dy + dv_dx).* cos_theta.* sin_theta + dv_dy.* (cos_theta).^2;
    tr_H = H_rr + H_theta_theta;
    
    % Generate new region boundaries starting from new epsilon
    region_boundaries_new = zeros(1, M+1);
    for k = 0:M
        region_boundaries_new(k+1) = epsilon + k * 100;  % Use new epsilon
    end
    
    % Initialize arrays to store new results
    rMean_new = zeros(M, 1);
    radial_strain_rate_Mean_new = zeros(M, 1);
    tangential_strain_rate_Mean_new = zeros(M, 1);
    divergence_Mean_new = zeros(M, 1);
    
    % Store masks for each region
    region_masks = cell(M, 1);
    
    % Calculate mean values for each new region using new epsilon
    for i = 1:M
        inner_distance = region_boundaries_new(i);
        outer_distance = region_boundaries_new(i+1);
        
        % Create mask for region between these two contours
        region_mask = (SDF >= inner_distance) & (SDF < outer_distance) & (SDF >= 0);
        region_masks{i} = region_mask;
        
        % Calculate mean values (ignoring NaN values)
        if sum(region_mask(:)) > 0
            valid_mask = region_mask & ~isnan(H_rr) & ~isnan(H_theta_theta) & ~isnan(div_velocity);
            if sum(valid_mask(:)) > 0
                rMean_new(i) = mean(R(valid_mask));
                radial_strain_rate_Mean_new(i) = mean(H_rr(valid_mask));
                tangential_strain_rate_Mean_new(i) = mean(H_theta_theta(valid_mask));
                divergence_Mean_new(i) = mean(div_velocity(valid_mask));
            else
                rMean_new(i) = NaN;
                radial_strain_rate_Mean_new(i) = NaN;
                tangential_strain_rate_Mean_new(i) = NaN;
                divergence_Mean_new(i) = NaN;
            end
        else
            rMean_new(i) = NaN;
            radial_strain_rate_Mean_new(i) = NaN;
            tangential_strain_rate_Mean_new(i) = NaN;
            divergence_Mean_new(i) = NaN;
        end
    end
    
    % Store results
    all_rMean_new{file_idx} = rMean_new;
    all_radial_strain_rate_Mean_new{file_idx} = radial_strain_rate_Mean_new;
    all_tangential_strain_rate_Mean_new{file_idx} = tangential_strain_rate_Mean_new;
    all_divergence_Mean_new{file_idx} = divergence_Mean_new;
    
    %% Calculate radial velocity
    v_r = U_grid .* cos_theta + V_grid .* sin_theta;
    
    % Create custom colormap for divergence and gradients
    N = 256; half = round(N/2);
    blue = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    red  = [ones(half,1),        linspace(1,0,half)', linspace(1,0,half)'];
    custom_cmap = [blue; red];
    
    % Region colors
    region_colors = jet(M);
    
    %% Load corresponding grayscale image for background
    gray_background = [];
    if file_idx <= length(tifFiles)
        try
            tifFilename = tifFiles(file_idx).name;
            fprintf('  Loading grayscale background: %s\n', tifFilename);
            gray_background = imread(tifFilename);
            
            % Convert to RGB if grayscale
            if ndims(gray_background) == 2
                gray_background = repmat(gray_background, [1, 1, 3]);
            end
            
            % Normalize and enhance contrast for better visibility
            gray_background = double(gray_background) / double(max(gray_background(:)));
            
            % Enhance brightness and contrast for better visibility
            gray_background = imadjust(gray_background, [0.1; 0.9], [0; 1]);
            
            % Calculate the physical coordinates that correspond to the background image
            [bg_height, bg_width, ~] = size(gray_background);
            x_min = min(X(:));
            x_max = max(X(:));
            y_min = min(Y(:));
            y_max = max(Y(:));
            
        catch ME
            fprintf('  Warning: Could not load grayscale image %s: %s\n', tifFilename, ME.message);
            gray_background = [];
        end
    else
        fprintf('  Warning: No grayscale image available for file index %d\n', file_idx);
    end
    
    %% Create individual Velocity Field with Contours plot for each file
    fprintf('  Creating velocity field with contours plot for file %d...\n', file_idx);
    
    fig_contours = figure('Name', sprintf('Velocity Field with Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                         preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                         'NumberTitle', 'off', 'Position', [150, 150, 1000, 800]);
    
    % Calculate velocity magnitude
    velocity_magnitude = sqrt(U_grid.^2 + V_grid.^2);
    
    % Plot velocity field
    h = pcolor(X, Y, velocity_magnitude);
    set(h, 'EdgeColor', 'none');
    shading interp; 
    colorbar;
    hold on;
    
    % Store contour handles and labels for legend
    contour_handles = [];
    contour_labels = {};
    
    % Generate and plot all contours using epsilon
    for i = 1:length(region_boundaries_new)
        distance = region_boundaries_new(i);
        
        % Extract contour at this SDF level
        C = contourc(unique_x, unique_y, SDF, [distance, distance]);
        
        % Process contour data
        [contour_x, contour_y] = processContourData(C);
        
        % Plot contour
        if ~isempty(contour_x)
            if abs(distance - epsilon) < 1e-6  % Compare with epsilon
                % First contour - use special color
                h_contour = plot(contour_x, contour_y, 'm-', 'LineWidth', 3);
                contour_handles = [contour_handles, h_contour];
                contour_labels{end+1} = sprintf('First Contour (%.1f μm)', distance);
            else
                % Other contours
                color_idx = mod(i-1, M) + 1;
                h_contour = plot(contour_x, contour_y, 'Color', region_colors(color_idx,:), ...
                         'LineWidth', 2);
                contour_handles = [contour_handles, h_contour];
                contour_labels{end+1} = sprintf('Contour: %.1f μm', distance);
            end
        end
    end
    
    % Plot boundary (SDF=0) for reference
    C_boundary = contourc(unique_x, unique_y, SDF, [0, 0]);
    [boundary_contour_x, boundary_contour_y] = processContourData(C_boundary);
    if ~isempty(boundary_contour_x)
        h_boundary = plot(boundary_contour_x, boundary_contour_y, 'r-', 'LineWidth', 3);
        contour_handles = [contour_handles, h_boundary];
        contour_labels{end+1} = 'Wound Boundary (0 μm)';
    end
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    
    % Only show contours in legend
    if ~isempty(contour_handles)
        legend(contour_handles, contour_labels, 'Location', 'eastoutside');
    end
    axis equal;
    
    %% Save velocity field with contours image as .fig format only
    fprintf('  Saving velocity field with contours plot for file %d...\n', file_idx);
    
    % Generate filename
    [~, name_only, ~] = fileparts(preprocessed_data(file_idx).filename);
    filename_fig = sprintf('%s_VelocityField_WithContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    
    % Full path
    fullpath_fig = fullfile(main_output_folder, filename_fig);
    
    % Save as FIG file only
    saveas(fig_contours, fullpath_fig, 'fig');
    
    fprintf('    Saved: %s\n', filename_fig);
    
    %% Create individual Velocity Field WITHOUT Contours plot for each file
    fprintf('  Creating velocity field WITHOUT contours plot for file %d...\n', file_idx);
    
    fig_no_contours = figure('Name', sprintf('Velocity Field WITHOUT Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                         preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                         'NumberTitle', 'off', 'Position', [150, 150, 1000, 800]);
    
    % Calculate velocity magnitude
    velocity_magnitude = sqrt(U_grid.^2 + V_grid.^2);
    
    % Plot velocity field
    h = pcolor(X, Y, velocity_magnitude);
    set(h, 'EdgeColor', 'none');
    shading interp; 
    colorbar;
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('Velocity Field','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save velocity field WITHOUT contours image in multiple formats to categorized folder
    fprintf('  Saving velocity field WITHOUT contours plot for file %d...\n', file_idx);
    
    % Generate filenames for multiple formats
    filename_fig_no_contours = sprintf('%s_VelocityField_NoContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    filename_jpeg_no_contours = sprintf('%s_VelocityField_NoContours_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
    filename_eps_no_contours = sprintf('%s_VelocityField_NoContours_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
    
    % Full paths to categorized folder
    velocity_folder = fullfile(main_output_folder, 'VelocityField_NoContours');
    fullpath_fig_no_contours = fullfile(velocity_folder, filename_fig_no_contours);
    fullpath_jpeg_no_contours = fullfile(velocity_folder, filename_jpeg_no_contours);
    fullpath_eps_no_contours = fullfile(velocity_folder, filename_eps_no_contours);
    
    % Save in multiple formats
    saveas(fig_no_contours, fullpath_fig_no_contours, 'fig');
    saveas(fig_no_contours, fullpath_jpeg_no_contours, 'jpeg');
    saveas(fig_no_contours, fullpath_eps_no_contours, 'eps');
    
    fprintf('    Saved: %s\n', filename_fig_no_contours);
    fprintf('    Saved: %s\n', filename_jpeg_no_contours);
    fprintf('    Saved: %s\n', filename_eps_no_contours);
    
    %% Create individual Divergence 2D plot with contours for each file
    fprintf('  Creating divergence 2D plot with contours for file %d...\n', file_idx);
    
    fig_divergence_contours = figure('Name', sprintf('Divergence Field with Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [350, 150, 1000, 800]);
    
    % Plot divergence field
    h = pcolor(X, Y, div_velocity);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    hold on;
    
    % Store contour handles and labels for legend
    div_contour_handles = [];
    div_contour_labels = {};
    
    % Generate and plot all contours using epsilon
    for i = 1:length(region_boundaries_new)
        distance = region_boundaries_new(i);
        
        % Extract contour at this SDF level
        C = contourc(unique_x, unique_y, SDF, [distance, distance]);
        
        % Process contour data
        [contour_x, contour_y] = processContourData(C);
        
        % Plot contour
        if ~isempty(contour_x)
            if abs(distance - epsilon) < 1e-6  % Compare with new epsilon
                % First contour - use special color
                h_contour = plot(contour_x, contour_y, 'm-', 'LineWidth', 3);
                div_contour_handles = [div_contour_handles, h_contour];
                div_contour_labels{end+1} = sprintf('First Contour (%.1f μm)', distance);
            else
                % Other contours
                color_idx = mod(i-1, M) + 1;
                h_contour = plot(contour_x, contour_y, 'Color', region_colors(color_idx,:), ...
                         'LineWidth', 2);
                div_contour_handles = [div_contour_handles, h_contour];
                div_contour_labels{end+1} = sprintf('Contour: %.1f μm', distance);
            end
        end
    end
    
    % % Plot boundary (SDF=0) for reference
    % C_boundary = contourc(unique_x, unique_y, SDF, [0, 0]);
    % [boundary_contour_x, boundary_contour_y] = processContourData(C_boundary);
    % if ~isempty(boundary_contour_x)
    %     h_boundary = plot(boundary_contour_x, boundary_contour_y, 'r-', 'LineWidth', 3);
    %     div_contour_handles = [div_contour_handles, h_boundary];
    %     div_contour_labels{end+1} = 'Wound Boundary (0 μm)';
    % end
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('Divergence','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save divergence with contours image as .fig format only
    fprintf('  Saving divergence with contours plot for file %d...\n', file_idx);
    
    filename_div_contours_fig = sprintf('%s_DivergenceField_WithContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    fullpath_div_contours_fig = fullfile(main_output_folder, filename_div_contours_fig);
    saveas(fig_divergence_contours, fullpath_div_contours_fig, 'fig');
    fprintf('    Saved: %s\n', filename_div_contours_fig);
    
    %% Create individual Divergence 2D plot WITHOUT contours for each file
    fprintf('  Creating divergence 2D plot WITHOUT contours for file %d...\n', file_idx);
    
    fig_divergence_no_contours = figure('Name', sprintf('Divergence Field WITHOUT Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [350, 150, 1000, 800]);
    
    % Plot divergence field
    h = pcolor(X, Y, div_velocity);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('Divergence','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save divergence WITHOUT contours image in multiple formats to categorized folder
    fprintf('  Saving divergence WITHOUT contours plot for file %d...\n', file_idx);
    
    filename_div_no_contours_fig = sprintf('%s_DivergenceField_NoContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    filename_div_no_contours_jpeg = sprintf('%s_DivergenceField_NoContours_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
    filename_div_no_contours_eps = sprintf('%s_DivergenceField_NoContours_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
    
    % Full paths to categorized folder
    divergence_folder = fullfile(main_output_folder, 'Divergence_NoContours');
    fullpath_div_no_contours_fig = fullfile(divergence_folder, filename_div_no_contours_fig);
    fullpath_div_no_contours_jpeg = fullfile(divergence_folder, filename_div_no_contours_jpeg);
    fullpath_div_no_contours_eps = fullfile(divergence_folder, filename_div_no_contours_eps);
    
    saveas(fig_divergence_no_contours, fullpath_div_no_contours_fig, 'fig');
    saveas(fig_divergence_no_contours, fullpath_div_no_contours_jpeg, 'jpeg');
    saveas(fig_divergence_no_contours, fullpath_div_no_contours_eps, 'eps');
    
    fprintf('    Saved: %s\n', filename_div_no_contours_fig);
    fprintf('    Saved: %s\n', filename_div_no_contours_jpeg);
    fprintf('    Saved: %s\n', filename_div_no_contours_eps);
    
    %% Create Divergence WITH Grayscale Background - ENHANCED VERSION
    if ~isempty(gray_background)
        fprintf('  Creating divergence WITH grayscale background for file %d...\n', file_idx);
        
        fig_div_gray = figure('Name', sprintf('Divergence with Gray Background - %s (M=%d, n=%d, eps=%.1f)', ...
                             preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                             'NumberTitle', 'off', 'Position', [350, 150, 1000, 800]);
        
        % Display enhanced background
        imagesc([x_min, x_max], [y_min, y_max], gray_background);
        hold on;
        
        % Overlay divergence field with transparency
        h_div = pcolor(X, Y, div_velocity);
        set(h_div, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        shading interp; 
        colormap(custom_cmap);
        caxis([-0.002, 0.002]);
        colorbar;
        
        % ADD WOUND BOUNDARY AND FIRST CONTOUR TO GRAYSCALE BACKGROUND PLOTS
        % Plot first contour (SDF=epsilon_new) using new epsilon
        C_first = contourc(unique_x, unique_y, SDF, [epsilon, epsilon]);
        [first_contour_x, first_contour_y] = processContourData(C_first);
        if ~isempty(first_contour_x)
            plot(first_contour_x, first_contour_y, 'm-', 'LineWidth', 3, 'DisplayName', 'First Contour');
        end
        
        xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
        ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
        title('Divergence','fontsize',16,Interpreter="latex");
        axis equal;
        
        % Save divergence WITH grayscale background
        filename_div_gray_fig = sprintf('%s_Divergence_WithGrayBackground_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
        filename_div_gray_jpeg = sprintf('%s_Divergence_WithGrayBackground_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
        filename_div_gray_eps = sprintf('%s_Divergence_WithGrayBackground_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
        
        % Full paths to grayscale background folder
        div_gray_folder = fullfile(main_output_folder, 'Divergence_WithGrayBackground');
        fullpath_div_gray_fig = fullfile(div_gray_folder, filename_div_gray_fig);
        fullpath_div_gray_jpeg = fullfile(div_gray_folder, filename_div_gray_jpeg);
        fullpath_div_gray_eps = fullfile(div_gray_folder, filename_div_gray_eps);
        
        saveas(fig_div_gray, fullpath_div_gray_fig, 'fig');
        saveas(fig_div_gray, fullpath_div_gray_jpeg, 'jpeg');
        saveas(fig_div_gray, fullpath_div_gray_eps, 'eps');
        
        fprintf('    Saved with enhanced grayscale background: %s\n', filename_div_gray_jpeg);
        
        close(fig_div_gray);
    end

     %% Create individual H_rr 2D plot with contours for each file
    fprintf('  Creating H_rr 2D plot with contours for file %d...\n', file_idx);
    
    fig_H_rr_contours = figure('Name', sprintf('H_rr Field with Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
    
    % Plot H_rr field
    h = pcolor(X, Y, H_rr);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    hold on;
    
    % Store contour handles and labels for legend
    H_rr_contour_handles = [];
    H_rr_contour_labels = {};
    
    % Generate and plot all contours using new epsilon
    for i = 1:length(region_boundaries_new)
        distance = region_boundaries_new(i);
        
        % Extract contour at this SDF level
        C = contourc(unique_x, unique_y, SDF, [distance, distance]);
        
        % Process contour data
        [contour_x, contour_y] = processContourData(C);
        
        % Plot contour
        if ~isempty(contour_x)
            if abs(distance - epsilon) < 1e-6  % Compare with new epsilon
                % First contour - use special color
                h_contour = plot(contour_x, contour_y, 'm-', 'LineWidth', 3);
                H_rr_contour_handles = [H_rr_contour_handles, h_contour];
                H_rr_contour_labels{end+1} = sprintf('First Contour (%.1f μm)', distance);
            else
                % Other contours
                color_idx = mod(i-1, M) + 1;
                h_contour = plot(contour_x, contour_y, 'Color', region_colors(color_idx,:), ...
                         'LineWidth', 2);
                H_rr_contour_handles = [H_rr_contour_handles, h_contour];
                H_rr_contour_labels{end+1} = sprintf('Contour: %.1f μm', distance);
            end
        end
    end
    
    % % Plot boundary (SDF=0) for reference
    % C_boundary = contourc(unique_x, unique_y, SDF, [0, 0]);
    % [boundary_contour_x, boundary_contour_y] = processContourData(C_boundary);
    % if ~isempty(boundary_contour_x)
    %     h_boundary = plot(boundary_contour_x, boundary_contour_y, 'r-', 'LineWidth', 3);
    %     H_rr_contour_handles = [H_rr_contour_handles, h_boundary];
    %     H_rr_contour_labels{end+1} = 'Wound Boundary (0 μm)';
    % end
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('$\partial v_{r} / \partial r$','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save H_rr with contours image as .fig format only
    fprintf('  Saving H_rr with contours plot for file %d...\n', file_idx);
    
    filename_H_rr_contours_fig = sprintf('%s_H_rr_WithContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    fullpath_H_rr_contours_fig = fullfile(main_output_folder, filename_H_rr_contours_fig);
    saveas(fig_H_rr_contours, fullpath_H_rr_contours_fig, 'fig');
    fprintf('    Saved: %s\n', filename_H_rr_contours_fig);
    
    %% Create individual H_rr 2D plot WITHOUT contours for each file
    fprintf('  Creating H_rr 2D plot WITHOUT contours for file %d...\n', file_idx);
    
    fig_H_rr_no_contours = figure('Name', sprintf('H_rr Field WITHOUT Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
    
    % Plot H_rr field
    h = pcolor(X, Y, H_rr);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('$\partial v_{r} / \partial r$','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save H_rr WITHOUT contours image in multiple formats to categorized folder
    fprintf('  Saving H_rr WITHOUT contours plot for file %d...\n', file_idx);
    
    filename_H_rr_no_contours_fig = sprintf('%s_H_rr_NoContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    filename_H_rr_no_contours_jpeg = sprintf('%s_H_rr_NoContours_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
    filename_H_rr_no_contours_eps = sprintf('%s_H_rr_NoContours_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
    
    % Full paths to categorized folder
    H_rr_folder = fullfile(main_output_folder, 'H_rr_NoContours');
    fullpath_H_rr_no_contours_fig = fullfile(H_rr_folder, filename_H_rr_no_contours_fig);
    fullpath_H_rr_no_contours_jpeg = fullfile(H_rr_folder, filename_H_rr_no_contours_jpeg);
    fullpath_H_rr_no_contours_eps = fullfile(H_rr_folder, filename_H_rr_no_contours_eps);
    
    saveas(fig_H_rr_no_contours, fullpath_H_rr_no_contours_fig, 'fig');
    saveas(fig_H_rr_no_contours, fullpath_H_rr_no_contours_jpeg, 'jpeg');
    saveas(fig_H_rr_no_contours, fullpath_H_rr_no_contours_eps, 'eps');
    
    fprintf('    Saved: %s\n', filename_H_rr_no_contours_fig);
    fprintf('    Saved: %s\n', filename_H_rr_no_contours_jpeg);
    fprintf('    Saved: %s\n', filename_H_rr_no_contours_eps);
    
    %% Create H_rr WITH Grayscale Background (no shading)
    if ~isempty(gray_background)
        fprintf('  Creating H_rr WITH grayscale background for file %d...\n', file_idx);
        
        fig_H_rr_gray = figure('Name', sprintf('H_rr with Gray Background - %s (M=%d, n=%d, eps=%.1f)', ...
                             preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                             'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
        
        % Display background
        imagesc([x_min, x_max], [y_min, y_max], gray_background);
        hold on;
        
        % Overlay H_rr field with transparency - NO SHADING
        h_H_rr = pcolor(X, Y, H_rr);
        set(h_H_rr, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        shading interp; 
        colormap(custom_cmap);
        caxis([-0.002, 0.002]);
        colorbar;
        
        % ADD WOUND BOUNDARY AND FIRST CONTOUR TO GRAYSCALE BACKGROUND PLOTS
        % Plot first contour (SDF=epsilon_new) using new epsilon
        C_first = contourc(unique_x, unique_y, SDF, [epsilon, epsilon]);
        [first_contour_x, first_contour_y] = processContourData(C_first);
        if ~isempty(first_contour_x)
            plot(first_contour_x, first_contour_y, 'm-', 'LineWidth', 3, 'DisplayName', 'First Contour');
        end
        
        xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
        ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
        title('$\partial v_{r} / \partial r$','fontsize',16,Interpreter="latex");
        axis equal;
        
        % Save H_rr WITH grayscale background
        filename_H_rr_gray_fig = sprintf('%s_H_rr_WithGrayBackground_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
        filename_H_rr_gray_jpeg = sprintf('%s_H_rr_WithGrayBackground_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
        filename_H_rr_gray_eps = sprintf('%s_H_rr_WithGrayBackground_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
        
        % Full paths to grayscale background folder
        H_rr_gray_folder = fullfile(main_output_folder, 'H_rr_WithGrayBackground');
        fullpath_H_rr_gray_fig = fullfile(H_rr_gray_folder, filename_H_rr_gray_fig);
        fullpath_H_rr_gray_jpeg = fullfile(H_rr_gray_folder, filename_H_rr_gray_jpeg);
        fullpath_H_rr_gray_eps = fullfile(H_rr_gray_folder, filename_H_rr_gray_eps);
        
        saveas(fig_H_rr_gray, fullpath_H_rr_gray_fig, 'fig');
        saveas(fig_H_rr_gray, fullpath_H_rr_gray_jpeg, 'jpeg');
        saveas(fig_H_rr_gray, fullpath_H_rr_gray_eps, 'eps');
        
        fprintf('    Saved with enhanced grayscale background: %s\n', filename_H_rr_gray_jpeg);
        
        close(fig_H_rr_gray);
    end
    
    %% Create individual H_theta_theta 2D plot with contours for each file
    fprintf('  Creating H_theta_theta 2D plot with contours for file %d...\n', file_idx);
    
    fig_H_theta_theta_contours = figure('Name', sprintf('H_theta_theta Field with Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
    
    % Plot H_theta_theta field
    h = pcolor(X, Y, H_theta_theta);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    hold on;
    
    % Store contour handles and labels for legend
    H_theta_theta_contour_handles = [];
    H_theta_theta_contour_labels = {};
    
    % Generate and plot all contours using new epsilon
    for i = 1:length(region_boundaries_new)
        distance = region_boundaries_new(i);
        
        % Extract contour at this SDF level
        C = contourc(unique_x, unique_y, SDF, [distance, distance]);
        
        % Process contour data
        [contour_x, contour_y] = processContourData(C);
        
        % Plot contour
        if ~isempty(contour_x)
            if abs(distance - epsilon) < 1e-6  % Compare with new epsilon
                % First contour - use special color
                h_contour = plot(contour_x, contour_y, 'm-', 'LineWidth', 3);
                H_theta_theta_contour_handles = [H_theta_theta_contour_handles, h_contour];
                H_theta_theta_contour_labels{end+1} = sprintf('First Contour (%.1f μm)', distance);
            else
                % Other contours
                color_idx = mod(i-1, M) + 1;
                h_contour = plot(contour_x, contour_y, 'Color', region_colors(color_idx,:), ...
                         'LineWidth', 2);
                H_theta_theta_contour_handles = [H_theta_theta_contour_handles, h_contour];
                H_theta_theta_contour_labels{end+1} = sprintf('Contour: %.1f μm', distance);
            end
        end
    end
    
    % % Plot boundary (SDF=0) for reference
    % C_boundary = contourc(unique_x, unique_y, SDF, [0, 0]);
    % [boundary_contour_x, boundary_contour_y] = processContourData(C_boundary);
    % if ~isempty(boundary_contour_x)
    %     h_boundary = plot(boundary_contour_x, boundary_contour_y, 'r-', 'LineWidth', 3);
    %     H_theta_theta_contour_handles = [H_theta_theta_contour_handles, h_boundary];
    %     H_theta_theta_contour_labels{end+1} = 'Wound Boundary (0 μm)';
    % end
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('$v_{r}/r$','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save H_theta_theta with contours image as .fig format only
    fprintf('  Saving H_theta_theta with contours plot for file %d...\n', file_idx);
    
    filename_H_theta_theta_contours_fig = sprintf('%s_H_theta_theta_WithContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    fullpath_H_theta_theta_contours_fig = fullfile(main_output_folder, filename_H_theta_theta_contours_fig);
    saveas(fig_H_theta_theta_contours, fullpath_H_theta_theta_contours_fig, 'fig');
    fprintf('    Saved: %s\n', filename_H_theta_theta_contours_fig);
    
    %% Create individual H_theta_theta 2D plot WITHOUT contours for each file
    fprintf('  Creating H_theta_theta 2D plot WITHOUT contours for file %d...\n', file_idx);
    
    fig_H_theta_theta_no_contours = figure('Name', sprintf('H_theta_theta Field WITHOUT Contours - %s (M=%d, n=%d, eps=%.1f)', ...
                           preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                           'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
    
    % Plot H_theta_theta field
    h = pcolor(X, Y, H_theta_theta);
    set(h, 'EdgeColor', 'none');
    colormap(custom_cmap);
    caxis([-0.002, 0.002]);
    colorbar;
    
    xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
    ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
    title('$v_{r}/r$','fontsize',16,Interpreter="latex");
    axis equal;
    
    %% Save H_theta_theta WITHOUT contours image in multiple formats to categorized folder
    fprintf('  Saving H_theta_theta WITHOUT contours plot for file %d...\n', file_idx);
    
    filename_H_theta_theta_no_contours_fig = sprintf('%s_H_theta_theta_NoContours_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
    filename_H_theta_theta_no_contours_jpeg = sprintf('%s_H_theta_theta_NoContours_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
    filename_H_theta_theta_no_contours_eps = sprintf('%s_H_theta_theta_NoContours_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
    
    % Full paths to categorized folder
    H_theta_theta_folder = fullfile(main_output_folder, 'H_theta_theta_NoContours');
    fullpath_H_theta_theta_no_contours_fig = fullfile(H_theta_theta_folder, filename_H_theta_theta_no_contours_fig);
    fullpath_H_theta_theta_no_contours_jpeg = fullfile(H_theta_theta_folder, filename_H_theta_theta_no_contours_jpeg);
    fullpath_H_theta_theta_no_contours_eps = fullfile(H_theta_theta_folder, filename_H_theta_theta_no_contours_eps);
    
    saveas(fig_H_theta_theta_no_contours, fullpath_H_theta_theta_no_contours_fig, 'fig');
    saveas(fig_H_theta_theta_no_contours, fullpath_H_theta_theta_no_contours_jpeg, 'jpeg');
    saveas(fig_H_theta_theta_no_contours, fullpath_H_theta_theta_no_contours_eps, 'eps');
    
    fprintf('    Saved: %s\n', filename_H_theta_theta_no_contours_fig);
    fprintf('    Saved: %s\n', filename_H_theta_theta_no_contours_jpeg);
    fprintf('    Saved: %s\n', filename_H_theta_theta_no_contours_eps);
    
    %% Create H_theta_theta WITH Grayscale Background
    if ~isempty(gray_background)
        fprintf('  Creating H_theta_theta WITH grayscale background for file %d...\n', file_idx);
        
        fig_H_theta_theta_gray = figure('Name', sprintf('H_theta_theta with Gray Background - %s (M=%d, n=%d, eps=%.1f)', ...
                             preprocessed_data(file_idx).filename, M, diff_order, epsilon), ...
                             'NumberTitle', 'off', 'Position', [450, 150, 1000, 800]);
        
        % Display background
        imagesc([x_min, x_max], [y_min, y_max], gray_background);
        hold on;
        
        % Overlay H_theta_theta field with transparency - NO SHADING
        h_H_theta_theta = pcolor(X, Y, H_theta_theta);
        set(h_H_theta_theta, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        shading interp; 
        colormap(custom_cmap);
        caxis([-0.002, 0.002]);
        colorbar;
        
        % ADD WOUND BOUNDARY AND FIRST CONTOUR TO GRAYSCALE BACKGROUND PLOTS 
        % Plot first contour (SDF=epsilon_new) using new epsilon
        C_first = contourc(unique_x, unique_y, SDF, [epsilon, epsilon]);
        [first_contour_x, first_contour_y] = processContourData(C_first);
        if ~isempty(first_contour_x)
            plot(first_contour_x, first_contour_y, 'm-', 'LineWidth', 3, 'DisplayName', 'First Contour');
        end
        
        xlabel('$x\;(\mu m)$','fontsize',14,Interpreter="latex"); 
        ylabel('$y\;(\mu m)$','fontsize',14,Interpreter="latex");
        title('$v_{r}/r$','fontsize',16,Interpreter="latex");
        axis equal;
        
        % Save H_theta_theta WITH grayscale background
        filename_H_theta_theta_gray_fig = sprintf('%s_H_theta_theta_WithGrayBackground_M%d_n%d_eps%.1f.fig', name_only, M, diff_order, epsilon);
        filename_H_theta_theta_gray_jpeg = sprintf('%s_H_theta_theta_WithGrayBackground_M%d_n%d_eps%.1f.jpeg', name_only, M, diff_order, epsilon);
        filename_H_theta_theta_gray_eps = sprintf('%s_H_theta_theta_WithGrayBackground_M%d_n%d_eps%.1f.eps', name_only, M, diff_order, epsilon);
        
        % Full paths to grayscale background folder
        H_theta_theta_gray_folder = fullfile(main_output_folder, 'H_theta_theta_WithGrayBackground');
        fullpath_H_theta_theta_gray_fig = fullfile(H_theta_theta_gray_folder, filename_H_theta_theta_gray_fig);
        fullpath_H_theta_theta_gray_jpeg = fullfile(H_theta_theta_gray_folder, filename_H_theta_theta_gray_jpeg);
        fullpath_H_theta_theta_gray_eps = fullfile(H_theta_theta_gray_folder, filename_H_theta_theta_gray_eps);
        
        saveas(fig_H_theta_theta_gray, fullpath_H_theta_theta_gray_fig, 'fig');
        saveas(fig_H_theta_theta_gray, fullpath_H_theta_theta_gray_jpeg, 'jpeg');
        saveas(fig_H_theta_theta_gray, fullpath_H_theta_theta_gray_eps, 'eps');
        
        fprintf('    Saved with enhanced grayscale background: %s\n', filename_H_theta_theta_gray_jpeg);
        
        close(fig_H_theta_theta_gray);
    end
    
    % Close all figures for this file iteration (except summary plots which are created later)
    close all;
    
    % Display summary information for this file
    fprintf('  File %d summary (M=%d, n=%d, eps=%.1f):\n', file_idx, M, diff_order, epsilon);
    for i = 1:M
        if ~isnan(rMean_new(i))
            fprintf('    Region %d: r=%.1f, H_rr=%.4e, H_θθ=%.4e, div=%.4e\n', ...
                    i, rMean_new(i), radial_strain_rate_Mean_new(i), ...
                    tangential_strain_rate_Mean_new(i), divergence_Mean_new(i));
        end
    end
end

%% Save gradient data to file
fprintf('\nSaving gradient data to file...\n');
gradient_filename = fullfile(main_output_folder, sprintf('gradient_data_M%d_n%d_eps%.1f.mat', M, diff_order, epsilon));
save(gradient_filename, 'gradient_data', 'M', 'diff_order', 'epsilon');
fprintf('  Saved gradient data to: %s\n', gradient_filename);

%% Create summary plots
% fprintf('\nCreating summary plots with M = %d, n = %d, epsilon = %.1f...\n', M_new, diff_order, epsilon_new);
% 
% % Store handles for summary figures
% summary_fig_handles = [];
% 
% % Plot 1: Radial Strain Rate vs Distance to Center
% fig1 = figure('Name', sprintf('Radial Strain Rate vs Distance (M=%d, n=%d, eps=%.1f)', M_new, diff_order, epsilon_new), ...
%        'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
% hold on;
% 
% for file_idx = 1:num_files
%     if ~isempty(all_rMean_new{file_idx}) && ~isempty(all_radial_strain_rate_Mean_new{file_idx})
%         % Remove NaN values
%         valid_idx = ~isnan(all_rMean_new{file_idx}) & ~isnan(all_radial_strain_rate_Mean_new{file_idx});
%         if any(valid_idx)
%             plot(all_rMean_new{file_idx}(valid_idx), all_radial_strain_rate_Mean_new{file_idx}(valid_idx), ...
%                  'o-', 'Color', colors(file_idx,:), 'LineWidth', 2, 'MarkerSize', 6, ...
%                  'DisplayName', legendLabels{file_idx});
%         end
%     end
% end
% 
% % Add horizontal zero line
% plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
% 
% xlim([0 1500]);
% xlabel("distance to center $(\mu m)$", 'fontsize', 16, Interpreter="latex")
% ylabel("${{\partial v_{r}(r)}/{\partial r}}\;({min}^{-1})$", 'fontsize', 16, Interpreter="latex")
% legend('show', 'Location', 'best')
% grid on
% 
% summary_fig_handles = [summary_fig_handles, fig1];
% 
% % Plot 2: Tangential Strain Rate vs Distance to Center
% fig2 = figure('Name', sprintf('Tangential Strain Rate vs Distance (M=%d, n=%d, eps=%.1f)', M_new, diff_order, epsilon_new), ...
%        'NumberTitle', 'off', 'Position', [200, 100, 800, 600]);
% hold on;
% 
% for file_idx = 1:num_files
%     if ~isempty(all_rMean_new{file_idx}) && ~isempty(all_tangential_strain_rate_Mean_new{file_idx})
%         % Remove NaN values
%         valid_idx = ~isnan(all_rMean_new{file_idx}) & ~isnan(all_tangential_strain_rate_Mean_new{file_idx});
%         if any(valid_idx)
%             plot(all_rMean_new{file_idx}(valid_idx), all_tangential_strain_rate_Mean_new{file_idx}(valid_idx), ...
%                  's-', 'Color', colors(file_idx,:), 'LineWidth', 2, 'MarkerSize', 6, ...
%                  'DisplayName', legendLabels{file_idx});
%         end
%     end
% end
% 
% % Add horizontal zero line
% plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
% 
% xlim([0 1500]);
% xlabel("distance to center $(\mu m)$", 'fontsize', 16, Interpreter="latex")
% ylabel("${{v_{r}}/r}\;({min}^{-1})$", 'fontsize', 16, Interpreter="latex")
% legend('show', 'Location', 'best')
% grid on
% 
% summary_fig_handles = [summary_fig_handles, fig2];
% 
% % Plot 3: Divergence vs Distance to Center
% fig3 = figure('Name', sprintf('Divergence vs Distance (M=%d, n=%d, eps=%.1f)', M_new, diff_order, epsilon_new), ...
%        'NumberTitle', 'off', 'Position', [300, 100, 800, 600]);
% hold on;
% 
% for file_idx = 1:num_files
%     if ~isempty(all_rMean_new{file_idx}) && ~isempty(all_divergence_Mean_new{file_idx})
%         % Remove NaN values
%         valid_idx = ~isnan(all_rMean_new{file_idx}) & ~isnan(all_divergence_Mean_new{file_idx});
%         if any(valid_idx)
%             plot(all_rMean_new{file_idx}(valid_idx), all_divergence_Mean_new{file_idx}(valid_idx), ...
%                  '^-', 'Color', colors(file_idx,:), 'LineWidth', 2, 'MarkerSize', 6, ...
%                  'DisplayName', legendLabels{file_idx});
%         end
%     end
% end
% 
% % Add horizontal zero line
% plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
% 
% xlim([0 1500]);
% xlabel("distance to center $(\mu m)$", 'fontsize', 16, Interpreter="latex")
% ylabel("div$(\mbox{\boldmath $v$)}" + "\;({min}^{-1})$", 'fontsize', 16, Interpreter="latex")
% legend('show', 'Location', 'best')
% grid on
% 
% summary_fig_handles = [summary_fig_handles, fig3];
% 
% %% Save summary plots as .fig format
% fprintf('Saving summary plots...\n');
% 
% summary_filenames = {
%     sprintf('Radial_Strain_Rate_M%d_n%d_eps%.1f.fig', M_new, diff_order, epsilon_new),
%     sprintf('Tangential_Strain_Rate_M%d_n%d_eps%.1f.fig', M_new, diff_order, epsilon_new),
%     sprintf('Divergence_M%d_n%d_eps%.1f.fig', M_new, diff_order, epsilon_new)
% };
% 
% for i = 1:length(summary_fig_handles)
%     fig = summary_fig_handles(i);
%     filename = fullfile(main_output_folder, summary_filenames{i});
%     saveas(fig, filename, 'fig');
%     fprintf('  Saved: %s\n', summary_filenames{i});
% end

fprintf('\nReanalysis completed successfully!\n');
fprintf('Created and saved %d individual velocity field plots with contours in main folder: %s\n', num_files, main_output_folder);
fprintf('Created and saved all 2D plots WITHOUT contours in categorized subfolders\n');
fprintf('Created and saved all 2D plots (except velocity field) WITH ENHANCED grayscale background in categorized subfolders\n');
fprintf('Enhanced grayscale background plots now include wound boundary and first contour\n');
fprintf('NOTE: Removed du_dx, du_dy, dv_dx, dv_dy plots as requested\n');
fprintf('Created and saved gradient data to: %s\n', gradient_filename);
% fprintf('Created and saved 3 summary plots\n');
fprintf('All plots use M = %d, differential order n = %d, and epsilon = %.1f\n', M, diff_order, epsilon);

%% Helper function to process contour data
function [x_out, y_out] = processContourData(C)
    % Extract contour lines from contour matrix C
    x_out = [];
    y_out = [];
    
    i = 1;
    while i < size(C, 2)
        level = C(1, i);
        n_points = C(2, i);
        
        if n_points > 0
            % Extract this contour line
            contour_x = C(1, i+1:i+n_points);
            contour_y = C(2, i+1:i+n_points);
            
            % Keep the longest contour (in case there are multiple)
            if length(contour_x) > length(x_out)
                x_out = contour_x;
                y_out = contour_y;
            end
        end
        
        i = i + n_points + 1;
    end
    
    % Ensure output is column vectors
    x_out = x_out(:);
    y_out = y_out(:);
end

%% Natural sorting function
function sortedFiles = sort_natural(files)
    names = {files.name};
    % Extract numbers from filenames
    nums = zeros(length(names), 1);
    for i = 1:length(names)
        num_str = regexp(names{i}, '\d+', 'match');
        if ~isempty(num_str)
            nums(i) = str2double(num_str{1});
        else
            nums(i) = inf;
        end
    end
    [~, indices] = sort(nums);
    sortedFiles = files(indices);
end