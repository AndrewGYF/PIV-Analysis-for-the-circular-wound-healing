# PIV Analysis for the circular wound healing
This repository includes the MATLAB script, the raw data, and the TIF. images that are required to reproduce PIV results regarding the circular wound healing. 
## Numerical Implementation
Provided the position data, i.e., $(x,y)$, the velocity data at each position, i.e., $(v_{x},v_{y})$, as well as the wound boundary in terms of the sign distance function (SDF) in the MATLAB file, we want to create the tangential and circumferential $2D$ field plots such that we can characterize and visualize its both local and global behaviors; that said, stretching or compression. In this regard, we can compute the velocity gradient $\mathbf{\nabla_{\mathbf{x}}}\mathbf{v}$, and hence the radial strain rate and the tangential strain rate:
```matlab
% Initialize gradient matrices
    [ny, nx] = size(U_grid);
    du_dx = zeros(ny, nx);
    du_dy = zeros(ny, nx);
    dv_dx = zeros(ny, nx);
    dv_dy = zeros(ny, nx);

    % Calculate gradients with First Contour SDF-aware method using epsilon
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
```
One can refer to the PIV guide in the repository (coming soon). In addition, we compute the divergence plots for each time frame in order to observe the wound healing situation as time evolves.
## User Manual for the MATLAB Script
One should make sure that both the MATLAB file (summarizing position data, velocity data, as well as the wound boundary information in terms of the sign distance function (SDF)), as well as the TIF. images are in the current folder while running the script. The `dir` function is going to detect and hence employ those TIF. images in the current file:
```matlab
tifFiles = dir('*.tif');
```
Users will need to define the parameters in the command window:
```matlab
%% User-defined input for new parameters
M = input('Enter the bin number M: ');
diff_order = input(['Enter the differential order n for' ...
    ' central difference (e.g., 1, 2, 3): ']);
epsilon = input('Enter the epsilon value: ');
```
In our simulation, we set the bin number to be nine, the differential order to be 4, and $\varepsilon$ to be 20 $\mu m$, where the wound region is divided into `M` regions of the same width, `diff_order` suggests the increment for the numerical derivatives, and $\varepsilon$ suggests the distance by which a contour is pushed back from the wound boundary. As users finish defining those parameters, the script will create those $2D$ plots, the tangential strain rate field, the circumferential strain rate field, as well as the divergence field for each time frame. Those TIF. images will be used as a background for those $2D$ plots. 
