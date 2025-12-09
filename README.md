# PIV Analysis for the circular wound healing
This repository includes the MATLAB script, the raw data, and the TIF. images that are required to reproduce PIV results regarding the circular wound healing. 
## Numerical Implementation
Provided the position data, i.e., $(x,y)$, the velocity data at each position, i.e., $(v_{x},v_{y})$, as well as the wound boundary in terms of the sign distance function (SDF) in the MATLAB file, we want to create the tangential and circumferential $2D$ field plots such that we can characterize and visualize its both local and global behaviors; that said, stretching or compression. In this regard, we can compute the velocity gradient $\mathbf{\nabla_{\mathbf{x}}}\mathbf{v}$, and hence the radial strain rate and the tangential strain rate. One can refer to the PIV guide in the repository (coming soon). In addition, we compute the divergence plots for each time frame in order to observe the wound healing situation as time evolves.
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
In our simulation, we set the bin number to be nine, the differential order to be 4, and $\varepsilon$ to be 20 $\mu m$. As users finish defining those parameters, the script will create those $2D$ plots, the tangential strain rate field, the circumferential strain rate field, as well as the divergence field for each time frame. Those TIF. images will be used as a background for those $2D$ plots. 
