%  SETPLOT3 sets user defined plotting parameters
%
%      User defined Matlab script for setting various Clawpack plotting
%      parameters.  This script is called by PLOTCLAW3.  A default
%      version of this script can be found in claw/matlab/setplot3.m and
%      copied to users working directory and modifed to set things up
%      differently.
%
%      Parameters that can be set with SETPLOT3
%
%        OutputFlag        - set to 'ascii' (default) to read ascii output
%                            files, and 'hdf' to read hdf files.
%        PlotType          - type of plot to produce:
% 			     - 1 = pcolor on slices (with optional contours
% 			         and isosurfaces)
% 			     - 2 = contour lines in 3d on white slices
% 			     - 3 = Schlieren plot on slices
% 			     - 4 = scatter plot of q vs. r
%
%        mq                  - which component of q to plot
%        UserVariable        - Set to 1 to specify a user defined variable.
%        UserVariableFile    - name of m-file mapping data to q
%        MappedGrid          - set to 1 if mapc2p.m exists for nonuniform grid
%        MaxFrames           - max number of frames
%        MaxLevels           - max number of AMR levels
%        PlotData            - Data on refinement level k is plotted only if
%                              PlotData(k) == 1
%        PlotGrid            - PLot grid lines on level k is PlotGrid(k) /= 0
%        PlotGridEdges       - Plot 2d patch borders if PlotGridEdges(k) /= 0
%        PlotCubeEdges       - Plot 3d patch cubes if PlotCubeEdges(k) /= 0
%        ContourValues       - Set to desired contour values, or [] for no ...
% 	                     lines.
%        xSliceCoords        - vector of x slice constants
%        ySliceCoords        - vector of y slice constants
%        zSliceCoords        - vector of z slice constants
%        x0,y0,z0            - center for scatter plots.
%        ScatterStyle        - symbols for scatter plots.
%        LineStyle           - same as ScatterStyle.
%        IsosurfValues       - constants for isosurfaces
%        IsosurfColors       - colors for isosurfaces.
%        UserMap1d           - set to 1 if 'map1d' file exists.
%
%      All parameters can be modified by typing 'k' at the PLOTCLAW3 prompt.
%
%      See also PLOTCLAW3, setPlotGrid, setPlotGridEdges.


setopengl;
setviews;  % set viewpoints so that view(xSlice), for example, can be used.

% OutputFlag = 'ascii';      % default.

PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours,
			     % and isosurfaces)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % max number of AMR levels

PlotData =  [1 1 1 0 0 0];       % Data on refinement level k is plotted only
			         % if k'th component is nonzero
PlotGrid =  [0 0 0 0 0 0];       % Plot grid lines on each level?
PlotGridEdges = [0 0 0 0 0 0];   % Plot edges of patches of each grid at
                                 % this level on slices?  Shows edges of MPI
                                 % subdomains with ASCII output.
PlotCubeEdges = [0 0 0 0 0 0];   % Plot edges of cube of refinement patch at
                                 % this level?


% ---------------------------------------------------------------------
% The next three parameters are vectors of x,y,z coordinates of 2d slices
% to be displayed for PlotType = 1,2,3.
% Empty ==> no slices in that direction.

  xSliceCoords = 0.6;
  ySliceCoords = 0.6;
  zSliceCoords = 0.425;


% ---------------------------------------------------------------------
% ContourValues is a vector of values used to draw contour lines.  If
% The valid settings for this parameter are identical to those used by the
% Matlab contour plotting routine.  See also CONTOUR.
% If ContourValues is the empty matrix, no contour lines will be drawn.
  ContourValues = [.1:.1:.9];

% Isosurfaces.  If empty, no isosurfaces will be drawn.
% Note:  Lengths of SurfTransparency and SurfColors must greater than or
%        equal to the length of SurfValues..
IsosurfValues =  [];     % Plot surfaces at q = surfValue(i).

IsosurfAlphas =  [];     % Transparency of each surface
                                     % (0=clear; 1=opaque)
                                     % NOTE: Your system must be able to
                                     % use the OpenGL Renderer.

IsosurfColors =  [];                 % Colors for each surface.
                                     % Set to 'q' to get colors from the
                                     % current colormaps.  Use STRVCAT
                                     % to get multiple colors, i.e.
                                     % strvcat('b','r','g','y');



% ---------------------------------------------------------------------
% For PlotType = 4 (Scatter plot)
% plot q(r) vs. r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2);
% Use symbol ScatterStyle{k} at refinement level k.
  x0 = 0.5;
  y0 = 0.5;
  z0 = 0.5;
  ScatterStyle = setplotstyle('o','x','.','s','v','^');
