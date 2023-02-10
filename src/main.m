clc();
clear;
close all;

% == User Settings ==

% Parameters:
% |
% |__ sample_number (int, [1, 2, 3, 4, 5]): The Default Samples.
% |       - 1: Tip-loaded Cantilever Beam
% |       - 2: Mid-loaded Cantilever Beam
% |       - 3: MBB Beam
% |       - 4: Michell
% |       - 4: Bridge
% |
% |__ nel_x (int): Number of element @ x-axis.
% |
% |__ nel_y (int): Number of element @ y-axis.
% |
% |__ volume_fraction (float): The volume fraction of final topology
% |                            structure you want.
% |
% |__ penalization (float): The penalization can reduce the effect of gray
% |                         element on the stiffness matrix.
% |
% |__ r_min (float): The minimum sensitivity (or density) filter radius.
% |
% |__ filter_class(int, [1, 2]):
% |      - 1: Sensitivity filter
% |      - 2: Density filter
% |
% |__ E_0 (float): The Young's modulus of material
% |
% |__ E_min (float): The minimum Young's modulus of stiffness matrix
% |
% |__ nu (float): The Poisson's ratio of material
% |
% |__ is_plot_result (boolean): It controls whether topology structure are
% |                             are drawn.
% |
% |__ is_show_iteration_result (boolean): It controls whether to print the
% |                                       compliance, gray index, etc.
% |
% |__ is_use_projection_function (boolean): It controls whether use the
% |                                         projection function.
% |
% |__ is_timekeeping (boolean): It controls whether to is_timekeeping.
% |
% |__ update_method (str, ["oc", "mma"]): The update method used by the
%                                         topology optimization algorithm.

sample_number = 1;

nel_x = 60;
nel_y = 30;

volume_fraction = 0.3;
penalization = 3;
r_min = 2;
filter_class = 1;

E_0 = 1; 
E_min = 1e-9;
nu = 0.3;

is_plot_result = true;
is_show_iteration_result = true;
is_use_projection_function = true;
is_timekeeping = true;

update_method = "mma";

% == User Settings ==

addpath(genpath("update_methods/"));
addpath(genpath("assets/samples/"));
addpath(genpath("utils"));

switch sample_number
    case 1
        [F, U, freedofs] = tip_loaded_cantilever_beam(nel_x, nel_y);
    case 2
        [F, U, freedofs] = mid_loaded_cantilever_beam(nel_x, nel_y);
    case 3
        [F, U, freedofs] = MBB_beam(nel_x, nel_y);
    case 4
        [F, U, freedofs] = michell(nel_x, nel_y);
    case 5
        [F, U, freedofs] = bridge(nel_x, nel_y);
end


x = topology_optimization( ...
    nel_x, ...
    nel_y, ...
    volume_fraction, ...
    3, ...
    2, ...
    1, ...
    E_0, ...
    E_min, ...
    nu, ...
    F, ...
    U, ...
    freedofs, ...
    is_plot_result, ...
    is_show_iteration_result, ...
    update_method, ...
    is_use_projection_function ...
);

if is_timekeeping
    disp(toc);
end