clc();
clear;
close all;

% == User Settings ==

% Parameters:
% |
% |__ sample_number (int, [1, 2, 3, 4, 5, 6]): The Default Samples.
% |       - 1: Tip-loaded Cantilever Beam
% |       - 2: Mid-loaded Cantilever Beam
% |       - 3: MBB Beam
% |       - 4: Michell
% |       - 5: Bridge
% |       - 6: Inverter
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
% |__ sensitivity_filter_class (int, [1, 2]):
% |      - 1: Sensitivity filter
% |      - 2: Density filter
% |
% |__ densign_variable_filter_class (int, [1, 2]):
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
% |                                       topology optimization algorithm.
% |
% |__ max_iteration_number (int): The max iteration number.

sample_number = 6;

nel_x = 120;
nel_y = 60;

volume_fraction = 0.25;
penalization = 3;
r_min = 2.5;
sensitivity_filter_class = 2;
densign_variable_filter_class = 2;

E_0 = 1;
E_min = 1e-9;
nu = 0.3;
 
is_plot_result = false;
is_show_iteration_result = true;
is_use_projection_function = true;
is_timekeeping = false;

update_method = "mma";

max_iteration_number = 1000;

% == User Settings ==

addpath(genpath("update_methods/"));
addpath(genpath("assets/samples/"));
addpath(genpath("utils"));

is_spring = false;
spring_k = [];
input_element = [];
output_element = [];
switch sample_number
    case 1
        [F, U, free_dofs] = tip_loaded_cantilever_beam(nel_x, nel_y);
    case 2
        [F, U, free_dofs] = mid_loaded_cantilever_beam(nel_x, nel_y);
    case 3
        [F, U, free_dofs] = MBB_beam(nel_x, nel_y);
    case 4
        [F, U, free_dofs] = michell(nel_x, nel_y);
    case 5
        [F, U, free_dofs] = bridge(nel_x, nel_y);
    case 6
        [F, U, free_dofs, spring_k, input_element, output_element] = inverter(nel_x, nel_y);
        is_spring = true;
end

[x, GA, MA] = topology_optimization( ...
    nel_x, ...
    nel_y, ...
    volume_fraction, ...
    penalization, ...
    r_min, ...
    sensitivity_filter_class, ...
    densign_variable_filter_class, ...
    E_0, ...
    E_min, ...
    nu, ...
    F, ...
    U, ...
    free_dofs, ...
    is_plot_result, ...
    is_show_iteration_result, ...
    update_method, ...
    is_use_projection_function, ...
    max_iteration_number, ...
    spring_k, ...
    is_spring, ...
    input_element, ...
    output_element ...
);

if is_timekeeping
    disp(toc);
end

if is_spring
    fprintf("GA: %6.4f, MA: %6.4f \n", GA, MA);
end