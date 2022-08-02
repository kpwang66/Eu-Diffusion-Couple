% Solves the 2-domain composite material 1D diffusion equation with an implicit finite difference scheme
% 
% There is in an initial stationary reference frame (SRF) followed by a moving reference frame (MRF) when the interface begins to move
% 
% FD4 updates: uses 1D data that is extracted from inverse-radon prism % technique

clear, clc, close all;

this_script = matlab.desktop.editor.getActiveFilename;
this_dir = split(this_script, '\');
this_dir = strjoin(this_dir(1:end-1),'\');
cd(this_dir);

Eu_data_old = load('.\Eu_data.mat');
Eu_data_new = load('.\IRad_export_a-90_merged_cyl_b\IRad_export_a-90_merged_cyl_b_1D.mat');


Eu_1D_new = Eu_data_new.Eu1d_dataset./100;
Eu0_new = Eu_1D_new(:, 1);

n_px_flux = 10;
l_flux_raw = zeros(12, 1);


% C0_tot: total averaged Eu concentration calculated from extrapolated initial data
% x_full: vector of x positions from approx 0
% Eu0_full: vector of Eu0 from time 0, including extrapolated portion
% t_vec: vector of time points (0 to 242)
% Eu_t_mat_raw: raw Eu 
% x_vec: vector of x positions from Anton's data (mm)
% int_vel_mm_min: interface velocity from time point 5 (88 min) to end (220 mm)

MRF_t_start_min = 110;

% Physical parameters
shared_params.Ds   =   1.9e-10;       %   Mass diffusivity of solid [m2/s]
shared_params.Dl   =  2.5e-10*1.1;       %   Mass diffusivity of liquid [m2/s]
shared_params.K = .45; % Partition coefficient
shared_params.Kfit = Eu_data_old.Kfit;


%% SRF Parameters
% SRF Spacing
int_loc_0_mm = 13; % 13mm
int_loc_0_m = int_loc_0_mm*1e-3;
SRF_params.x_div      =   200;        %   Approx number of gridpoints in x-direction
SRF_params.gamma = int_loc_0_mm/Eu_data_old.x_vec(end);
SRF_x_vec_A = linspace(1e-3*Eu_data_old.x_vec(1), int_loc_0_m, round(SRF_params.x_div*SRF_params.gamma));
shared_params.dx = mean(unique(diff(SRF_x_vec_A)));
SRF_x_vec_B = int_loc_0_m : shared_params.dx : fix(1e-3*Eu_data_old.x_vec(end)/shared_params.dx)*shared_params.dx;
SRF_params.nx = length(SRF_x_vec_A) + length(SRF_x_vec_B) - 1;

% SRF initial condition (importing)
Eu0_interp_A = 1e-2*interp1(Eu_data_old.x_vec./1e3, Eu0_new, SRF_x_vec_A);
Eu0_interp_B = 1e-2*interp1(Eu_data_old.x_vec./1e3, Eu0_new, SRF_x_vec_B);

% SRF special nodes
SRF_params.nx_gamma_s = length(SRF_x_vec_A);
SRF_params.nx_gamma_sr = length(SRF_x_vec_A) + 1;
SRF_params.nx_gamma_ll = length(SRF_x_vec_A) + 2;
SRF_params.nx_gamma_l = length(SRF_x_vec_A) + 3;

SRF_params.int_vel = 0;

%% Time parameters
SRF_params.t_end = MRF_t_start_min*60; % 
MRF_params.t_end = 242*60;

shared_params.dt      =   22*60/10;      %   Timestep [s]

SRF_params.nt      =   fix(SRF_params.t_end/shared_params.dt);        %   Number of timesteps to compute
MRF_params.nt      =   fix((MRF_params.t_end-SRF_params.t_end)/shared_params.dt);

%% MRF Parameters 
% MRF Spacing
origin_start_m = int_loc_0_m;
total_frame_displc_mm = (Eu_data_old.t_vec(end) - SRF_params.t_end/60)*Eu_data_old.int_vel_mm_min;
total_frame_displc_m = total_frame_displc_mm*1e-3;

RF_l_abs_m = SRF_x_vec_A(1);
RF_r_abs_m = SRF_x_vec_B(end) - total_frame_displc_m;

MRF_params.RF_l_m = -abs(origin_start_m - RF_l_abs_m);
MRF_params.RF_r_m = abs(origin_start_m - RF_r_abs_m);
MRF_params.int_vel = Eu_data_old.int_vel_mm_min./60000; % conversion from mm/min to m/s

MRF_params.dom_nx_A = fix((origin_start_m - RF_l_abs_m)/shared_params.dx); % number of elements on solid side NOT INCLUDING FICTIOUS ELEMENTS
MRF_params.dom_nx_B = fix((RF_r_abs_m - origin_start_m)/shared_params.dx); % number of elements on liq side NOT INCLUDING FICTIOUS ELEMENTS
MRF_params.dom_lb = -MRF_params.dom_nx_A*shared_params.dx;
MRF_params.dom_rb = MRF_params.dom_nx_B*shared_params.dx;


MRF_x_vec_A = SRF_x_vec_A - origin_start_m;
MRF_x_vec_B = SRF_x_vec_B(1:MRF_params.dom_nx_B) - origin_start_m;

MRF_params.nx = length(MRF_x_vec_A) + length(MRF_x_vec_B) - 2;

MRF_params.nx_gamma_s = length(MRF_x_vec_A);
MRF_params.nx_gamma_sr = length(MRF_x_vec_A) + 1;
MRF_params.nx_gamma_ll = length(MRF_x_vec_A) + 2;
MRF_params.nx_gamma_l = length(MRF_x_vec_A) + 3;

for k = 1:6
    smoothed = 1e-2*smooth(Eu_1D_new(:, k), 10);
    l_D_raw(k) = mean(smoothed(1:5));
end

for k = 7:12
    smoothed = 1e-2*Eu_1D_new(:, k);
    mig_time = (k-1)*22*60 - SRF_params.t_end;
    l_bound_x  = 1000*shift_x(MRF_x_vec_A(1), origin_start_m, MRF_params.int_vel, mig_time);
    [minValue,closestIndex] = min(abs(Eu_data_old.x_vec-l_bound_x));
    l_D_raw(k) = smoothed(closestIndex);
end

l_bound_Dirichlet_spline = spline([0:11]*22*60, l_D_raw);

% keyboard


% for k = 7:12
%     smoothed = smooth(Eu_1D_new(:, k));
%     mig_time = k*22 - MRF_t_start_min;
%     l_bound_x = shift_x(0, origin_start_m, MRF_params.int_vel, mig_time);
%     
%     keyboard
%     l_flux_raw(k) = (smoothed(n_px_flux) - smoothed(1))./Eu_data_old.x_vec(n_px_flux)*1e-3 - Eu_data_old.x_vec(1)*1e-3;
% end

%%
shared_params.s2s = shared_params.Ds*shared_params.dt/shared_params.dx^2;
shared_params.s2l = shared_params.Dl*shared_params.dt/shared_params.dx^2;
shared_params.s1 = MRF_params.int_vel*shared_params.dt/shared_params.dx/2;

% Setup initial SRF temperature profile
T       =   [Eu0_interp_A'; 0; 0; Eu0_interp_B'];
T(SRF_params.nx_gamma_sr) = Eu0_interp_B(2);
T(SRF_params.nx_gamma_ll) = Eu0_interp_A(end-1);
T(SRF_params.nx_gamma_l)  = Eu0_interp_B(1);

% T(find(abs(x)<=params.W/2)) = Tmagma;
SRF_Toft = cell(SRF_params.nt, 2);
SRF_Toft{1, 1} = Eu0_interp_A;
SRF_Toft{1, 2} = Eu0_interp_B;

% SRF Boundary conditions
B1_SRF.type     =   'D'; 
B1_SRF.val      =   0; 

B2_SRF.type     =   'N';
B2_SRF.val      =   0;

B1_SRF.loc = 'L';
B2_SRF.loc = 'R';

% MRF Boundary conditions
B1_MRF.type     =   'D'; 
B1_MRF.val      =   0; 

B2_MRF.type     =   'N';
B2_MRF.val      =   0;


B1_MRF.loc = 'L';
B2_MRF.loc = 'R';

SRF_params = appendstruct(SRF_params, shared_params);
MRF_params = appendstruct(MRF_params, shared_params);


% First timeslice
time    =   0;
% time = MRF_t_start_min*60 + params.dt;
SRF_t_vec = zeros(1, SRF_params.nt);
MRF_t_vec = zeros(1, MRF_params.nt);

% nx = SRF_params.nx;
% s2s = shared_params.s2s;
% s2l = SRF_params.s2l;
% s1 = SRF_params.s1;
SRF_nt = SRF_params.nt;
dt = shared_params.dt;

MRF_nt = MRF_params.nt;

% Plotting and outputting video
M(SRF_params.nt + MRF_params.nt) = struct('cdata',[],'colormap',[]);
m_i = 1;
V = VideoWriter('C:\Users\Kerry\Dropbox\UMN\Research\Eu Segregation\Vid_FD4.avi');
V.FrameRate = 24;
open(V);

figure(1), clf
plot(1000*SRF_x_vec_A, 100*T(1:length(SRF_x_vec_A)), 'r', 'Linewidth', 2); hold on;
plot(1000*SRF_x_vec_B, 100*T(end-length(SRF_x_vec_B)+1:end), 'r', 'Linewidth', 2);
%     plot(SRF_x_vec_B(1), Tnew(SRF_params.nx_gamma_l), 'ro');

hXLabel = xlabel('Position [mm]');
hYLabel = ylabel('Eu Mole %');
hTitle = title(['Concentration evolution after ',num2str(time/60, '%3.1f'),' minutes']);

set([hXLabel, hYLabel]  , ...
'FontSize'   , 15          );

set( hTitle                    , ...
'FontSize'   , 17          , ...
'FontWeight' , 'bold'      );

ylim([1 5]);
xlim([7 27]);
drawnow;
frame = getframe(gcf);
writeVideo(V,frame);

time = time + SRF_params.dt;

%% SRF Timestep loop
for n=2: SRF_nt  

    % Assemble A for SRF
    A = Assemble_A(SRF_x_vec_A, SRF_x_vec_B, SRF_params);
    Af = full(A);
    
    b = T;

    [A, b] = inputIC(A, b, T, SRF_params);

    % Implement boundary conditions
    B1_SRF.val = ppval(l_bound_Dirichlet_spline, time);
    [A, b] = inputBC(B1_SRF, A, b, T, SRF_params);
    [A, b] = inputBC(B2_SRF, A, b, T, SRF_params);

    % Implicit solve step
    Tnew = A\b;
    
    % Update temperature and time
    SRF_Toft{n, 1}  = Tnew(1:length(SRF_x_vec_A));
    SRF_Toft{n, 2}  = Tnew(end-length(SRF_x_vec_B)+1:end);
    T           =   Tnew;
    time        =   time+dt;
    SRF_t_vec (n) = time;
        
    % Plot solution
    figure(1), clf
    plot(1000*SRF_x_vec_A, 100*Tnew(1:length(SRF_x_vec_A)), 'r', 'Linewidth', 2); hold on;
    plot(1000*SRF_x_vec_B, 100*Tnew(end-length(SRF_x_vec_B)+1:end), 'r', 'Linewidth', 2);
    
    hXLabel = xlabel('Position [mm]');
    hYLabel = ylabel('Eu Mole %');
    hTitle = title(['Concentration evolution after ',num2str(time/60, '%3.1f'),' minutes']);
    
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 15          );

    set( hTitle                    , ...
    'FontSize'   , 17          , ...
    'FontWeight' , 'bold'      );

    ylim([1 5]);
    xlim([7 27]);
    drawnow;

    frame = getframe(gcf);
    writeVideo(V,frame);
end

MRF_Toft = cell(MRF_params.nt, 2);
MRF_Toft{1, 1} = SRF_Toft{end, 1};
MRF_Toft{1, 2} = SRF_Toft{end, 2};

MRF_T = T;
MRF_T_A = T(1:MRF_params.nx_gamma_s);
MRF_T_B = T(MRF_params.nx_gamma_l : (MRF_params.nx_gamma_l + MRF_params.dom_nx_B));

% Setup initial temperature profile for MRF
T       =   T(1: (length(MRF_x_vec_A) + length(MRF_x_vec_B) + 2));

fprintf('Entering MRF loop\n')
%% MRF Timestep loop
for n = 2: MRF_params.nt

    % Assemble A for MRF
    A = Assemble_A(MRF_x_vec_A, MRF_x_vec_B, MRF_params);
    Af = full(A);
    b = T;

    % Implement interface conditions
    [A, b] = inputIC(A, b, T, MRF_params);

    % Implement boundary conditions
    B1_MRF.val = ppval(l_bound_Dirichlet_spline, time);
%     B1_MRF.val = .1;
    [A, b] = inputBC(B1_MRF, A, b, T, MRF_params);
    [A, b] = inputBC(B2_MRF, A, b, T, MRF_params);

    % Implicit solve step
    Tnew = A\b;
    
    mig_time = time - SRF_params.t_end;
    
        
    % Plot solution
    figure(1), clf;
    
    shsc_x_vec_A = 1000*shift_x(MRF_x_vec_A, origin_start_m, MRF_params.int_vel, mig_time);
    shsc_x_vec_B = 1000*shift_x(MRF_x_vec_B, origin_start_m, MRF_params.int_vel, mig_time);
    
    sc_MRF_Tnew_A = 100*Tnew(1:length(MRF_x_vec_A));
    sc_MRF_Tnew_B = 100*Tnew(end-length(MRF_x_vec_B)+1:end);
    plot(shsc_x_vec_A, sc_MRF_Tnew_A, 'r', 'Linewidth', 2); hold on;
    plot(shsc_x_vec_B, sc_MRF_Tnew_B, 'r', 'Linewidth', 2);
%     plot(shift_x(MRF_x_vec_B(1), origin_start_m, MRF_params.int_vel, time), 100*Tnew(MRF_params.nx_gamma_l), 'ro');


    hXLabel = xlabel('Position [mm]');
    hYLabel = ylabel('Eu Mole %');
    hTitle = title(['Concentration evolution after ',num2str(time/60, '%3.1f'),' minutes']);
    
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 15          );

    set( hTitle                    , ...
    'FontSize'   , 17          , ...
    'FontWeight' , 'bold'      );

    ylim([1 5]);
    xlim([7 27]);
    drawnow;
    
    frame = getframe(gcf);
    writeVideo(V,frame);
    
    % Update temperature and time
    MRF_Toft{n, 1}  = Tnew(1:length(MRF_x_vec_A));
    MRF_Toft{n, 2}  = Tnew(end-length(MRF_x_vec_B)+1:end);
    T           =   Tnew;
    time        =   time+dt;
    MRF_t_vec (n) = time;
    
%     keyboard
end

close(V);

% Multiplot
figure
for k = 1:length(Eu_data_old.t_vec)
%     subplot(fix((k-1)/4)+1, rem(k-1,4)+1, k);
    subplot(4, 3, k);
    current_time_min = Eu_data_old.t_vec(k);
    
    
    plot(Eu_data_old.x_vec, Eu_1D_new(:, k), 'b', 'Linewidth', 2); hold on;    
    if k < 5
        [found_t found_t_i] = min((SRF_t_vec./60 - current_time_min).^2);
        plot(SRF_x_vec_A*1000, 100*SRF_Toft{found_t_i, 1}, 'r', 'Linewidth', 2); hold on;
        plot(SRF_x_vec_B*1000, 100*SRF_Toft{found_t_i, 2}, 'r', 'Linewidth', 2);
    
    else
        mig_time = current_time_min*60 - SRF_params.t_end;
        [found_t found_t_i] = min((MRF_t_vec./60 - current_time_min).^2);
        plot(1000*shift_x(MRF_x_vec_A, origin_start_m, MRF_params.int_vel, mig_time), 100*MRF_Toft{found_t_i, 1}, 'r', 'Linewidth', 2); hold on;
        plot(1000*shift_x(MRF_x_vec_B, origin_start_m, MRF_params.int_vel, mig_time), 100*MRF_Toft{found_t_i, 2}, 'r', 'Linewidth', 2);
    end
    
    xlim([7 27])
    ylim([1 5])
    hTitle = title(['t = ' num2str(current_time_min) ' min']);
    hXLabel = xlabel('Position [mm]');
    hYLabel = ylabel('Eu mole%');
%     hLegend = legend('Exp', 'Solid', 'Liquid')
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 13          );

    set( hTitle                    , ...
    'FontSize'   , 15          , ...
    'FontWeight' , 'bold'      );
%     C_tot_vec(k) = max(cumtrapz(x_vec, Eu_t_mat_raw(:, k)));
    
%     xlim
end

% Assemble main finite difference matrix
function A = Assemble_A(x_vec_A, x_vec_B, params)
    n_rows = length(x_vec_A) + length(x_vec_B) + 2;
    n_cols = n_rows;

    D2 = sparse(n_rows, n_cols);
    
    s2s = params.s2s;
    s2l = params.s2l;
    s1 = params.s1;
    
    upper_diag = [-s2s*ones(1, length(x_vec_A) + 1), 0, 0, -s2l*ones(1, length(x_vec_B))]';
    diag = [(1+2*s2s)*ones(1, length(x_vec_A)), 0, 0, (1+2*s2l)*ones(1, length(x_vec_B)+1)]';
    lower_diag = [-s2s*ones(1, length(x_vec_A) -1), 0, 0,-s2l*ones(1, length(x_vec_B) + 1), 0]';
    D2 = spdiags(diag, 0, D2);
    D2 = spdiags(upper_diag, 1, D2);
    D2 = spdiags(lower_diag, -1, D2);
    
    D1 = sparse(n_rows, n_rows);
    upper_diag = [-s1*ones(1, length(x_vec_A) + 1), 0, 0, -s1*ones(1, length(x_vec_B))]';
    lower_diag = [s1*ones(1, length(x_vec_A) -1), 0, 0, s1*ones(1, length(x_vec_B) + 1), 0]';
    
    D1 = spdiags(upper_diag, 1, D1);
    D1 = spdiags(lower_diag, -1, D1);
    
    A = D2 + D1;
end

% Real x-position given MRF interface velocity and time
function x_out = shift_x(x_in, x0, int_vel, t)
    x_out = x_in + x0 + int_vel*t;
end

function struct_out = appendstruct(structA, structB)
% Appends structB to structA
    f = fieldnames(structB);
    struct_out = structA;
    for i = 1:length(f)
        struct_out = setfield(struct_out, f{i}, getfield(structB, f{i}));
    end
end

% Modifies main matrix to input initial conditions
function [A, b] = inputIC(A, b, T, params)
    K = params.K;

    A(params.nx_gamma_sr, :) = 0;
    A(params.nx_gamma_sr, params.nx_gamma_s - 1) = -params.Ds/2/params.dx;
    A(params.nx_gamma_sr, params.nx_gamma_sr) = params.Ds/2/params.dx;
    A(params.nx_gamma_sr, params.nx_gamma_ll) = params.Dl/2/params.dx;
    A(params.nx_gamma_sr, params.nx_gamma_l + 1) = -params.Dl/2/params.dx;
    A(params.nx_gamma_sr, params.nx_gamma_s) = params.int_vel;
    A(params.nx_gamma_sr, params.nx_gamma_l) = -params.int_vel;
    b(params.nx_gamma_sr) = 0;
    
    
    A(params.nx_gamma_ll, :) = 0;
    A(params.nx_gamma_ll, params.nx_gamma_s) = 1;
    A(params.nx_gamma_ll, params.nx_gamma_l) = -K;
    b(params.nx_gamma_ll) = 0;

end

% Modifies main matrix to input boundary conditions
function [A, b] = inputBC(BC, A, b, T, params)
    if BC.type == 'D'
        if BC.loc == 'L'
            A(1, :) = 0;
            A(1, 1) = 1;
            b(1) = BC.val;
        elseif BC.loc == 'R'
            A(end, :) = 0;
            A(end, end) = 1;
            b(end) = BC.val;
        end
    elseif BC.type == 'N'
        if BC.loc == 'L'
            A(1, :) = 0;
            A(1, 1) = 1 + 2*params.s2s;
            A(1, 2) = -2*params.s2s;
            
            b(1) = T(1) - 2*params.s2s*params.dx*BC.val/params.Dl;
        elseif BC.loc == 'R'
            A(end, :) = 0;
            A(end, end-1) = -2*params.s2l;
            A(end, end) = 1 + 2*params.s2l;
            
            b(end) = T(end) + 2*params.s2l*params.dx*BC.val/params.Dl;
        end
    end
end