% This script solves a 1D segregation problem in a diffusion couple with a solid and liquid phase
% using Finite Fourier transforms.  The solution is thus a summation of products of a
% temporal component exponential and spatial basis function.
% The initial conditions is given by Anton's neutron imaging data, which
% has been reduced to a 1D concentration profile using the inverse-radon technique.

clc; clear all;close all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
delete 'EuSegregationDiary.txt'
diary 'EuSegregationDiary.txt'

% Load Eu data
Eu_data = load('.\Eu_data.mat');
toload = 0;
gamma = 0.5;

% Physical parameters
K = .45; % Partition coefficient
D1 =  1.9e-10; % m^2/s  Solid diffusivity
D2 = 2.5e-10; % m^2/s   Liquid diffusivity
D = D2/D1;
Ctot = Eu_data.C0_tot;

% Time parameters
t0 = 0;
dt = 1320; % 22 min
t_final = 4000*dt;

% Max number of Fourier terms
n_max = 100;
RC_orthogonality = false;

% load past data if 1
% solve new if 0
if toload == 1
    load('.\Eigen_results.mat');
    fprintf('No need for symbolic math!  Loaded previous results!\n');

else
    
    % Eigen search parameters
    lambda_max_range = 6.2778*n_max*1.1;
    y_guess_lim = 10;
    n_sample_pts = 100000*n_max;

    % Characteristic scales
    char_L = max(Eu_data.x_full)*1e-3;
    char_t = char_L.^2/D1;
    char_C = 100*Ctot;

    % Nondimensionalization overwrite
    D1 = 1;
    D2 = D;

    % Plot the characteristic equation (to set initial guess)
    lambda_vec = cumtrapz(fliplr(diff(logspace(0, log10(lambda_max_range), n_sample_pts))))';

    t_vec = logspace(log10(1), log10(t_final), 24)./char_t;

    t_vec(1:7)=[]

    L1 = sqrt(D1);
    L2 = sqrt(D2);

    tot_dom_len = max(Eu_data.x_full)*1e-3/char_L;
    gamma = gamma*tot_dom_len;

    x_res = 8000;
    x_vec_A = linspace(0, gamma, round(x_res*gamma));
    x_vec_B = linspace(gamma, tot_dom_len, round(x_res*(1-gamma)));

    Eu0_interp_A = interp1(Eu_data.x_full./char_L/1e3, Eu_data.Eu0_full/100, x_vec_A);
    Eu0_interp_B = interp1(Eu_data.x_full./char_L/1e3, Eu_data.Eu0_full/100, x_vec_B);

    fprintf('D1 = %f\n', D1);
    fprintf('D2 = %f\n', D2);
    fprintf('gamma = %f\n', gamma);
    fprintf('n_max = %f\n', n_max);

    syms an An x_s lam_s gamma_s bn Bn D_s L1s L2s K_s

    % Symbolic f( ) and g( ) without gamma (only x_s)
    f_s =   cos(lam_s*x_s/L1s);
    g_s =   cos(lam_s.*x_s/L2s) + tan(lam_s./L2s).*sin(lam_s.*x_s/L2s);
    fp_s =  -sin(lam_s*x_s/L1s); % my version 
%     fp_s =  sin(lam_s*x_s/L1s); % JHP's
    gp_s =  -sin(lam_s./L2s*x_s) + tan(lam_s/L2s).*cos(lam_s./L2s*x_s); % my version 
%     gp_s =  sin(lam_s./L2s*x_s) - tan(lam_s/L2s).*cos(lam_s./L2s*x_s); % JHP's

    % Symbolic f( ) and g( ) at x_s = gamma
    f_sg =   subs(f_s, x_s, gamma_s);
    g_sg =   subs(g_s, x_s, gamma_s);
    fp_sg =  subs(fp_s, x_s, gamma_s);
    gp_sg =  subs(gp_s, x_s, gamma_s);

    % Function-handle of f( ) and g( ) at x_s = gamma
    f = matlabFunction(f_sg);
    g = matlabFunction(g_sg);
    fp = matlabFunction(fp_sg);
    gp = matlabFunction(gp_sg);
    
    intf2_s = int(f_s^2, x_s, 0, gamma_s);
    intg2_s = int(g_s^2, x_s, gamma_s, 1);

    fprintf('Eigenvalue pre-search\n')

    eigen_prob_fh = @(lam) (L1.*lam.*K.*fp(L1, gamma, lam).*g(L2, gamma, lam) - L2.*lam.*f(L1, gamma, lam).*gp(L2, gamma, lam));
%     eigen_prob_NA = @(lam) K*L1*cos(lam/L2*(1-gamma)).*sin(lam*gamma/L1) + L2*sin(lam/D2*(1-gamma)).*cos(lam*gamma/L1);

    eigen_prob = eigen_prob_fh(lambda_vec);
%     eigen_prob_NA = eigen_prob_NA(lambda_vec);

    % Filter function evaluations of the characteristic equation to those
    % between +/-y_guess_lim (near the root)
    guess_A = eigen_prob;
    guess_A(guess_A < -y_guess_lim) = 0;
    guess_A(guess_A > y_guess_lim) = 0;

    figure    
    plot(lambda_vec, eigen_prob, '.'); hold on;
    ylim([-30 30]);
    hXLabel = xlabel('Eigenvalue \lambda');
    hYLabel = ylabel('Characteristic equation');
    
    drawnow;

    % Counter loop to 
    k = 2;
    n_guesses = zeros(1);

    for i = 1:length(guess_A)-1
        if guess_A(i)*guess_A(i+1)<0
            n_guesses(k) = lambda_vec(i);
            k = k+1;
        end 
    end
    

    if length(n_guesses) < n_max
        fprintf('Insufficient pre-search.  Expand pre-search domain.  Exiting.');
        return;
    end

    fprintf('Calculating normalization constants symbolically...')

    % % Normalization equation
    eq1 = bn^2*int(f_s^2, x_s, 0, gamma_s) + K_s*Bn^2*int(g_s^2, x_s, gamma_s, 1) - 1;
    % Concentration (dis)continuity
    eq2 = bn*f_sg - K_s*Bn*g_sg;
    % Flux continuity
    eq3 = D1*bn*lam_s/L1*fp_sg - D2*Bn*lam_s/L2*gp_sg;
    % Total Mass conservation
    eq4 = bn*int(f_s, x_s, 0, gamma_s) + Bn*int(g_s, x_s, gamma_s, 1) - Ctot;

    eq5 = -D1*bn*lam_s/L1*fp_sg; % Check the sign this flux

    fprintf('done!\n')

    Bn_swan_s = (K_s^2*g_sg^2/f_sg^2*intf2_s + K_s*intg2_s)^(-1/2);
    Bn_swan_f = matlabFunction(Bn_swan_s, 'vars', [L1s,L2s,gamma_s,lam_s, K_s]);
    bn_swan_f = matlabFunction(K*Bn_swan_s*g_sg/f_sg, 'vars', [L1s,L2s,gamma_s,lam_s, K_s]);
    
    phi_nAs = bn*f_s;
    phi_nBs = Bn*g_s;

    phi_nAf = matlabFunction(phi_nAs, 'vars', [L1s, bn, lam_s,x_s]);
    phi_nBf = matlabFunction(phi_nBs, 'vars', [L2s, Bn, lam_s,x_s]);

    % Use the original symbolic equations as functions for reality checks
    % later
    RC1_f = matlabFunction(eq1, 'vars', [L1s,L2s,bn,Bn,gamma_s,lam_s, K_s]);
    RC2_f = matlabFunction(eq2, 'vars', [L1s,L2s,bn,Bn,gamma_s,lam_s, K_s]);
    RC3_f = matlabFunction(eq3, 'vars', [L1s,L2s,bn,Bn,gamma_s,lam_s, K_s]);
    RC4_f = matlabFunction(eq4, 'vars', [L1s,L2s,bn,Bn,gamma_s,lam_s, K_s]);
    RC5_f = matlabFunction(eq5, 'vars', [L1s,L2s,bn,Bn,gamma_s,lam_s, K_s]);


    % @(Bn, bn, lam) lam*D1*bn/L1*fp(lam) - lam*Bn*D2/K/L2*gp(lam);
    fprintf('done!\n')

    fprintf('Calculating weighting functions symbolically...')
%      f(x) = p1*x^3 + p2*x^2 + p3*x + p4
% Coefficients (with 95% confidence bounds):
%        p1 =    -0.05133  (-0.05623, -0.04643)
%        p2 =      0.1877  (0.1793, 0.1961)
%        p3 =     -0.2485  (-0.253, -0.2441)
%        p4 =      0.1422  (0.1414, 0.1429)
%     Eu0_s = (-.05133*x_s^3 + 0.1877*x_s^2 -0.2485*x_s + .1422);
    Eu0_s = eps*x_s + .05;
    theta_n0_s = int(phi_nAs*Eu0_s, x_s, 0, gamma) + K_s*int(phi_nBs*Eu0_s, x_s, gamma, 1);
    theta_n0_f = matlabFunction(theta_n0_s, 'vars',[L1s,L2s,bn,Bn,lam_s, K_s]);

    fprintf('done!\n')

    bn_vec = zeros(1, n_max);
    An_vec = zeros(1, n_max);
    Bn_vec = zeros(1, n_max);

    % Numerically solve for nth lambda
    lambdas = zeros(1, n_max);
    fsolve_options = optimset('Display','off');
    attempt = 1;
    n= 2;
    hold all;


    while n<=n_max

        if attempt ==1
            guess = n_guesses(n);
        end

        fprintf('n = %d, attempt %d, guess = %f\n', n, attempt, guess);
        [lambdas(n), check] = fsolve(eigen_prob_fh, guess, fsolve_options);

        fprintf('\tlambda_%d = %f\n', n, lambdas(n));

    %     check = abs(eigen_prob1(lambdas(n)));
    %     fprintf('\tCheck = %e', abs(check));

        if abs(check)<1e-3
            plot(lambdas(n), 0, 'r*');
            guess = lambdas(n) + lambdas(1);
            attempt = 1;

            fprintf('\t\tCompute normalization constants:\n')

            Bn_do = Bn_swan_f(L1,L2,gamma,lambdas(n), K);
            bn_do = bn_swan_f(L1,L2,gamma,lambdas(n), K);

            [bn_d Bn_d score_bnBn] = RCsuite2(L1, L2, bn_do, Bn_do, gamma, lambdas(n), K, RC1_f, RC2_f, RC3_f, RC4_f, RC5_f);

            bn_vec(n) = bn_do;
            Bn_vec(n) = Bn_do;
            An_vec(n) = Bn_do*tan(lambdas(n)/L2);

            % Reality check for orthogonality
            if RC_orthogonality
                if n > 2
                    RC_ortho_f = matlabFunction(bn_vec(n)*bn_vec(n-1)*int(subs(f_s, lam_s, lambdas(n))*subs(f_s, lam_s, lambdas(n-1)), x_s, 0, gamma_s) + ...
                        K*Bn_vec(n)*Bn_vec(n-1)*int(subs(g_s, lam_s, lambdas(n))*subs(g_s, lam_s, lambdas(n-1)), x_s, gamma_s, 1));

                    RC_ortho = RC_ortho_f(L1,L2,gamma);
                    fprintf('\t\tRC_ortho = %f\n', RC_ortho);
                    if RC_ortho > 1e-3
                        fprintf('\t\t!!!Current basis function is NOT orthogonal with the previous!!!\n')
                    end
                end
            end

            n = n+1;

            continue
        else
            guess = lambdas(n) + lambdas(1);
            attempt = attempt + 1;
            fprintf('\n')
        end


    end

    xlim([0, 1.1*max(lambdas)]);
    plot(0, 0, 'r*')
    set([hXLabel, hYLabel]  , ...
    'FontSize'   , 15          );
    fprintf('All eigenvalues found!\n')

    save('C:\Users\Kerry\Dropbox\UMN\Research\Eu Segregation\Eigen_results.mat', 'Eu0_s', 'bn_vec', 'Bn_vec', 'lambdas', 'x_vec_A', 'x_vec_B', 'hXLabel', 'hYLabel' ...
        , 'Eu0_interp_A', 'Eu0_interp_B');
end


figure1 = figure;hold all;
Eu0_f = matlabFunction(Eu0_s);
Eu0_v_A = Eu0_f(x_vec_A);
Eu0_v_B = Eu0_f(x_vec_B);

plot(x_vec_A, Eu0_v_A, 'g', 'Linewidth', 2, 'HandleVisibility','off');
plot(x_vec_B, Eu0_v_B, 'g', 'Linewidth', 2);

theta_n0_tv = zeros(length(t_vec), n_max); % this is the inner product term
theta_n_tv = zeros(length(t_vec), n_max); % this is the exponential decay term
phi_nAv = zeros(n_max, length(x_vec_A)); % 2D array of eigenfunctions n_max rows, x_res cols
phi_nBv = zeros(n_max, length(x_vec_B));

%% Assembly 
for n = 1:n_max
    
    % Handle lambda = 0 case separately
    if lambdas(n) == 0
        phi0_B = (sqrt(K^2*gamma + K*(1-gamma)))^(-1);
        phi0_A = K*phi0_B;
        
        phi_nAv(n, :) = ones(size(phi_nAv(n, :)))*phi0_A;
        phi_nBv(n, :) = ones(size(phi_nBv(n, :)))*phi0_B;
        
        theta_n0_v(n) = double(phi0_A*int(Eu0_s, x_s, 0, gamma) + K*phi0_B*int(Eu0_s, x_s, gamma, 1));
    else
        phi_nAv(n, :) = phi_nAf(L1, bn_vec(n), lambdas(n), x_vec_A);
        phi_nBv(n, :) = phi_nBf(L2, Bn_vec(n), lambdas(n), x_vec_B);
        
        theta_n0_v(n) = theta_n0_f(L1, L2, bn_vec(n), Bn_vec(n), lambdas(n), K);
        
        fprintf('\tK reality check: %f \n', phi_nAv(n, end)/phi_nBv(n, 1))
    end    
end

Q_n_Av = theta_n0_v'.*phi_nAv;
Q_n_Bv = theta_n0_v'.*phi_nBv;
Ex = exp(-t_vec'*lambdas.^2); % Rows of time, Cols of lambda

% The final concentration profiles...temporal part multiplied by spatial
% part
C_A = Ex(:, 1:n_max)*Q_n_Av(1:n_max, :);
C_B = Ex(:, 1:n_max)*Q_n_Bv(1:n_max, :);

% Plotting
c = cool(length(t_vec)); %colormap
legend_t_vec = 0;

for k = 1:length(t_vec)
    if rem(k, 3)==1
%     t = t_vec(k);
        plot(x_vec_A, C_A(k, :), 'Linewidth', 2, 'Color', c(k, :), 'HandleVisibility','off');
        plot(x_vec_B, C_B(k, :), 'Linewidth', 2, 'Color', c(k, :));
        
        fprintf('t = %i\n', t_vec(k))
        fprintf('\tK reality check: %f \n', C_A(k, end)/C_B(k, 1));
        
        if k > 1
            legend_t_vec = [legend_t_vec t_vec(k)*char_t/60];
        end
    end
end

legend_str = cell(1,length(legend_t_vec)+1);
legend_str{1} = 'Initial condition';
for k = 1:length(legend_t_vec)
    legend_str{k+1} = ['t = ' num2str(legend_t_vec(k), '%3.1e') ' min'];
end
legend(legend_str, 'Location', 'Southeast');
hXLabel1 = xlabel('Position');
hYLabel1 = ylabel('Eu mole fraction');
set([hXLabel1, hYLabel1]  , ...
'FontSize'   , 15          );
diary off;

% Suite of reality checks to ensure equations are solved correctly
function [bn_d Bn_d score] = RCsuite2(L1, L2, bn_d, Bn_d, gamma, lambda, K, RC1_f, RC2_f, RC3_f, RC4_f, RC5_f)
    fprintf('\t\tTesting bn = %f,\t Bn = %f\n', bn_d, Bn_d);
        
    RC1 = RC1_f(L1, L2, bn_d, Bn_d, gamma, lambda, K);
    fprintf('\t\t\tRC1, Normalization (should be 0):   %f\n', RC1);     

    RC2 = RC2_f(L1, L2, bn_d, Bn_d, gamma, lambda, K);
    fprintf('\t\t\tRC2, C discontinuity (should be 0): %f\n', RC2);

    RC3 = RC3_f(L1, L2, bn_d, Bn_d, gamma, lambda, K);
    fprintf('\t\t\tRC3, Flux continuity (should be 0): %f\n', RC3);
    
%     RC4 = RC4_f(L1, L2, bn_d, Bn_d, gamma, lambda);
%     fprintf('\t\t\tRC3, Ma conservation (should be 0): %f\n', RC4);
    
%     if(RC5_f(L1, L2, bn_d, Bn_d, gamma, lambda, K) > 0)
%         RC5 = 0;
%         fprintf('\t\t\tRC5, Flux direction (should be 0): %f\n', RC5);
%     else
%         RC5 = 1;
%         bn_d = -bn_d;
%         Bn_d - Bn_d;
%         fprintf('\t\t\tRC5, Flux direction (should be 0): %f (bn/Bn sign flipped!)\n', RC5);
%         
%     end
    
    score = (RC1)^2 + (RC2)^2 + (RC3)^2;
    
    fprintf('\t\t\tbn/RC score: %f\n', score);
end
