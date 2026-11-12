% run_kuramoto_model.m
%
% Kuramoto model simulation of interhemispheric neural entrainment.
%
% This is the computational model used in Okazaki et al., "Asymmetric neural
% entrainment at resonance frequencies underlies unilateral spatial neglect".
% Running this script once executes all simulation conditions reported in
% Figure 9 of the paper and computes, for every trial and every stimulation
% frequency, the order parameter of the left and the right hemisphere.
%
% The random seed is not fixed, as in the original analysis. Natural
% frequencies and initial phases are drawn anew on every run, so the
% numerical values vary slightly from run to run.
%
% Conditions:
%   interhemispheric  : asymmetric interhemispheric coupling (K_inter_RL > K_inter_LR)
%   nostim            : same parameters as 'interhemispheric' but without stimulation (b = 0)
%   intrahemispheric  : asymmetric intrahemispheric coupling (K_intra_L > K_intra_R)
%   natural_frequency : different natural frequencies between hemispheres (omegaL > omegaR)
%
% Output (left in the workspace):
%   results.<condition>.rho_L         [num_trials x num_stimulation_frequencies]
%   results.<condition>.rho_R         [num_trials x num_stimulation_frequencies]
%   results.<condition>.omega_f_list  stimulation frequencies (Hz)
%   results.<condition>.params        parameters used for that condition
%
% Each entry of rho_L / rho_R is the order parameter of one hemisphere,
% averaged over the stimulation period of a single trial.
%
% Requirements: MATLAB (base product only; no additional toolboxes).

%% Common parameters (identical for all conditions)
omega_f_list = 3:3:30;  % stimulation frequencies (Hz)
sigma = 2;              % SD of the natural frequency distribution, rad/s (approximately 0.32 Hz)

T_total  = 10;          % simulation duration (s)
stimtime = [0 10];      % stimulation period (s): 10 s of continuous stimulation

num_trials = 21;        % number of trials per stimulation frequency

num_oscillators = 100;              % total number of oscillators
num_L = num_oscillators / 2;        % number of oscillators in the left hemisphere
num_R = num_oscillators / 2;        % number of oscillators in the right hemisphere

%% Condition-specific parameters
% omegaL, omegaR       : natural frequencies of the left / right hemisphere (Hz)
% K_intra_L, K_intra_R : coupling strength within the left / right hemisphere
% K_inter_LR           : coupling strength from the left to the right hemisphere
% K_inter_RL           : coupling strength from the right to the left hemisphere
% b                    : strength of the external stimulation
params.interhemispheric = struct('omegaL', 9, 'omegaR', 9, ...
    'K_intra_L', 0.2, 'K_intra_R', 0.2, 'K_inter_LR', 0.2, 'K_inter_RL', 0.8, 'b', 2);

params.nostim = struct('omegaL', 9, 'omegaR', 9, ...
    'K_intra_L', 0.2, 'K_intra_R', 0.2, 'K_inter_LR', 0.2, 'K_inter_RL', 0.8, 'b', 0);

params.intrahemispheric = struct('omegaL', 9, 'omegaR', 9, ...
    'K_intra_L', 0.8, 'K_intra_R', 0.2, 'K_inter_LR', 0.2, 'K_inter_RL', 0.2, 'b', 2);

params.natural_frequency = struct('omegaL', 9, 'omegaR', 7, ...
    'K_intra_L', 0.2, 'K_intra_R', 0.2, 'K_inter_LR', 0.2, 'K_inter_RL', 0.2, 'b', 2);

%% Run all conditions
condition_names = {'interhemispheric', 'nostim', 'intrahemispheric', 'natural_frequency'};

results = struct();
for ci = 1:length(condition_names)
    condition = condition_names{ci};
    p = params.(condition);

    fprintf('Running condition %d/%d: %s\n', ci, length(condition_names), condition);

    [rho_L, rho_R] = run_condition(p, omega_f_list, sigma, T_total, stimtime, ...
        num_trials, num_L, num_R);

    results.(condition).rho_L        = rho_L;
    results.(condition).rho_R        = rho_R;
    results.(condition).omega_f_list = omega_f_list;
    results.(condition).params       = p;
end

%% Local functions
function [rho_all_L, rho_all_R] = run_condition(p, omega_f_list, sigma, T_total, ...
    stimtime, num_trials, num_L, num_R)
% Simulate all trials and all stimulation frequencies of one condition.

% Order parameters
rho_all_L = zeros(num_trials, length(omega_f_list));
rho_all_R = zeros(num_trials, length(omega_f_list));

for idx = 1:length(omega_f_list)
    omega_f = omega_f_list(idx); % current stimulation frequency

    for sim = 1:num_trials
        % Natural frequencies and initial phases
        omega_array_L = 2 * pi * p.omegaL + sigma * randn(num_L, 1);
        omega_array_R = 2 * pi * p.omegaR + sigma * randn(num_R, 1);
        phi_initial_L = 2 * pi * rand(num_L, 1) - pi;
        phi_initial_R = 2 * pi * rand(num_R, 1) - pi;

        % Combine the initial values of the state variables
        phi_initial = [phi_initial_L; phi_initial_R];

        % Solve the ODE
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        [t_out, phi_out] = ode45(@(t, phi) ode_function(t, phi, omega_array_L, omega_array_R, omega_f, p.b, p.K_intra_L, p.K_intra_R, p.K_inter_LR, p.K_inter_RL, num_L, num_R, stimtime), [0 T_total], phi_initial, options);

        % Extract the data during the stimulation period
        stim_indices = find(t_out >= stimtime(1) & t_out <= stimtime(2));
        phi_L_stim = phi_out(stim_indices, 1:num_L);
        phi_R_stim = phi_out(stim_indices, num_L+1:end);

        % Order parameter during the stimulation period
        zL = mean(exp(1i * phi_L_stim), 2);
        zR = mean(exp(1i * phi_R_stim), 2);
        rho_L = abs(zL);
        rho_R = abs(zR);

        % Store the mean values
        rho_all_L(sim, idx) = mean(rho_L);
        rho_all_R(sim, idx) = mean(rho_R);
    end
end
end

function dphi_dt = ode_function(t, phi, omega_array_L, omega_array_R, omega_f, b, K_intra_L, K_intra_R, K_inter_LR, K_inter_RL, num_L, num_R, stimtime)
    % Split the state variables
    phi_L = phi(1:num_L);
    phi_R = phi(num_L+1:end);

    % Apply the external stimulation (only during the stimulation period)
    if t >= stimtime(1) && t <= stimtime(2)
        external_force_L = b * sin(2*pi*omega_f * t - phi_L);
        external_force_R = b * sin(2*pi*omega_f * t - phi_R);
    else
        external_force_L = zeros(num_L, 1);
        external_force_R = zeros(num_R, 1);
    end

    % Intrahemispheric and interhemispheric coupling terms
    % Left hemisphere
    sum_left_L = sum(sin(phi_L' - phi_L), 2) / num_L;
    sum_right_L = sum(sin(phi_R' - phi_L), 2) / num_R;

    dphi_L = omega_array_L + K_intra_L * sum_left_L + K_inter_RL * sum_right_L + external_force_L;

    % Right hemisphere
    sum_right_R = sum(sin(phi_R' - phi_R), 2) / num_R;
    sum_left_R = sum(sin(phi_L' - phi_R), 2) / num_L;

    dphi_R = omega_array_R + K_intra_R * sum_right_R + K_inter_LR * sum_left_R + external_force_R;

    % Combine the derivatives
    dphi_dt = [dphi_L; dphi_R];
end
