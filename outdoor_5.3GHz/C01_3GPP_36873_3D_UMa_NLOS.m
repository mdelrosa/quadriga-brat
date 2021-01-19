%% Multi-frequency simulations

% This tutorial demonstrates how to perform simultaneous multi-frequency simulations at two carrier
% frequencies: 2.6 GHz and 28 GHz in an Urban-Macrocell deployment. The BS is equipped with two
% different array antennas. A conventional high-gain antenna operates at 2.6 GHz. The higher
% frequency band uses a massive-MIMO array antenna with in an 8x8 dual-polarized setup. The model is
% consistent in both, the spatial domain and the frequency domain. Simulation assumptions are in
% accordance with 3GPP 38.901 v14.1.0 (see Section 7.6.5 Correlation modeling for multi-frequency
% simulations).
%
% Identical parameters for each frequency:
%
% * LOS / NLOS state must be the same
% * BS and MT positions are the same (antenna element positions are different!)
% * Cluster delays and angles for each multi-path component are the same
% * Spatial consistency of the LSPs is identical
%
% Differences:
%
% * Antenna patterns are different for each frequency
% * Path-loss is different for each frequency
% * Path-powers are different for each frequency
% * Delay- and angular spreads are different
% * K-Factor is different
% * XPR of the NLOS components is different


%% Basic setup
% Multiple frequencies are set in the simulation parameters by providing a vector of frequency
% sample points. A new layout is created with one 25 m high BS positions and 100 MT positions. The
% MTs are placed in accordance with the 3GPP assumptions, where 80% of them are situated indoors at
% different floor levels.

% define batch_num/N_samples
debug_flag = 1;
if debug_flag==0
	N_samples = 2500;
else
	N_samples = 10;
end

try 
	if ~batch_num
		batch_num = 1;
	end
catch
	fprintf('Batch number undefined. Setting batch_num=1\n');
	batch_num = 1;
end

% seed random number generator; prevent identical datasets between different batches
rng(batch_num*10);

% for mobile MTs, we define the following:
sample_density = 1.2;   % [samples / half wavelength]
% track_distance = 0.25;  % [meters] - for 5.2GHz
track_distance = 4;  % [meters] - for 300MHz
track_speed = 0.9;  % [meters/second] - for 300MHz
area_len = 400;         % area of MT initial positions [meters]
area_half = area_len/2;
timeslots = 10;
best_shift = 84;

N = 1;
M = 32;
no_sc = 1024;
sc_bw = 2e4;
n_truncate = 96;

Hur_down = zeros(N_samples,timeslots,M,n_truncate);
Hur_up = zeros(N_samples,timeslots,M,n_truncate);
P_excess_down = zeros(N_samples,timeslots);
P_excess_up   = zeros(N_samples,timeslots);

% k -- shifting sweep
k_lim = 256;
k_list = (-k_lim:k_lim);
pow_ratio = zeros(1,length(k_list));
batch_times = zeros(1,N_samples);

fprintf(sprintf('Batch #%d -- %d samples. \n',batch_num,N_samples));

for i_sample = 1:N_samples
    
    strout=sprintf("Sample %03d / %03d -> ", i_sample, N_samples);
    fprintf(strout);
    tic
    
    set(0,'defaultTextFontSize', 18)                        % Default Font Size
    set(0,'defaultAxesFontSize', 18)                        % Default Font Size
    set(0,'defaultAxesFontName','Times')                    % Default Font Type
    set(0,'defaultTextFontName','Times')                    % Default Font Type
    set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
    set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type
    set(0,'DefaultFigurePaperSize',[14.5 7.7])            	% Default Paper Size
    
    s = qd_simulation_parameters;
    % s.center_frequency = [5.1e9 5.1e9 5.3e9];              % Assign two frequencies, the first one is redundent frequency
    s.center_frequency = [2.6e8 2.6e8 3.0e8];              % Assign two frequencies, the first one is redundent frequency
    freq_str = '300MHz';
    s.sample_density = sample_density;                     % from t_09_speed_profile_interpolation.m - samples per half-wavelength
    s.show_progress_bars = 0;                               % Disable progress bars
    no_user = 1;
    
    % from t09_speed_profile_interpolation.m - define a linear track for MT
    % to travel along based on random init position/direction 
    
    % randomize init point for MT in square area defined by area_len
    x_i = area_len*rand(1) - area_half;
    y_i = area_len*rand(1) - area_half;
    theta = pi*(2*rand(1) - 1);
    t = qd_track('linear', track_distance, theta);            % 20 m track, direction SE
    t.initial_position = [x_i; y_i; 1.5];                        % Start position
    t.interpolate_positions( s.samples_per_meter );         % Apply sample density
    t.set_speed( track_speed );         % Apply sample density
    % t.interpolate_positions( 128/20 );                      % Interpolate
    t.segment_index       = [1];                      % Assign segments
    % t.scenario            = {'BERLIN_UMa_LOS','BERLIN_UMa_NLOS','BERLIN_UMa_LOS'};
    t.scenario            = {'3GPP_3D_UMa_LOS'};
    t.name = 'Rx1';
    
    % from t09_speed_profile_interpolation.m - apply linear track to link
    % layout handle
    l = qd_layout( s );                                     % New QuaDRiGa layout
    % l.tx_array = qd_arrayant('dipole');                     % Set Dipole antenna
    % l.rx_array = qd_arrayant('dipole');                     % Set Dipole antenna
    % l.tx_position(3) = 25;                                  % BE height
    l.rx_track = t;                                         % Assign track
    
    % l = qd_layout( s );                                     % New QuaDRiGa layout
    l.tx_position = [0 0 20]';              % 20 m BS height
    % l.no_rx = no_user;                                      % 200 MTs
    
    % l.randomize_rx_positions( 100, 1.5, 1.5, 0 );           % Assign random user positions
    % l.rx_position(1,:) = l.rx_position(1,:) + 110;                %
    
    % floor = randi(5,1,l.no_rx) + 3;                         % Set random floor levels
    % for n = 1:l.no_rx
    %     floor( n ) =  randi(  floor( n ) );
    % end
    % l.rx_position(3,:) = 3*(floor-1) + 1.5;
    
    % indoor_rx = l.set_scenario('3GPP_3D_UMa',[],[],0);    % Set the scenario
    % l.rx_position(3,~indoor_rx) = 1.5;                      % Set outdoor-users to 1.5 m height
    % l.visualize
    
    %% Antenna set-up
    % Two different antenna configurations are used at the BS. The 2.6 GHz antenna is constructed from 8
    % vertically stacked patch elements with +/- 45 degree polarization. The electric downtilt is set to
    % 8 degree. The mm-wave antenna uses 64 dual-polarized elements in a 8x8 massive-MIMO array
    % configuration. The antennas are assigned to the BS by an array of "qd_arrayant" objects. Rows
    % correspond to the frequency, columns to the BS. There is only 1 BS in the layout. The mobile
    % terminal uses a vertically polarized omni-directional antenna for both frequencies.
    
    % a_35000_Mhz  = qd_arrayant( '3gpp-3d',  16, 8, s.center_frequency(1), 6, 8 );
    % a_35000_Mhz  = qd_arrayant.generate( '3gpp-mmw',  16, 8, s.center_frequency(1), 6, 8 );
    a_redundent  = qd_arrayant.generate( '3gpp-3d',  N, M, s.center_frequency(1), 1 ); % originally, last arg was 3
    a_uplink  = qd_arrayant.generate( '3gpp-3d',  N, M, s.center_frequency(2), 1 );    % yielded +/-45deg polarized
    a_downlink  = qd_arrayant.generate( '3gpp-3d',  N, M, s.center_frequency(3), 1 );  % elements, yielding 2x elements/antenna
    
    l.tx_array(1,1) = a_redundent;                           
    l.tx_array(2,1) = a_uplink;
    l.tx_array(3,1) = a_downlink;
    
    l.rx_array = qd_arrayant('omni');                       % Set omni-rx antenna
    
    
    %% Generate channel coefficients
    % Channel coefficients are generated by calling "l.get_channels". The output is an array of QuaDRiGa
    % channel objects. The first dimension corresponds to the MTs (100). The second dimension
    % corresponds to the number of BSs (1) and the third dimension corresponds to the number of
    % frequencies.

    c = l.get_channels;
    
    % get sample time
    t=toc;
    batch_times(i_sample) = t;
    avg_time = mean(batch_times(1:i_sample));
    proj_time = (N_samples - i_sample)*avg_time / 60;
    strout=sprintf("Sample Time: %03.1f sec - Sample Avg. Time: %03.1f sec - Projected Batch Time: %03.1f min \n", t, avg_time, proj_time);
    fprintf(strout);
    
    freq_response = zeros(2,timeslots,M,no_sc); % init freq response
    H_freq_down = zeros(timeslots,M,no_sc);
    H_freq_up = zeros(timeslots,M,no_sc);
    H_ang_down_sample = zeros(timeslots,M,no_sc);
    H_ang_up_sample = zeros(timeslots,M,no_sc);
    for t_i = 1:timeslots
        H_freq_down(t_i,:,:) = c(1,2).fr(no_sc*sc_bw,no_sc,t_i); %clear CFR; CFR(:) = abs(freq_response(1,:,3)); CFR = CFR(1:2:end);
        H_freq_up(t_i,:,:) = c(1,3).fr(no_sc*sc_bw,no_sc,t_i);
        H_ang_down_sample(t_i,:,:) = fft(squeeze(H_freq_down(t_i,:,:)), [], 2);
        H_ang_down_sample(t_i,:,:) = ifft(H_ang_down_sample(t_i,:,:), [], 1);
        H_ang_up_sample(t_i,:,:) = fft(squeeze(H_freq_up(t_i,:,:)), [], 2);
        H_ang_up_sample(t_i,:,:) = ifft(H_ang_up_sample(t_i,:,:), [], 1);
    end
    try
        H_ang_down_sample = circshift(H_ang_down_sample, best_shift, 3);
        H_ang_up_sample = circshift(H_ang_up_sample, best_shift, 3);
    catch
        for k = k_list
            a_shift = circshift(squeeze(H_ang_down_sample(1,:,:)),k,2); % for T = 1
            % a_shift = squeeze(circshift(a_at(timeslot,:,:),k,2)); % for T > 1
            temp = sum(a_shift(:,1:n_truncate).*conj(a_shift(:,1:n_truncate)), 'all');
            pow_ratio(k+k_lim+1) = pow_ratio(k+k_lim+1) + (temp / sum(a_shift.*conj(a_shift), 'all')) / N_samples;
        end
    end

    Hur_down(i_sample,:,:,:) = H_ang_down_sample(:,:,1:n_truncate);
    Hur_up(i_sample,:,:,:) = H_ang_up_sample(:,:,1:n_truncate);
    P_excess_down(i_sample,:) = sum(sum(conj(H_ang_down_sample(:,:,n_truncate+1:end)).*H_ang_down_sample(:,:,n_truncate+1:end), 2), 3);
    P_excess_up(i_sample,:)   = sum(sum(conj(H_ang_up_sample(:,:,n_truncate+1:end)).*H_ang_up_sample(:,:,n_truncate+1:end), 2), 3);
end

f_down = sprintf('H_ang_quadriga_%d_down.mat', batch_num);
f_up   = sprintf('H_ang_quadriga_%d_up.mat', batch_num);
P_down = sprintf('P_excess_quadriga_%d_down.mat', batch_num);
P_up   = sprintf('P_excess_quadriga_%d_up.mat', batch_num);
save(f_down,'Hur_down');
save(f_up,'Hur_up');
save(P_down,'P_excess_down');
save(P_up,'P_excess_up');
% clear all;

%% plot circ

if debug_flag
    try
        [m, idx] = max(pow_ratio);

        figure(1); clf; hold on;
        plot(k_list,pow_ratio);
        % title_str = sprintf('Circular Shift on Indoor 5.3GHz (max=%d)', idx - (k_lim+1));
        ylabel('P_{truncate} / P_{total}');
        title(sprintf('Quadriga Circular Shift on %s (max=%d, N=%d samples)', freq_str, idx-(k_lim+1), N_samples));
    catch
        fprintf("--- best_shift defined. Skipping pow_ratio plot. ---");
    end


    %% sanity checking - angular delay domain
    figure(2); clf; hold on;
    H_samp = squeeze(H_ang_down_sample);
    for t_i = 1:timeslots
        row = round(t_i / timeslots) + 1;
        col = round(t_i / timeslots) + 1;
        subplot(2,5,t_i);
        surf(10*log10(abs(squeeze(H_samp(t_i,:,1:n_truncate)))), 'EdgeColor', 'none');
        view(0,90);
        xlabel('delay');
        ylabel('angle');
        title(sprintf("t_{%d}", t_i));
    end
    % sgtitle("QuaDRiGa -- Outdoor 5.3GHz - 0.9m/s mobility - 40ms feedback interval");
    sgtitle("QuaDRiGa -- Outdoor 300MHz - 0.9m/s mobility - 40ms feedback interval");

    %% inspect power ratio for different truncation windows
%     truncate_lim = 128;
%     truncate_stride = 4;
%     l_list = (0:truncate_stride:truncate_lim);
%     pow_list = zeros(1,length(l_list));
%     
%     for l = l_list
%        a_drop = H_ang_down(l+1:end,:);
%        del, ang = size(a_drop);
%        pow_idx = round(l/truncate_stride)+1;
%        pow_list(pow_idx) = pow_list(pow_idx) + (sqrt(sum(a_drop .* conj(a_drop), 'all')) / N_samples);
%     end
%     
%     figure(3); clf; hold on;
%     plot(l_list, pow_list);
end