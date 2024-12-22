% Main script for wireless channel simulation
function main_wireless_simulation()
    % Initialize parameters
    params = initialize_parameters();

    % Run path loss analysis
    run_path_loss_analysis(params);

    % Run fading analysis
    run_fading_analysis(params);

    % Run multipath simulation
    run_multipath_simulation();

    % Run BER analysis
    run_ber_analysis();
end

function params = initialize_parameters()
    params.frequency = 2.4e9;           % Signal frequency in Hz
    params.speed_light = 3e8;           % Speed of electromagnetic waves in m/s
    params.reference_dist = 10;         % Reference distance in meters
    params.transmit_power_dBm = 0;      % Transmission power in dBm
    params.range_distances = 1:1:1000;  % Distance range in meters
    params.path_loss_exp = 3;           % Path loss exponent for urban areas
    params.shadow_std_dB = 4;           % Shadowing standard deviation in dB
end

function run_path_loss_analysis(params)
    % Calculate path loss components
    [path_loss_dB, total_path_loss, combined_path_loss] = calculate_path_loss(params);

    % Calculate received signal levels
    [received_signal_dBm, received_signal_combined, received_signal_no_shadow] = ...
        calculate_received_signals(params, path_loss_dB, total_path_loss, combined_path_loss);

    % Plot results
    plot_path_loss_results(params, path_loss_dB, total_path_loss, combined_path_loss, ...
        received_signal_dBm, received_signal_combined, received_signal_no_shadow);
end

function [path_loss_dB, total_path_loss, combined_path_loss] = calculate_path_loss(params)
    % Calculate reference path loss
    wavelength = params.speed_light / params.frequency;
    ref_path_loss_dB = 10 * log10((4 * pi * params.reference_dist / wavelength)^2);

    % Calculate path loss components
    path_loss_dB = ref_path_loss_dB + 10 * params.path_loss_exp * ...
        log10(params.range_distances / params.reference_dist);

    % Add shadowing effect
    shadow_effect_dB = params.shadow_std_dB * randn(size(params.range_distances));
    total_path_loss = path_loss_dB + shadow_effect_dB;

    % Add small-scale fading
    small_fading_dB = 20 * log10(abs(randn(size(params.range_distances)) + ...
        1i * randn(size(params.range_distances))) / sqrt(2));
    combined_path_loss = total_path_loss + small_fading_dB;
end

function [received_signal_dBm, received_signal_combined, received_signal_no_shadow] = ...
    calculate_received_signals(params, path_loss_dB, total_path_loss, combined_path_loss)

    received_signal_dBm = params.transmit_power_dBm - total_path_loss;
    received_signal_combined = params.transmit_power_dBm - combined_path_loss;
    received_signal_no_shadow = params.transmit_power_dBm - path_loss_dB;
end

function plot_path_loss_results(params, path_loss_dB, total_path_loss, combined_path_loss, ...
    received_signal_dBm, received_signal_combined, received_signal_no_shadow)

    % Plot 1: Signal Strength (Linear Scale)
    create_figure(params.range_distances, received_signal_dBm, ...
        params.transmit_power_dBm - path_loss_dB, ...
        'Signal Strength over Distance (Linear Scale)', ...
        'Distance (m)', 'Signal Level (dBm)', ...
        {'Shadowing', 'No Shadowing'});

    % Plot 2: Path Loss Comparison
    create_figure(params.range_distances, total_path_loss, path_loss_dB, ...
        'Path Loss vs Distance', ...
        'Distance (m)', 'Path Loss (dB)', ...
        {'Shadowing', 'No Shadowing'});

    % Plot 3: Combined Path Loss
    create_figure(params.range_distances, combined_path_loss, path_loss_dB, ...
        'Combined Path Loss vs No Shadowing', ...
        'Distance (m)', 'Path Loss (dB)', ...
        {'Large + Small-Scale Fading', 'No Shadowing'});

    % Plot 4: Signal Strength with Combined Fading
    create_figure(params.range_distances, received_signal_combined, received_signal_no_shadow, ...
        'Signal Strength over Distance (Linear Scale with Combined Fading)', ...
        'Distance (m)', 'Signal Level (dBm)', ...
        {'Combined Fading', 'No Shadowing'});
end

function run_fading_analysis(params)
    % Rayleigh fading analysis
    analyze_rayleigh_fading(100000);

    % Rician fading analysis
    analyze_rician_fading(10000, 10);
end

function analyze_rayleigh_fading(num_samples)
    % Generate and plot Rayleigh fading
    complex_signal = generate_complex_gaussian(num_samples);
    rayleigh_magnitude = abs(complex_signal / sqrt(2));

    plot_distribution(rayleigh_magnitude, 'PDF of Rayleigh Distribution', ...
        'Signal Amplitude', 'Probability Density', 'Empirical PDF');
end

function analyze_rician_fading(num_samples, k_factor)
    % Generate Rician fading parameters
    [rician_magnitude, amplitude_range, pdf_baseline] = generate_rician_fading(num_samples, k_factor);

    % Plot Rician distribution with fixed legend
    figure;

    % Plot histogram first
    [counts, centers] = hist(rician_magnitude, 50);
    pdf_values = counts / (sum(counts) * (centers(2) - centers(1)));
    bar(centers, pdf_values, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'none');
    hold on;

    % Plot theoretical Rayleigh PDF
    plot(amplitude_range, pdf_baseline, '-r', 'LineWidth', 2);

    % Configure plot
    title(['Rician Fading with Comparison (K = ' num2str(k_factor) ')']);
    xlabel('Signal Amplitude');
    ylabel('Probability Density Function (PDF)');
    legend('Rician Distribution', 'Theoretical Rayleigh PDF');
    grid on;
end


function complex_signal = generate_complex_gaussian(num_samples)
    real_part = randn(1, num_samples);
    imaginary_part = randn(1, num_samples);
    complex_signal = real_part + 1i * imaginary_part;
end

function [rician_magnitude, amplitude_range, pdf_baseline] = generate_rician_fading(num_samples, k_factor)
    std_rayleigh = 1 / sqrt(2);
    amplitude_range = linspace(0, 5, 1000);

    % Calculate baseline PDF (Rayleigh theoretical)
    pdf_baseline = (amplitude_range / std_rayleigh^2) .* ...
        exp(-amplitude_range.^2 / (2 * std_rayleigh^2));

    % Generate Rician signal
    line_of_sight = sqrt(k_factor / (k_factor + 1));
    std_rician = sqrt(1 / (2 * (k_factor + 1)));

    % Generate complex Rician signal
    complex_signal = line_of_sight + ...
        (std_rician * randn(1, num_samples) + 1i * std_rician * randn(1, num_samples));
    rician_magnitude = abs(complex_signal);
end

function plot_distribution(magnitude, title_text, xlabel_text, ylabel_text, legend_text)
    figure;
    [counts, centers] = hist(magnitude, 50);
    pdf_values = counts / (sum(counts) * (centers(2) - centers(1)));
    bar(centers, pdf_values, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'none');
    title(title_text);
    xlabel(xlabel_text);
    ylabel(ylabel_text);
    legend(legend_text);
end

function create_figure(x, y1, y2, title_text, xlabel_text, ylabel_text, legend_labels)
    figure;
    plot(x, y1, 'b', 'LineWidth', 1.5); hold on;
    plot(x, y2, 'g-', 'LineWidth', 1.5);
    grid on;
    xlabel(xlabel_text);
    ylabel(ylabel_text);
    title(title_text);
    legend(legend_labels);
end

function run_multipath_simulation()
    % Parameters
    params.N_paths = 5;
    params.path_delays = [0, 1e-6, 2e-6, 3e-6, 4e-6];
    params.Fs = 1e6;
    params.t = (0:1/params.Fs:1e-3)';
    params.x = cos(2 * pi * 1e3 * params.t);
    params.fading_type = 'rayleigh';
    params.K = 3;

    % Run simulation
    [y_total, y_paths] = simulate_multipath(params);

    % Plot results
    plot_multipath_results(params, y_total, y_paths);
end

function [y_total, y_paths] = simulate_multipath(params)
    y_total = zeros(size(params.t));
    y_paths = zeros(length(params.t), params.N_paths);

    for i = 1:params.N_paths
        [y_paths(:,i), y_total] = process_path(params, i, y_total);
    end
end

function [delayed_signal, y_total] = process_path(params, path_index, y_total)
    delay_samples = round(params.path_delays(path_index) * params.Fs);
    h_i = generate_fading_coefficient(params.fading_type, params.K);

    faded_signal = real(h_i) * params.x;
    delayed_signal = [zeros(delay_samples, 1); faded_signal(1:end-delay_samples)];
    y_total = y_total + delayed_signal;
end

function h_i = generate_fading_coefficient(fading_type, K)
    if strcmp(fading_type, 'rayleigh')
        h_i = (randn + 1j * randn) / sqrt(2);
    else % rician
        LOS = sqrt(K / (K + 1));
        NLOS = sqrt(1 / (K + 1)) * (randn + 1j * randn) / sqrt(2);
        h_i = LOS + NLOS;
    end
end

function plot_multipath_results(params, y_total, y_paths)
    % Plot individual paths
    figure; hold on;
    for i = 1:params.N_paths
        plot(params.t, y_paths(:,i), 'DisplayName', ['Path ' num2str(i)]);
    end
    plot(params.t, y_total, 'k', 'LineWidth', 1.5, 'DisplayName', 'Combined Signal');
    configure_plot('Multipath Components with Fading', 'Time (s)', 'Amplitude');

    % Plot original vs received signal
    figure;
    plot(params.t, params.x, 'b', 'DisplayName', 'Original Signal'); hold on;
    plot(params.t, y_total, 'r', 'DisplayName', 'Received Signal (Combined)');
    configure_plot('Original vs Received Signal with Multipath Fading', ...
        'Time (s)', 'Amplitude');
end

function configure_plot(title_text, xlabel_text, ylabel_text)
    xlabel(xlabel_text);
    ylabel(ylabel_text);
    title(title_text);
    legend;
    grid on;
end

function run_ber_analysis()
    % Parameters
    ebno_dB_range = 0:1:10;
    ebno_linear = 10.^(ebno_dB_range / 10);

    % Calculate BER
    ber_bpsk = 0.5 * erfc(sqrt(ebno_linear));
    ber_qpsk = erfc(sqrt(ebno_linear / 2)) - 0.25 * erfc(sqrt(ebno_linear)).^2;

    % Plot results
    figure;
    semilogy(ebno_dB_range, ber_bpsk, '-o', 'LineWidth', 1.5, 'DisplayName', 'BPSK');
    hold on;
    semilogy(ebno_dB_range, ber_qpsk, '-x', 'LineWidth', 1.5, 'DisplayName', 'QPSK');

    grid on;
    title('Bit Error Rate (BER) Performance in AWGN');
    xlabel('E_b/N_0 (dB)');
    ylabel('Bit Error Rate (BER)');
    legend('Location', 'southwest');
end
