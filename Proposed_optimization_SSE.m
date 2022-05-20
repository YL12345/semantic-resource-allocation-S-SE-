clc
clear

% the proposed method

tic

n_devices = 5; % the number of users
n_channels = 3;  % the number of channels
radius = 500;  % radius of the cell
p_noise = 180000 * 10^(-17.4); % noise power in mW,about -121dBm
shadow_factor = 6; % shadowing factor

p = 10; % transmit power in dBm，20dBw=100mW,0dBm=1mW,-10dBm=0.1mW

n_sym = 1:1:20; % the number of transmitted semantic symbols per word
snr_range = -10:1:20; % the SNR range of the performance table of DeepSC
load('sem_table.mat');% the performance table of DeepSC: row--n_sym; column--snr_range

%香农系统中SE为1bits/s/Hz；按transforming facor=40转换，应为1/40=0.025; 
t_factor = 40;  % the transforming facor

f_th = 0.9; % semantic similarity threshold
s_th = 1/t_factor; % semantic spectral efficiency threshold

mento=1000; % simulation times

SE_results=[]; % save the S-SE for different numbers of channels
for n_channels = 1:1:10 % the number of channels
    SE_mento=[]; % save the optimal S-SE in every simulation
    for sim_times = 1:1:mento % for each simulation
        %% Large-scale fading, including pathloss and shadowing
        d = zeros(n_devices, 1); % save the distance from n_devices devices to the BS
        h_large_scale = zeros(n_devices, 1); % save large-scale channel gain

        position = zeros(n_devices, 2); % save the locations of users
        % generate the locations of users randomly with the BS at (0,0)
        radius_dev = radius * sqrt(rand(n_devices, 1)); 
        phase = rand(n_devices, 1) * 2 * pi; 
        position(:, 1) = radius_dev.*cos(phase);  % x-coordinate
        position(:, 2) = radius_dev.*sin(phase);  % y-coordinate
        for i = 1:n_devices
            d(i) = sqrt(position(i,1)^2 + position(i,2)^2); % the distance from the i-th user to the BS
            pl = 128.1 + 37.6 * log10(d(i)/1000); % pathloss
            h_large_scale(i) = 10^(-(pl + shadow_factor)/10); % large-scale channel gain
        end

        %% Small-scale fading: Rayleigh fading
        h_real = randn(n_devices, n_channels); % Gaussion distribution
        h_image = randn(n_devices, n_channels); % Gaussion distribution
        h_small_scale = (h_real.^2 + h_image.^2)/2; % small-scale fading channel gain of all users over all channels
        
       %% Calculate the maximum S-SE of every user over all channels
        SE = zeros(n_devices, n_channels); % save the maximum S-SE of all users over all channels
        sym_results = zeros(n_devices, n_channels); % save the corresponding number of transmitted semantic symbols
        f_n_results = zeros(n_devices, n_channels); % save the corresponding semantic similarity

        snr_temp = [];

        SE_temp = zeros(1,length(length(n_sym))); % save the S-SE of each user with different number of transmitted semantic symbols
        f_temp = zeros(1,length(length(n_sym))); % save the corresponding semantic similarity

        for i = 1:n_devices  % for each user
            h_large = h_large_scale(i);  % the large-scale fading channel gain of this user
            for j = 1: n_channels  % for each channel
                h_small = h_small_scale(i, j);  % obtain the small-scale fading channel gain of user i over channel j
                for p_index = 1:length(p)  % for each transmit power; here, we set the transmit power as a fixed value
                    snr = 10* log10(10^(p(p_index)/10) * h_large * h_small / p_noise); % calculate SNR
                    snr_temp = [snr_temp, snr];
                    if snr < min(snr_range) % if SNR is smaller than the minimum value of snr_range, set it as the minimum value, the semantic similarity with which will be certainly lower than the threshold
                       snr = min(snr_range); 
                    elseif snr > max(snr_range) % if SNR is larger than the maximum value, set it as the maximum value as the semantic similarity keeps stable with snr keeping increasing
                       snr = max(snr_range); 
                    end
                    % use the exhaustive searching method to solve the optimal number of transmitted semantic symbols, obtaining the corresponding maximum S-SE
                    for sym_index = 1:length(n_sym) 
                        snr_index = round(snr) - min(snr_range) + 1; % obtain the index of snr in sem_table
                        f_n = sem_table(n_sym(sym_index), snr_index); % obtain the semantic similairy through the looking-up table method
                        if f_n < f_th % if the semantic similarity is smaller than the threshold, set it to 0
                            f_n = 0;
                        end
                        f_temp(p_index, sym_index) = f_n;
                        SE_n = f_n / n_sym(sym_index); % calculate the S-SE
                        if SE_n < s_th % if the S-SE is smaller than the threshold, set it to 0.
                           SE_n = 0; 
                        end
                        SE_temp(p_index, sym_index) = SE_n;
                    end
                end
                [temp, index_row] = max(SE_temp); % find the maximum S-SE 
                SE(i, j) = temp; % the maximum SE of user i in channel j

                f_n_results(i, j) = f_temp(1, index_row); % save the corresponding semantic similarity
                sym_results(i, j) = n_sym(index_row); % save the corresponding optimal number of transmitted semantic symbols
            end
        end

        %% use the Hungarian algorithm to solve the optimal channel assignment varibles, obtaining the maximum sum S-SE
        [alpha, cost] = Hungarian(-SE);
        SE_mento = [SE_mento, -cost]; % save the maximum sum S-SE
    end
    SE_results = [SE_results, mean(SE_mento)]; % obtain the average sum S-SE
end

figure;
plot(SE_results);

toc