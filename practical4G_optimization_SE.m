clc
clear

%  the 4G system
tic

n_devices = 5; % the number of users
n_channels = 3;  % the number of channels
radius = 500;  % radius of the cell
p_noise = 180000 * 10^(-17.4); % noise power in mW,about -121dBm
shadow_factor = 6; % shadowing factor

p = 10; % transmit power in dBm£¬20dBw=100mW,0dBm=1mW,-10dBm=0.1mW

SE_th = 1; % SE threshold, 1bit/s/Hz

%% the SNR-CQI mapping table [17] in 4G system (QPSK 16QAM 64QAM turboÂë)
se_actual_RE = [0.15, 0.23, 0.38, 0.6, 0.88, 1.18, 1.48, 1.91, 2.41, 2.73, 3.32, 3.9, 4.52, 5.12, 5.55]; % the SE specified in TS 236.213
% In [17], the unit of SE is bits/s/Hz while that in 3GPP standard is bits/RE, so a transformation is needed.
% a RE consists of 15KHz bandwidth and 1 OFDM symbols, and 14 OFDM symbols occupy 1ms
% se_actual = se_actual * 15 /14;
se_actual = se_actual_RE * 14 /15; % convert SE in bits/RE to SE in bits/s/Hz
% obtain the SNR according to [17]
snr_th = 10*log10((2.^(se_actual*1.25)-1)/0.66);

BLER=0; % BLER should be 0.1 in this case. However, as the error bits are random, we set it as 0 for comparison, so the performance will be a little worse than BLER=0.1

mento=1000; % simulation times

SE_results=[]; % save the SE for different numbers of channels
snr=[];
SE_temp=[];

for n_channels = 1:1:10
    SE_mento=[]; % save the optimal SE in every simulation
    for sim_times = 1:1:mento
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

        %% Calculate the maximum SE of every user over all channels
        SE = zeros(n_devices, n_channels);  % save the maximum S-SE of all users over all channels
        sym_results = zeros(n_devices, n_channels); % save the corresponding number of transmitted semantic symbols
        f_n_results = zeros(n_devices, n_channels); % save the corresponding semantic similarity

        for i = 1:n_devices 
            h_large = h_large_scale(i); 
            for j = 1: n_channels  
                h_small = h_small_scale(i, j); 
                snr_temp=10^(p/10) * (h_large * h_small) / p_noise;  % calculate SNR
                snr_temp = 10 * log10(snr_temp); 
                snr=[snr,snr_temp];
                
                if snr_temp < snr_th(1) % if snr is smaller than the minimum value, SE should be 0
                   SE(i, j) = 0; 
                elseif snr_temp > snr_th(length(snr_th)) % if snr is larger than the maximum value, SE is equal to the largest one
                    SE(i, j) = se_actual(length(snr_th));
                else
                    SE(i, j) = se_actual(length(find(snr_temp - snr_th >0))); % obtain the correponding SE
                end
                
                if SE(i, j) < SE_th % if the SE is smaller than the threshold, set it to 0.
                   SE(i, j) = 0; 
                end
                
            end
        end

        %% use the Hungarian algorithm to solve the optimal channel assignment varibles, obtaining the maximum sum SE
        [alpha, cost] = Hungarian(-SE);
        SE_mento = [SE_mento, -cost]; % obtain the average sum SE
        
    end
    SE_results = [SE_results, mean(SE_mento)];
end


figure; 
plot(SE_results*(1-BLER)/40); % plot the equavalent S-SE of ideal system



toc