clc
clear

% the ideal system

tic

n_devices = 5; % the number of users
n_channels = 3;  % the number of channels
radius = 500;  % radius of the cell
p_noise = 180000 * 10^(-17.4); % noise power in mW,about -121dBm
shadow_factor = 6; % shadowing factor

p = 10; % transmit power in dBm£¬20dBw=100mW,0dBm=1mW,-10dBm=0.1mW

SE_th = 1; % SE threshold, 1bit/s/Hz

mento=5000; % simulation times

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
                snr_temp=10^(p/10) * (h_large * h_small) / p_noise; % calculate SNR

                SE(i, j) = log2(1 + snr_temp); % calculate the SE
                if SE(i, j) < SE_th  % if the SE is smaller than the threshold, set it to 0.
                   SE(i,j) = 0;
                end

            end
        end

        %% use the Hungarian algorithm to solve the optimal channel assignment varibles, obtaining the maximum sum SE
        [alpha, cost] = Hungarian(-SE);

        SE_mento = [SE_mento, -cost]; % save the maximum sum SE
        
    end
    SE_results = [SE_results, mean(SE_mento)]; % obtain the average sum SE
end


figure; 
plot(SE_results/40); % plot the equavalent S-SE of ideal system




toc