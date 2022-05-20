clc
clear

% S-SE versus the transforming facor

compression_ratio = 0.45:0.05:1;
t_factor = 5*8*compression_ratio; %  transforming factor

%% semantic system
semantic = ones(1, length(compression_ratio))*1.17779;
figure;
plot(t_factor, semantic);

%% Ideal system
shannon=0.79299*40./t_factor;
figure;
plot(t_factor, shannon);

%% 4G
G4 = 0.46323*40./t_factor;
figure;
plot(t_factor, G4);

%% 5G
G5 = 0.55966*40./t_factor;
figure;
plot(t_factor, G5);



