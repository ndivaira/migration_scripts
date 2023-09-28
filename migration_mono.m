clear all
%close all

case_name = "validation_mono";
directory = 'validation';
file_name = 'particles';

% set timestep variables for looping
rate = 10000;
stepNum = 10000:rate:10000000; 

dt_LBM = 2.6666667e-8; %[0.16666667e-6, 0.16666667e-6, 0.16666667e-6];
dt_DEM = 2.6666667e-8; %[0.16666667e-6, 0.16666667e-6, 0.16666667e-6];
dx = 0.4e-3; %lattice spacing, mm
nu = 1/6*dx^2./dt_LBM;
aw =  2.1333e-3; % wall particle radius
aw_factor = 1.75;
X = 38.4e-3; % mm
Y = 25.6e-3; % mm
Z = 104.8e-3; % mm
W = Z - 2*dx - 2*aw_factor*aw;

% set up spatial bins for measuring particle quantities
N_spat_bins = 19;
bins = linspace(dx+aw_factor*aw+0.00001,Z-dx-aw_factor*aw-0.00001,N_spat_bins+1);
bin_width = W/N_spat_bins;
bin_vol = bin_width*X*Y;
total_vol = bin_vol*N_spat_bins;
bin_centres = linspace(dx+aw_factor*aw+bin_width/2,Z-dx-aw_factor*aw-bin_width/2,N_spat_bins);
bin_centres_norm = (bin_centres-dx-aw_factor*aw)/((max(bins)-dx-aw_factor*aw));

a = 2.24e-3;

rs = [1.23e-3 1.4e-3 1.34e-3];
rl = [2.3e-3 3.11e-3 3.66e-3];

%set up plots
markers = ["x","x","x","x","x","x","x"];
colours = ["b", "r", "g", "c", "m", "y", "k"];
%figure('units','normalized','outerposition',[0 0 1 1])
figure(1)
subplot(2,2,1)
hold on
ylabel('Cj')
xlabel('d (\mum)')
box on
hold off
subplot(2,2,2)
hold on
ylabel('C')
xlabel('Lx/W')
%xlabel('Lx')
box on
hold off
subplot(2,2,3)
hold on
ylabel('Lx/Lx_{max}')
xlabel('d (\mum)')
box on
%xlim([0 400])
subplot(2,2,4)
hold on
ylabel('Ly (mm)')
xlabel('Lx (mm)')
box on
hold off

% initialise lists for particle counters, must be within case loop
count = zeros(1,length(stepNum));
bin_count = zeros(N_spat_bins,1,length(stepNum));
bin_count_vol = zeros(N_spat_bins,1,length(stepNum));
bin_vel_x = zeros(N_spat_bins,1,length(stepNum));
bin_vel_y = zeros(N_spat_bins,1,length(stepNum));
inner_Cj = zeros(1,length(stepNum)); % for storing the inner sum of scalar dispersion function (C) for each part size at each timestep
median_vels = zeros(1,length(stepNum));
% iterate through timesteps

for t = 1:length(stepNum)
    velY_total = 0;
    % read in data from LIGGGHTS dump file
    FID=fopen(strcat(directory,'/',case_name,'/',file_name,'_',num2str(stepNum(t))));
    datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 9, 'CollectOutput', 1);
    fclose(FID);
    A = datacell{1};
    
    for p = 1:length(A)
        if A(p,5) > min(bins) && A(p,5) < max(bins)
            if A(p,12) == a
                count(t) = count(t) + 1;
            end
        end
    end

    %---------------- median particle velocity ---------------
    part_vels = [];
    for p = 1:length(A)
        p_rad = A(p,12);
        p_vel_x = A(p,6);
        if p_rad == a && p_vel_x ~= 0
            part_vels = [part_vels p_vel_x];
        end
    end
    median_vels(t) = median(part_vels);
    Q1_vels(t) = quantile(part_vels,0.35);
    Q3_vels(t) = quantile(part_vels,0.65);
    
    %iterate through particles
    for p = 1:length(A)
        p_pos = A(p,5);
        p_rad = A(p,12);
        p_vel_x = A(p,6);
        p_vel_y = A(p,7);
        % calculate inner scalar dispersion function of current particle j,
        % and sum with appropriate particle size
        if p_pos > bins(1) && p_pos < bins(end)
            inner_Cj(t) = inner_Cj(t) + abs(p_pos/W - 0.5);
        end
        %iterate through bins
        for j = 1:N_spat_bins
            %set lower and upper threshold of current bin
            bin_min = bins(j);
            bin_max = bins(j+1);
            % check if current particle within current bin (bin k)
            if p_pos > bin_min && p_pos < bin_max
                % calculate particle velocities in each bin
                bin_vel_x(j,t) = bin_vel_x(j,t) + p_vel_x;
                bin_vel_y(j,t) = bin_vel_y(j,t) + p_vel_y;
                %add to particle count in each bin
                bin_count(j,t) = bin_count(j,t) + 1;
                %spatial bin j+1
                vol_j1 = 0;
                % ensure bin k isn't last bin
                if j < N_spat_bins
                    % calculate particle overlap into bin k+1 (h_k1)
                    h_j1 = p_pos + p_rad - bin_max;
                    % estimate volume of current particle occupying bin k+1
                    if h_j1 > 0
                        vol_j1 = pi*h_j1^3*(p_rad/h_j1-1/3);
                    end
                    bin_count_vol(j+1,t) = bin_count_vol(j+1,t) + vol_j1;
                end
                %bin k-1
                vol_j_1 = 0;
                % ensure bin k isn't first bin
                if j > 1    
                    % calculate particle overlap into bin k-1 (h_k_1)
                    h_j_1 = bin_min - (p_pos - p_rad);
                    % estimate volume of current particle occupying bin k-1
                    if h_j_1 > 0
                        vol_j_1 = pi*h_j_1^3*(p_rad/h_j_1-1/3);
                    end
                    bin_count_vol(j-1,t) = bin_count_vol(j-1,t) + vol_j_1;
                end
                % calculate volume of current particle in current bin
                bin_count_vol(j,t) = bin_count_vol(j,t) + 4/3*pi*p_rad^3 - vol_j1 - vol_j_1;
            end
        end
    end
    % mean settling velocity
    part_total = sum(sum(bin_count(:,t)));
    velY_mean(t) = velY_total/part_total;
end

% simplified phi calculation for checking
total_particle_volume = sum(4/3*pi*(A(A(:,6)>0,12)).^3);
phi_total_check = total_particle_volume/(X*Y*(Z-2*dx-2*1.75*aw));

% calculate mean bin and channel velocities, then mean channel Re and Q
bin_vel_x_mean = bin_vel_x(:,:)./bin_count(:,:);
channel_vel_mean = squeeze(mean(bin_vel_x_mean, 'omitnan'));
vel_profile_norm = mean(bin_vel_x_mean(:,end), 2, 'omitnan')/channel_vel_mean(end);
Re_c = channel_vel_mean*W/nu;
Q = channel_vel_mean*X*W;
Re_c_end = Re_c(end);
%Re_p = 4/3/nu*part_bis(1)^3/(H/2)^2*channel_vel_mean;
Q_end_m3s = Q(end)/1000^3/(X/1000);

%calculate number of particles in each bin
part_count = (sum(bin_count(:,end)));
%calculate phi of each size species j
phi = count.*(4/3*pi*a^3)/total_vol;
phi = phi(end)

% phi for each spatial bin to plot phi distribution
bin_phi = sum(bin_count_vol(:,end)/bin_vol,2);

% calculate characteristic length particles have travelled along channel at each timestep,
% defined as in Chun et al. (2019)
part_mean_vel_x =  sum(bin_vel_x,1)./part_count;
part_mean_vel_y =  sum(bin_vel_y,1)./part_count;
Lx = zeros(1,length(stepNum));
for t = 1:length(stepNum)
    Lx(t) = sum(channel_vel_mean(1:t))*rate*dt_DEM;
end
Lx_W = Lx/W;

%scalar dispersion functions
Cj = inner_Cj./count;

figure(1)
subplot(2,2,2)
hold on
plot(Lx_W(1:5:end), Cj(1:5:end), 'bo', 'markerfacecolor', 'b')
legend({'Total', 'j = 1-2\mum', 'j = 2-3\mum', 'j = 3-4\mum', 'j = 4-5\mum'}, 'Location', 'southwest')
hold off

subplot(2,2,3)
hold on
plot(bin_centres_norm, vel_profile_norm, 'ko', 'markerfacecolor', 'k')
legend()
hold off

subplot(2,2,4)
hold on
plot(bin_centres_norm, bin_phi, 'bo', 'markerfacecolor', 'b')
hold off

