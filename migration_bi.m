clear all
%close all

addpath('~/Documents/00_PhD/02_Simulations_Results/plotting');
colours = linspecer(5);

bi_num = 2;
col_num = 1;
case_name = "bi2_phi20";
case_range = 1:1;
cases = [];
for c = case_range
    cases = [cases strcat(case_name,"_",string(c))];
end
directory = 'bi';
file_name = 'particles';

% set timestep variables for looping
rate = 10000;
stepNum = 10000:rate:90000;
plot_interval = 28;

dt_LBM = 2.6666667e-8;
dt_DEM = 2.6666667e-8;
dx = 0.4e-3; %lattice spacing, mm
a_dx = 1.1e6; % acceleration gradient mm-g (mm/s^2) = pressure gradient m-kg (kg/m^2s^2)
nu = 1/6*dx^2./dt_LBM;
aw =  1.0e-3; % wall particle radius
aw_factor = 2;
X = 64e-3; % mm
Y = 38e-3; % mm
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

%rs = [0.827e-3 0.782e-3 0.743e-3];
%rl = [1.692e-3 1.530e-3 1.306e-3];
%rs = [0.827e-3 0.72e-3 0.62e-3];
%rl = [1.692e-3 1.715e-3 1.75e-3];
rs = [1.3e-3 1.47e-3 1.83e-3 1.405e-3 1.53e-3];
rl = [2.7e-3 2.88e-3 3.23e-3 3.05e-3 2.66e-3];
species_sizes = [rs(bi_num), rl(bi_num)]
%species_sizes = [rs(1), rl(1)];
%species_sizes = [0 0];

%set up plots
markers = ["o","s","^","d","v","+","o","o","o","o"];
%figure('units','normalized','outerposition',[0 0 1 1])
figure(10)
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
ylim([0.21 0.28])
box on
hold off
subplot(2,2,3)
hold on
ylabel('u/u_{max}')
xlabel('Lx/W')
box on
%xlim([0 400])
subplot(2,2,4)
hold on
xlabel('Lx/W')
ylabel('\phi')
ylim([0.1 0.55])
box on
hold off

Lx_median = zeros(length(cases),2,length(stepNum));
Lx_Q1 = zeros(length(cases),2,length(stepNum));
Lx_Q3 = zeros(length(cases),2,length(stepNum));
C_seg = zeros(length(cases),length(stepNum));

for c = 1:length(cases)
    case_name = cases(c)
    % initialise lists for particle counters, must be within case loop
    bin_count = zeros(N_spat_bins,2,length(stepNum));
    bin_count_vol = zeros(N_spat_bins,2,length(stepNum));
    bin_vel_x = zeros(N_spat_bins,2,length(stepNum));
    bin_vel_y = zeros(N_spat_bins,2,length(stepNum));
    inner_Cj = zeros(2,length(stepNum)); % for storing the inner sum of scalar dispersion function (C) for each part size at each timestep
    median_vels = zeros(2,length(stepNum));
    % iterate through timesteps

    %count actual numbers of particles @ t=0
    %FID=fopen(strcat(directory,'/',cases(c),'/',file_name,'_10000'));
    FID=fopen(strcat(directory,'/',case_name,'/',file_name,'_1'));
    datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 9, 'CollectOutput', 1);
    fclose(FID);
    A = datacell{1};
    count_j = zeros(1,2);
    for j = 1:2
        for p = 1:length(A)
            %if A(p,5) > min(bins) && A(p,5) < max(bins)
            if A(p,2) == 1
                %if ismember(A(p,12),species_sizes) == 0
                %    species_sizes(j) = A(p,12);
                %end
                if A(p,12) == species_sizes(j)
                    count_j(j) = count_j(j) + 1;
                end
            end
        end
    end


    for t = 1:length(stepNum)
        velY_total = 0;
        % read in data from LIGGGHTS dump file
        FID=fopen(strcat(directory,'/',case_name,'/',file_name,'_',num2str(stepNum(t))));
        datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 9, 'CollectOutput', 1);
        fclose(FID);
        A = datacell{1};
        
        %---------------- median particle velocity of each size ---------------
        for r = 1:2
            part_vels = [];
            for p = 1:length(A)
                p_rad = A(p,12);
                p_vel_x = A(p,6);
                %if p_rad > part_sizes(r)-0.00001 && p_rad < part_sizes(r)+0.00001 && p_vel_x ~= 0
                if p_rad == species_sizes(r) && p_vel_x ~= 0
                    part_vels = [part_vels p_vel_x];
                end
            end
            median_vels(r,t) = median(part_vels);
            Q1_vels(r,t) = quantile(part_vels,0.35);
            Q3_vels(r,t) = quantile(part_vels,0.65);
        end
        
        %iterate through particles
        for p = 1:length(A)
            p_pos = A(p,5);
            p_rad = A(p,12);
            p_vel_x = A(p,6);
            p_vel_y = A(p,7);
            % calculate inner scalar dispersion function of current particle j,
            % and sum with appropriate particle size
            for k = 1:2
                if p_rad == species_sizes(k) 
                    if p_pos > bins(1) && p_pos < bins(end)
                        inner_Cj(k,t) = inner_Cj(k,t) + abs(p_pos/W - 0.5);
                    end
                end
            end
            %iterate through bins
            for j = 1:N_spat_bins
                %set lower and upper threshold of current bin
                bin_min = bins(j);
                bin_max = bins(j+1);
                % check if current particle within current part bin (bin k)
                if p_pos > bin_min && p_pos < bin_max
                    % calculate particle velocities in each bin
                    for k = 1:2
                        if p_rad == species_sizes(k) 
                            bin_vel_x(j,k,t) = bin_vel_x(j,k,t) + p_vel_x;
                            bin_vel_y(j,k,t) = bin_vel_y(j,k,t) + p_vel_y;
                        end
                    end
                    %add to particle count in each bin
                    for k = 1:2
                        if p_rad == species_sizes(k)
                            bin_count(j,k,t) = bin_count(j,k,t) + 1;
                        end
                    end
                    %spatial bin j+1
                    vol_j1 = 0;
                    % ensure bin j isn't last bin
                    if j < N_spat_bins
                        % calculate particle overlap into bin j+1 (h_j1)
                        h_j1 = p_pos + p_rad - bin_max;
                        % estimate volume of current particle occupying bin j+1
                        if h_j1 > 0
                            vol_j1 = pi*h_j1^3*(p_rad/h_j1-1/3);
                        end
                        for k = 1:2
                            if p_rad == species_sizes(k)
                                bin_count_vol(j+1,k,t) = bin_count_vol(j+1,k,t) + vol_j1;
                            end
                        end
                    end
                    %bin j-1
                    vol_j_1 = 0;
                    % ensure bin j isn't first bin
                    if j > 1    
                        % calculate particle overlap into bin j-1 (h_j_1)
                        h_j_1 = bin_min - (p_pos - p_rad);
                        % estimate volume of current particle occupying bin j-1
                        if h_j_1 > 0
                            vol_j_1 = pi*h_j_1^3*(p_rad/h_j_1-1/3);
                        end
                        for k = 1:2
                            if p_rad == species_sizes(k)
                                bin_count_vol(j-1,k,t) = bin_count_vol(j-1,k,t) + vol_j_1;
                            end
                        end
                    end
                    % calculate volume of current particle in current bin
                    for k = 1:2
                        if p_rad == species_sizes(k)
                            bin_count_vol(j,k,t) = bin_count_vol(j,k,t) + 4/3*pi*p_rad^3 - vol_j1 - vol_j_1;
                        end
                    end
                end
            end
        end
        part_total = sum(sum(bin_count(:,:,t)));
    end
  
    % simplified phi calculation for checking
    total_particle_volume = sum(4/3*pi*(A(A(:,6)>0,12)).^3);
    phi_total_check = total_particle_volume/(X*Y*(Z-2*dx-2*1.75*aw));

    % calculate mean bin and channel velocities, then mean channel Re and Q
    bin_vel_x_mean = bin_vel_x(:,:,:)./bin_count(:,:,:);
    channel_vel_mean = squeeze(mean(mean(bin_vel_x_mean, 'omitnan')));
    vel_profile_norm(:,c) = mean(bin_vel_x_mean(:,:,end), 2, 'omitnan')/channel_vel_mean(end);
    Re_c = channel_vel_mean*W/nu;
    Q = channel_vel_mean*X*W;
    Re_c_end = Re_c(end);
    %Re_p = 4/3/nu*part_bins(1)^3/(H/2)^2*channel_vel_mean;
    Q_end_m3s = Q(end)/1000^3/(X/1000);

    nu_eff = a_dx*W^2/12./(channel_vel_mean); % mm^2/s
    Re_c_eff = channel_vel_mean*W./nu_eff;
    Re_c_eff_alt = 12*channel_vel_mean.^2/a_dx/W;

    %calculate number of particles in each bin
    part_count = (sum(bin_count(:,:,end)));
    %calculate phi of each size species j
    phi_j = count_j.*(4/3*pi*species_sizes.^3)/total_vol;
    phi_total = sum(phi_j);

    % phi for each spatial bin to plot phi distribution
    bin_phi_part(:,:,c) = bin_count_vol(:,:,end)/bin_vol;
    bin_phi_total(:,c) = sum(bin_count_vol(:,:,end)/bin_vol,2);

    % calculate characteristic length particles have travelled along channel at each timestep,
    % defined as in Chun et al. (2019)
    part_mean_vel_x =  sum(bin_vel_x,1)./part_count;
    part_mean_vel_y =  sum(bin_vel_y,1)./part_count;
    Lx = zeros(1,length(stepNum));
    Lx_part = zeros(2,length(stepNum));
    Ly_part = zeros(2,length(stepNum));
    for t = 1:length(stepNum)
        Lx(t) = sum(channel_vel_mean(1:t))*rate*dt_DEM;
        Lx_part(:,t) = sum(part_mean_vel_x(:,:,1:t),3)*rate*dt_DEM;
        Ly_part(:,t) = sum(part_mean_vel_y(:,:,1:t),3)*rate*dt_DEM;
        Lx_median(c,:,t) = sum(median_vels(:,1:t),2)*rate*dt_DEM;
        Lx_Q1(c,:,t) = sum(Q1_vels(:,1:t),2)*rate*dt_DEM;
        Lx_Q3(c,:,t) = sum(Q3_vels(:,1:t),2)*rate*dt_DEM;
    end
    Lx_W = Lx/W;
    Lx_W_part = Lx_part/W;
    Ly_mean = mean(Ly_part);
    Ly_Lx = Ly_part./Lx_part;
    Lx_W_part_norm = Lx_W_part(:,end)/max(Lx_W_part(:,end)); %dividing by max travelled distance of a particle size, rather than mean travelled distance of suspension, is more useful, because it is the furthest travelling particles which will enter/block cleats far from the wellbore first

    %scalar dispersion functions
    %calculate scalar dispersion function for each particle size species
    Cj_current = inner_Cj./count_j';
    Cj(c,:,:) = Cj_current;
    %calculate total scaler dispersion function
    for t = 1:length(stepNum)
        C(c,t) = sum(Cj_current(:,t).*phi_j')/phi_total;
    end
    C_seg_inner = zeros(2,length(stepNum));
    for i = 1:2
        for j = 1:2
            C_seg_inner(i,:) = abs(Cj_current(i,:)-Cj_current(j,:)) + C_seg_inner(i,:);
        end
    end
    C_seg_inner = C_seg_inner/2;
    for i = 1:2
        C_seg(c,:) = C_seg(c,:) + C_seg_inner(i,:);
    end
    C_seg(c,:) = 4*C_seg(c,:)/2;
end

bin_phi_part_mean = mean(bin_phi_part,3);
bin_phi_total_mean = mean(bin_phi_total,2);
vel_profile_mean = mean(vel_profile_norm,2);

Cj_median = squeeze(median(Cj(:,:,:),1));
Cj_mean = squeeze(mean(Cj(:,:,:),1));
C_median = squeeze(median(C(:,:,:),1));
C_mean = squeeze(mean(C(:,:,:),1));
%migration function
phi_max = [0.783 0.7395 0.6804 0.783 0.7395 0.6804];
C_sigma = (1-4*C_mean)/(1-phi_total/phi_max(bi_num));
C_seg_median = median(C_seg,1);

subplot(2,2,1)
hold on
plot(2*1000*species_sizes(:), Cj_median(:,end), markers(col_num), 'color', colours(col_num,:), 'MarkerFaceColor', colours(col_num,:))
%plot(2*1000*part_bins(1:end-1), Cj_mean(:,:,end),'ko','MarkerFaceColor','k')
hold off

figure(10)
subplot(2,2,2)
hold on
%plot(Lx_W(1:5:end), C_median(:,1:5:end), 'ko', 'markerfacecolor', 'k')
plot(Lx_W(1:plot_interval:end), C_mean(1:plot_interval:end), markers(col_num), 'color', colours(col_num,:), 'MarkerFaceColor', colours(col_num,:))
%for c = 1:length(cases)
%    plot(Lx_W(1:5:end), C(c,1:5:end), strcat(colours(col_num),'-'))
%end
lines_bi = ["-" "--"];
for i = 1:2
    plot(Lx_W(1:plot_interval:end), Cj_median(i,1:plot_interval:end), lines_bi(i), 'color', colours(col_num,:), 'MarkerFaceColor', colours(col_num,:))
end
legend({'Total', 'j = 1-2\mum', 'j = 2-3\mum', 'j = 3-4\mum', 'j = 4-5\mum'}, 'Location', 'southwest')
hold off

Lx_median_mean = mean(Lx_median(:,:,:),1);
Lx_Q1_mean = mean(Lx_Q1(:,:,:),1);
Lx_Q3_mean = mean(Lx_Q3(:,:,:),1);
Lx_median_med = median(Lx_median(:,:,:),1);
Lx_Q1_med = median(Lx_Q1(:,:,:),1);
Lx_Q3_med = median(Lx_Q3(:,:,:),1);

subplot(2,2,3)
hold on
%errorbar(2*1000*part_bins(1:end-1), Lx_median_mean(:,:,end), Lx_median_mean(:,:,end)-Lx_Q1_mean(:,:,end), abs(Lx_median_mean(:,:,end)-Lx_Q3_mean(:,:,end)), 'ko', 'MarkerEdgeColor',colours(mark),'MarkerFaceColor',colours(mark))
%errorbar(2*1000*part_bins(1:end-1), Lx_median_med(:,:,end), Lx_median_med(:,:,end)-Lx_Q1_med(:,:,end), abs(Lx_median_med(:,:,end)-Lx_Q3_med(:,:,end)), 'bo', 'MarkerEdgeColor','k','MarkerFaceColor','k')
%plot(2*1000*species_sizes(:), Lx_median_med2(:,:,end)/max(Lx_median_med2(:,:,end)),strcat(colours(col_num),markers(col_num)), 'markerfacecolor', colours(col_num))
plot(bin_centres_norm, vel_profile_mean, markers(col_num), 'color', colours(col_num,:), 'MarkerFaceColor', colours(col_num,:));
hold off

subplot(2,2,4)
hold on
Lx_W(1:plot_interval:end)
C_seg_median(1:plot_interval:end)
%plot(Lx_W(1:plot_interval:end), C_seg_median(1:plot_interval:end), markers(col_num), 'color', colours(col_num,:), 'markerfacecolor', colours(col_num,:))
plot(bin_centres_norm, bin_phi_total_mean, markers(col_num), 'color', colours(col_num,:), 'MarkerFaceColor', colours(col_num,:))
%for i = 1:2
%    plot(bin_centres_norm, bin_phi_part(:,i,1), strcat(colours(i),markers(i)), 'markerfacecolor', colours(i))
%end
legend({string(species_sizes(1)),string(species_sizes(2))})
hold off

