clear all
%close all

addpath('~/Documents/00_PhD/02_Simulations_Results/plotting');
colours = linspecer(5);

pd_num = 1;
case_name = "pd2";
case_range = 1:1;
cases = [];
for c = case_range
    cases = [cases strcat(case_name,"_",string(c))];
end
directory = 'pd';
file_name = 'particles';

% set timestep variables for looping
rate = 20000;
stepNum = 20000:rate:5960000;
%stepNum = [1];
plot_interval = 10;

dt_LBM = 2.6666667e-8; %[0.16666667e-6, 0.16666667e-6, 0.16666667e-6];
dt_DEM = 2.6666667e-8; %[0.16666667e-6, 0.16666667e-6, 0.16666667e-6];
dx = 0.4e-3; %lattice spacing, mm
a_dx = 2e6; % acceleration gradient mm-g (mm/s^2) = pressure gradient m-kg (kg/m^2s^2)
nu = 1/6*dx^2./dt_LBM; % mm^2/s
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

% set particle bin sizes
N_part_bins = 6; % must be factor of species
part_min = 1.0e-3;
%part_max = 5.0e-3;
part_max = 4.0e-3;
part_bins = linspace(part_min, part_max, N_part_bins+1);
part_bins_mean = (part_bins(1:end-1) + part_bins(2:end))/2;

%set up plots
%figure('units','normalized','outerposition',[0 0 1 1])
figure(1)
subplot(2,2,1)
hold on
ylabel('C_j')
xlabel('a (\mum)')
box on
hold off
subplot(2,2,2)
hold on
ylabel('C')
xlabel('L/W')
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
ylabel('\phi_{local}')
xlabel('y/W')
box on
hold off

Cj_bins = zeros(length(cases),N_part_bins,length(stepNum));
Lx_median = zeros(length(cases),N_part_bins,length(stepNum));
Lx_Q1 = zeros(length(cases),N_part_bins,length(stepNum));
Lx_Q3 = zeros(length(cases),N_part_bins,length(stepNum));
C_seg = zeros(length(cases),length(stepNum)); 
Re_p = zeros(length(cases),N_part_bins,length(stepNum));

%set up bins for actual discreet particle counting
species=30;
species_sizes = linspace(part_min, part_max, species+1);
species_sizes = species_sizes(1:end-1);

for c = 1:length(cases)
    % initialise lists for particle counters, must be within case loop
    bin_count = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_count_vol = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_vel_x = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_vel_y = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_vel_z = zeros(N_spat_bins,N_part_bins,length(stepNum));
    inner_Cj = zeros(species,length(stepNum)); % for storing the inner sum of scalar dispersion function (C) for each part size at each timestep
    inner_Cj_g = zeros(species,length(stepNum));
    median_vels = zeros(length(part_bins)-1,length(stepNum));
  
    % particle pressures sigma 
    bin_sigmax = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_sigmaz = zeros(N_spat_bins,N_part_bins,length(stepNum));
    bin_mu = zeros(N_spat_bins,length(stepNum));
    
    
    % iterate through timesteps

    %count actual numbers of particles @ t=0
    FID=fopen(strcat(directory,'/',cases(c),'/',file_name,'_1'));
    %FID=fopen(strcat(directory,'/',cases(c),'/',file_name,'_40000'));
    datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 9, 'CollectOutput', 1);
    fclose(FID);
    A = datacell{1};
    count_j = zeros(1,species);
    for j =1:species
        for p = 1:length(A)
            if A(p,5) > min(bins) && A(p,5) < max(bins)
                if A(p,12) > species_sizes(j) - 0.000001 && A(p,12) < species_sizes(j) + 0.000001
                    count_j(j) = count_j(j) + 1;
                end
            end
        end
    end

    for t = 1:length(stepNum)
        velY_total = 0;
        % read in data from LIGGGHTS dump file
        FID=fopen(strcat(directory,'/',cases(c),'/',file_name,'_',num2str(stepNum(t))));
        datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 9, 'CollectOutput', 1);
        fclose(FID);
        A = datacell{1};
        
        %---------------- median particle velocity of each size ---------------
        for r = 1:(length(part_bins)-1)
            part_vels = [];
            for p = 1:length(A)
                p_rad = A(p,12);
                p_vel_x = A(p,6);
                %if p_rad > part_sizes(r)-0.00001 && p_rad < part_sizes(r)+0.00001 && p_vel_x ~= 0
                if p_rad > part_bins(r)-0.00001 && p_rad < part_bins(r+1)+0.00001 && p_vel_x ~= 0
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
            p_vel_z = A(p,8); % fluctuation velocity
            % particle forces and stresses
            p_Fx = A(p,9);
            p_Fz = A(p,11);
            p_sigmax = p_Fx/(pi*p_rad^2);
            p_sigmaz = p_Fz/(pi*p_rad^2);
            % calculate inner scalar dispersion function of current particle j,
            % and sum with appropriate particle size
            for k = 1:species
                if p_rad > species_sizes(k) - 0.000001 && p_rad < species_sizes(k) + 0.000001
                    if p_pos > bins(1) && p_pos < bins(end)
                        inner_Cj(k,t) = inner_Cj(k,t) + abs((p_pos-dx-2*aw)/W - 0.5);
                        inner_Cj_g(k,t) = inner_Cj_g(k,t) + (p_pos-dx-2*aw)/W;
                    end
                end
            end
            %iterate through bins
            for j = 1:N_spat_bins
                %set lower and upper threshold of current bin
                bin_min = bins(j);
                bin_max = bins(j+1);
                % check if current particle within current bin (bin k)
                if p_pos > bin_min && p_pos < bin_max
                    % calculate particle velocities in each bin
                    for k = 1:N_part_bins
                        if p_rad > part_bins(k)-0.000001 && p_rad < part_bins(k+1)-0.000001 %float correction added
                            bin_vel_x(j,k,t) = bin_vel_x(j,k,t) + p_vel_x;
                            bin_vel_y(j,k,t) = bin_vel_y(j,k,t) + p_vel_y;
                            bin_vel_z(j,k,t) = bin_vel_z(j,k,t) + p_vel_z;
                        end
                    end
                    %add to particle count in each bin
                    for k = 1:N_part_bins
                        if p_rad > part_bins(k)-0.000001 && p_rad < part_bins(k+1)-0.000001 %float correction added
                            bin_count(j,k,t) = bin_count(j,k,t) + 1;
                            if bin_min < Z/2 && j > 1
                                bin_sigmaz(j,k,t) = bin_sigmaz(j,k,t) + p_sigmaz;
                            end
                            else if bin_max > Z/2 && j < (N_spat_bins)
                                bin_sigmaz(j,k,t) = bin_sigmaz(j,k,t) - p_sigmaz;
                            end
                        end
                    end

                    % add to stresses in each bin
                    bin_sigmax(j,t) = bin_sigmax(j,t) + p_sigmax;

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
                        for k = 1:N_part_bins
                            if p_rad > part_bins(k)-0.000001 && p_rad < part_bins(k+1)-0.000001 %float correction added
                                bin_count_vol(j+1,k,t) = bin_count_vol(j+1,k,t) + vol_j1;
                            end
                        end
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
                        for k = 1:N_part_bins
                            if p_rad > part_bins(k)-0.000001 && p_rad < part_bins(k+1)-0.000001 %float correction added
                                bin_count_vol(j-1,k,t) = bin_count_vol(j-1,k,t) + vol_j_1;
                            end
                        end
                    end
                    % calculate volume of current particle in current bin
                    for k = 1:N_part_bins
                        if p_rad > part_bins(k)-0.000001 && p_rad < part_bins(k+1)-0.000001 %float correction added
                            bin_count_vol(j,k,t) = bin_count_vol(j,k,t) + 4/3*pi*p_rad^3 - vol_j1 - vol_j_1;
                        end
                    end
                end
            end
        end
        % mean settling velocity
        part_total = sum(sum(bin_count(:,:,t)));
        velY_mean(t) = velY_total/part_total;
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

    nu_eff = a_dx*(W)^2/12./(channel_vel_mean) % mm^2/s
    Re_c_eff = channel_vel_mean*W./nu_eff;
    Re_c_eff_alt = 12*channel_vel_mean.^2/a_dx/W;

    %calculate number of particles in each bin
    part_count = (sum(bin_count(:,:,end)));
    %calculate phi of each size species j
    phi_j = count_j.*(4/3*pi*species_sizes.^3)/total_vol;
    phi_total = sum(phi_j)
    %calculate phi of each particle bin
    species_per_bin = species/N_part_bins;
    for N = 1:N_part_bins
        phi_j_bins(N) = sum(phi_j((N-1)*species_per_bin+1:N*species_per_bin));
    end

    % phi for each spatial bin to plot phi distribution
    bin_phi_part(:,:,c) = bin_count_vol(:,:,end)/bin_vol;  
    bin_phi_total(:,c) = sum(bin_count_vol(:,:,end)/bin_vol,2);


    % calculate stresses and fluxes
    sigmaz = bin_sigmaz./bin_count;
 
    sigmaz_part = zeros(N_part_bins,length(stepNum)) ;
    for i = 1:N_part_bins
        for t = 1:length(stepNum)
            sigmaz_sub = sigmaz(:,i,t);
            sigmaz_part(i,t) = sum(sigmaz_sub(~isinf(sigmaz_sub)),'omitnan');
        end
    end

    tavg_steps = 10;
    sigmaz_part_tavg = zeros(N_part_bins,floor(length(stepNum)/tavg_steps));
    for t = 1:(floor(length(stepNum)/tavg_steps))
        for i = 1:N_part_bins
            sigmaz_part_tavg(i,t) = mean(sigmaz_part(i,t:t+tavg_steps));
        end
    end


    sigmaz_bins = sum(sigmaz,2,'omitnan');



    % calculate characteristic length particles have travelled along channel at each timestep,
    % defined as in Chun et al. (2019)
    part_mean_vel_x =  sum(bin_vel_x,1)./part_count;
    part_mean_vel_y =  sum(bin_vel_y,1)./part_count;
    part_mean_vel_z =  sum(bin_vel_z,1)./part_count;
    Lx = zeros(1,length(stepNum));
    Lx_part = zeros(N_part_bins,length(stepNum));
    Ly_part = zeros(N_part_bins,length(stepNum));
    for t = 1:length(stepNum)
        Lx(t) = sum(channel_vel_mean(1:t))*rate*dt_DEM;
        Lx_part(:,t) = sum(part_mean_vel_x(:,:,1:t),3)*rate*dt_DEM;
        Ly_part(:,t) = sum(part_mean_vel_y(:,:,1:t),3)*rate*dt_DEM;
        Lx_median(c,:,t) = sum(median_vels(:,1:t),2)*rate*dt_DEM;
        Lx_Q1(c,:,t) = sum(Q1_vels(:,1:t),2)*rate*dt_DEM;
        Lx_Q3(c,:,t) = sum(Q3_vels(:,1:t),2)*rate*dt_DEM;
    end
    Lx_median2(c,:) = sum(median_vels(:,round(length(stepNum)*0.6):t),2)*rate*dt_DEM;
    %L = sum(channel_vel_mean*rate*dt)
    Lx_W = Lx/W;
    Lx_W_part = Lx_part/W;
    Ly_mean = mean(Ly_part);
    Ly_Lx = Ly_part./Lx_part;
    Ly_Lx(:,end);
    Lx_W(end)
    Lx_W_part_norm = Lx_W_part(:,end)/max(Lx_W_part(:,end)); %dividing by max travelled distance of a particle size, rathan mean travelled distance of suspension, is more useful, because it is the furthest travelling particles which with enter/block cleats far from the wellbore first
    %Lx_W_part_norm = Lx_W_part(:,end)./Lx_W(end)

    %part_mean_vel = zeros(N_part_bins,length(part_mean_vel_x))
    %for i = length(part_mean_vel_x)
    %    part_mean_vel = squeeze(part_mean_vel_x);
    %    Re_p = part_mean_vel(:,i).*part_bins.^3;
    %end
    Re_p(c,:,:) = 4/3*squeeze(part_mean_vel_x).*reshape(part_bins_mean.^3,6,1)/(W/2)^2./reshape(nu_eff, 1, length(nu_eff));
    
    
    %scalar dispersion functions
    %calculate scalar dispersion function for each particle size species
    Cj = inner_Cj./count_j';
    Cj_g = inner_Cj_g./count_j';
    %calculate Cj for each particle bin
    for i = 1:N_part_bins
        Cj_bins(c,i,:) = sum(inner_Cj((i-1)*species_per_bin+1:i*species_per_bin, :))/part_count(i);
        Cj_g_bins(c,i,:) = sum(inner_Cj_g((i-1)*species_per_bin+1:i*species_per_bin, :))/part_count(i);
    end
    %calculate total scaler dispersion function
    for t = 1:length(stepNum)
        C(c,t) = sum(Cj(:,t).*phi_j')/phi_total;
        C_g(c,t) = sum(Cj_g(:,t).*phi_j')/phi_total;
    end

    
    C_seg_inner = zeros(species,length(stepNum));
    for i = 1:species
        for j = 1:species
            C_seg_inner(i,:) = abs(Cj(i,:)-Cj(j,:)) + C_seg_inner(i,:);
        end
    end
    C_seg_inner = C_seg_inner/species;
    for i = 1:species
        C_seg(c,:) = C_seg(c,:) + C_seg_inner(i,:);
    end
    C_seg(c,:) = 4*C_seg(c,:)/species;
    %{
    C_seg_inner = zeros(N_part_bins,length(stepNum));
    for i = 1:N_part_bins
        for j = 1:N_part_bins
            C_seg_inner(i,:) = squeeze(abs(Cj_bins(c,i,:)-Cj_bins(c,j,:)))' + C_seg_inner(i,:);
        end
    end
    C_seg_inner = C_seg_inner/N_part_bins;
    for i = 1:N_part_bins
        C_seg(c,:) = C_seg(c,:) + C_seg_inner(i,:);
    end
    C_seg(c,:) = 4*C_seg(c,:)/N_part_bins;
    %} 
end

bin_phi_part_mean = mean(bin_phi_part,3);
bin_phi_total_mean = mean(bin_phi_total,2);
vel_profile_mean = mean(vel_profile_norm,2);

%averaging over cases
Cj_median = squeeze(median(Cj_bins(:,:,:),1));
Cj_mean = squeeze(mean(Cj_bins(:,:,:),1));
Cj_g_median = squeeze(median(Cj_g_bins(:,:,:),1));
C_median = squeeze(median(C(:,:,:),1));
C_mean = squeeze(mean(C(:,:,:),1));
C_g_mean = squeeze(mean(C_g(:,:,:),1));
%migration function
%phi_max = [0.783 0.7395 0.6804];
phi_max = [0.7395 0.7395 0.7395 0.7395 0.7395 0.7395 0.7395];
C_sigma = (1-4*C_mean)/(1-phi_total/phi_max(pd_num));
C_seg_median = median(C_seg,1);

Re_p_median = median(Re_p,1);

%Cj statistics for boxplot
Cj_bins
Cj_quantiles = quantile(Cj_bins(:,:,:),[0, 0.25, 0.5, 0.75, 1],1)

markers = ["o","s","^","d","v","+","o","o","o","o"];

subplot(2,2,1)
hold on
plot(1000*(part_bins(1:end-1)+(part_bins(end)-part_bins(end-1))/2), Cj_median(:,end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
%plot(2*1000*(part_bins(1:end-1)+(part_bins(end)-part_bins(end-1))/2), Cj_g_median(:,end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
legend('Location', 'northeast')
xlim([1 4])
hold off

subplot(2,2,2)
hold on
plot(Lx_W(1:plot_interval:end), C_median(1:plot_interval:end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
%plot(Lx_W(1:5:end), C_mean(1:5:end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
%plot(Lx_W(1:5:end), C_g_mean(1:5:end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
for c = 1:length(cases)
    %plot(Lx_W(1:5:end), C(c,1:5:end), strcat(colours(pd_num),'-'))
end
for i = 1:N_part_bins
    %plot(Lx_W(1:plot_interval:end), Cj_median(i,1:plot_interval:end),'-')
    %plot(Lx_W, Cj_median(i,1:end),'-','color', colours(i,:))
end
%legend({'Distribution 1', 'Distribution 2', 'Distribution 3'}, 'Location', 'northeast')
%legend({'G=0.5MPa/m','G=2MPa/m','G=4MPa/m'},'location','northeast')
legend('location','northeast')
hold off

%{
Lx_median_mean = mean(Lx_median(:,:,:),1);
Lx_Q1_mean = mean(Lx_Q1(:,:,:),1);
Lx_Q3_mean = mean(Lx_Q3(:,:,:),1);
Lx_median_med = median(Lx_median(:,:,:),1);
Lx_median_med2 = median(Lx_median2(:,:,:),1);
Lx_Q1_med = median(Lx_Q1(:,:,:),1);
Lx_Q3_med = median(Lx_Q3(:,:,:),1);

subplot(2,2,3)
hold on
errorbar(2*1000*part_bins(1:end-1), Lx_median_mean(:,:,end), Lx_median_mean(:,:,end)-Lx_Q1_mean(:,:,end), abs(Lx_median_mean(:,:,end)-Lx_Q3_mean(:,:,end)), 'ko', 'MarkerEdgeColor',colours(pd_num),'MarkerFaceColor',colours(pd_num))
errorbar(2*1000*part_bins(1:end-1), Lx_median_med(:,:,end), Lx_median_med(:,:,end)-Lx_Q1_med(:,:,end), abs(Lx_median_med(:,:,end)-Lx_Q3_med(:,:,end)), 'bo', 'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(2*1000*part_bins(1:end-1), Lx_median_med2(:,:,end)/max(Lx_median_med2(:,:,end)),strcat(colours(pd_num),markers(pd_num)), 'markerfacecolor', colours(pd_num))
hold off
%}

subplot(2,2,3)
hold on
y = linspace(0,1,100);
u = y.*(1-y);
%plot(bin_centres_norm, vel_profile_mean, markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
%for i = 1:N_part_bins
%    plot(bin_centres_norm, bin_vel_x_mean(:,i,end), strcat(colours(i),markers(i)), 'markerfacecolor', colours(i))
%end
plot(y,u/mean(u),'k-')
%plot([0.4 0.4],[-2 2],'c--')
%plot([0.6 0.6],[-2 2],'c--')
%plot([0.47 0.47],[-2 2],'g--')
%plot([0.53 0.53],[-2 2],'g--')
ylabel('u/\langle u \rangle')
xlabel('y/W')
legend('Location', 'northeast')
%legend({'\phi=0.2','\phi=0.3','\phi=0.4','\phi=0.5'},'location','northeast')
ylim([0 1.6])
hold off

subplot(2,2,4)
hold on
%plot(Lx_W(1:10:end), C_seg_median(1:10:end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
plot(bin_centres_norm, bin_phi_total_mean, markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))
for i = 1:N_part_bins
    %plot(bin_centres_norm, bin_phi_part_mean(:,i), markers(i), 'color', colours(i,:), 'MarkerFaceColor', colours(i,:))
    %plot(Lx_W(1:10:end), sigmaz(1:10:end,i,), markers(i), 'color', colours(i,:), 'MarkerFaceColor', colours(i,:))
end
%legend({'Distribution 1','Distribution 2','Distribution 3'},'Location','northeast')
legend('location', 'northeast')
hold off

%{
subplot(2,2,4)
hold on
plot(Lx_W(1:5:end), C_sigma(1:5:end), strcat(colours(pd_num),markers(pd_num)), 'markerfacecolor', colours(pd_num))
plot(Lx, Ly_mean, strcat(markers(pd_num),colours(pd_num)))
legend({'Total', 'j = 1-2\mum', 'j = 2-3\mum', 'j = 3-4\mum', 'j = 4-5\mum'}, 'Location', 'southwest')
hold off
%}

%{
figure(2)
subplot(2,2,1)
plot(bin_centres_norm, sigmaz_bins(:,end), markers(pd_num), 'color', colours(pd_num,:), 'MarkerFaceColor', colours(pd_num,:))

subplot(2,2,2)
hold on
for i = 1:N_part_bins
    %plot(Lx_W(1:10:end), sigmaz_part(i,1:10:end),'-', 'color', colours(i,:))
    %plot(Lx_W(1:10:end), sigmaz_part_tavg(i,:),'-', 'color', colours(i,:))
end
hold off
%}
