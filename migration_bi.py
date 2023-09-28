import numpy as np
from matplotlib import pyplot as plt

bi_num = 1 # for selecting particle sizes, defined in rs, rl
col_num = 1 # colour for plotting
case_name = 'bi2_phi20'
case_range = range(1,1+1)
cases = []
for c in case_range:
    cases.append(f'{case_name}_{c}')
directory = 'bi'
file_name = 'particles'

# set timestep variables for looping
rate = 10000
stepNum = range(10000,100000,rate)
plot_interval = 28

dt_LBM = 2.6666667e-8
dt_DEM = 2.6666667e-8
dx = 0.4e-3 # lattice spacing
a_dx = 1.1e6 # acceleration gradient for effective viscosity calc, mm-g (mm/s^2) = pressure gradient m-kg (kg/m^2s^2)
nu = 1/6*dx**2./dt_LBM
aw =  1.0e-3 # wall particle radius
aw_factor = 2
X = 64e-3
Y = 38e-3
Z = 104.8e-3
W = Z - 2*dx - 2*aw_factor*aw

# set up spatial bins for measuring particle quantities
N_spat_bins = 19
bins = np.linspace(dx+aw_factor*aw+0.00001,Z-dx-aw_factor*aw-0.00001,N_spat_bins+1)
bin_width = W/N_spat_bins
bin_vol = bin_width*X*Y
total_vol = bin_vol*N_spat_bins
bin_centres = np.linspace(dx+aw_factor*aw+bin_width/2,Z-dx-aw_factor*aw-bin_width/2,N_spat_bins)
bin_centres_norm = (bin_centres-dx-aw_factor*aw)/((max(bins)-dx-aw_factor*aw))

# possible particle sizes (depending on case, selected by bi_num)
rs = [1.3e-3, 1.47e-3, 1.83e-3, 1.405e-3, 1.53e-3]
rl = [2.7e-3, 2.88e-3, 3.23e-3, 3.05e-3, 2.66e-3]
phi_max = [0.783,0.7395,0.6804,0.783,0.7395,0.6804] # pre-calculate and hard code these based on particle sizes
species_sizes = np.array([rs[bi_num], rl[bi_num]])

# initialise plots
markers = ['o','s','^','d']
colours = ['b','r','g','o']
fig, ax = plt.subplots(2,2)
ax[0,0].set_xlabel(r'd ($\mu m$)')
ax[0,0].set_ylabel('$C_j$')
ax[0,1].set_xlabel('$L_x/W$')
ax[0,1].set_ylabel('C')
ax[1,0].set_xlabel('$z/W$')
ax[1,0].set_ylabel(r'$u/u_{max}$')
ax[1,1].set_xlabel('$z/W$')
ax[1,1].set_ylabel('$\phi$')

# initialise arrays for storing data
Lx_median = np.zeros((len(stepNum),2,len(cases)))
Lx_Q1 = np.zeros((len(stepNum),2,len(cases)))
Lx_Q3 = np.zeros((len(stepNum),2,len(cases)))
C_seg = np.zeros((len(cases),len(stepNum)))
Cj = np.zeros((len(stepNum),2,len(cases)))
C = np.zeros((len(stepNum),len(cases)))
vel_profile_norm = np.zeros((N_spat_bins,len(cases)))
bin_phi_part = np.zeros((2,N_spat_bins,len(cases)))
bin_phi_total = np.zeros((N_spat_bins,len(cases)))

for c in range(len(cases)):
    case_name = cases[c]
    # initialise lists for particle counters, must be within case loop
    bin_count = np.zeros((len(stepNum),2,N_spat_bins))
    bin_count_vol = np.zeros((len(stepNum),2,N_spat_bins))
    bin_vel_x = np.zeros((len(stepNum),2,N_spat_bins))
    bin_vel_y = np.zeros((len(stepNum),2,N_spat_bins))
    inner_Cj = np.zeros((len(stepNum),2)) # for storing the inner sum of scalar dispersion function (C) for each part size at each timestep
    median_vels = np.zeros((len(stepNum),2))
    Q1_vels = np.zeros((len(stepNum),2))
    Q3_vels = np.zeros((len(stepNum),2))

    # count actual numbers of particles @ t=0
    with open(f'{directory}/{case_name}/{file_name}_1') as f:
        lines = f.readlines()[9:]
        count_j = np.zeros(2)
        for j in range(2):
            for p in range(len(lines)):
                data = lines[p].split()
                if int(data[1]) == 1 and float(data[11]) == species_sizes[j]:
                    count_j[j] += 1

    for t in range(len(stepNum)):
        print(f'{t} / {len(stepNum)}')
        velY_total = 0
        # read in data from LIGGGHTS dump file
        with open(f'{directory}/{case_name}/{file_name}_{stepNum[t]}') as f:
            lines = f.readlines()[9:]

            # median particle velocity of each size
            for r in range(2):
                part_vels = []
                for p in range(len(lines)):
                    data = lines[p].split()
                    p_rad = float(data[11])
                    p_vel_x = float(data[5])
                    if p_rad == species_sizes[r] and p_vel_x != 0:
                        part_vels.append(p_vel_x)
                median_vels[t,r] = np.median(part_vels)
                Q1_vels[t,r] = np.quantile(part_vels,0.25)
                Q3_vels[t,r] = np.quantile(part_vels,0.75)

            # iterate through particles
            for p in range(len(lines)):
                data = lines[p].split()
                p_pos = float(data[4])
                p_rad = float(data[11])
                p_vel_x = float(data[5])
                p_vel_y = float(data[6])
                # calculate inner scalar dispersion function of current particle j, and sum with appropriate particle size
                for k in range(2):
                    if p_rad == species_sizes[k]:
                        if p_pos > bins[0] and p_pos < bins[-1]:
                            inner_Cj[t,k] += abs(p_pos/W - 0.5)
                # iterate through bins
                for j in range(N_spat_bins):
                    # set lower and upper threshold of current bin
                    bin_min = bins[j]
                    bin_max = bins[j+1]
                    # check if current particle within current part bin (bin k)
                    if p_pos > bin_min and p_pos < bin_max:
                        # calculate particle velocities in each bin
                        for k in range(2):
                            if p_rad == species_sizes[k]:
                                bin_vel_x[t,k,j] += p_vel_x
                                bin_vel_y[t,k,j] += p_vel_y
                                bin_count[t,k,j] += 1
                        # spatial bin j+1
                        vol_j1 = 0
                        # ensure bin j isn't last bin
                        if j < N_spat_bins - 1:
                            # calculate particle overlap into bin j+1 (h_j1)
                            h_j1 = p_pos + p_rad - bin_max
                            # estimate volume of current particle occupying bin j+1
                            if h_j1 > 0:
                                vol_j1 = np.pi*h_j1**3*(p_rad/h_j1-1/3)
                            for k in range(2):
                                if p_rad == species_sizes[k]:
                                    bin_count_vol[t,k,j+1] += vol_j1
                        # spatial bin j-1
                        vol_j_1 = 0
                        # ensure bin j isn't first bin
                        if j > 0:
                            # calculate particle overlap into bin j-1 (h_j_1)
                            h_j_1 = bin_min - (p_pos - p_rad)
                            # estimate volume of current particle occupying bin j-1
                            if h_j_1 > 0:
                                vol_j_1 = np.pi*h_j_1**3*(p_rad/h_j_1-1/3)
                            for k in range(2):
                                if p_rad == species_sizes[k]:
                                    bin_count_vol[t,k,j-1] += vol_j_1
                        # calculate volume of current particle in current bin
                        for k in range(2):
                            if p_rad == species_sizes[k]:
                                bin_count_vol[t,k,j] += 4/3*np.pi*p_rad**3 - vol_j1 - vol_j_1

        part_total = sum(sum(bin_count[t,:,:]))

    # calculate mean bin and channel velocities, then mean channel Re and Q
    bin_vel_x_mean = bin_vel_x/bin_count
    channel_vel_mean = np.nanmean(np.nanmean(bin_vel_x_mean,axis=2),axis=1)
    vel_profile_norm[:,c] = np.nanmean(bin_vel_x_mean[-1,:,:],axis=0)/channel_vel_mean[-1]
    Re_c = channel_vel_mean*W/nu
    Q = channel_vel_mean*X*W
    Re_c_end = Re_c[-1]
    Q_end_m3s = Q[-1]/1000**3/(X/1000)

    nu_eff = a_dx*W**2/12./(channel_vel_mean) # mm^2/s
    Re_c_eff = channel_vel_mean*W/nu_eff
    Re_c_eff_alt = 12*channel_vel_mean**2/a_dx/W

    # calculate number of particles in each bin
    part_count = np.sum(bin_count[-1,:,:],axis=1)
    # calculate phi of each size species j
    phi_j = count_j*(4/3*np.pi*np.power(species_sizes,3))/total_vol
    phi_total = sum(phi_j)

    # phi for each spatial bin to plot phi distribution
    bin_phi_part[:,:,c] = np.squeeze(bin_count_vol[-1,:,:])/bin_vol
    bin_phi_total[:,c] = np.squeeze(np.sum(bin_phi_part,axis=0))

    # calculate characteristic length particles have travelled along channel at each timestep,
    # defined as in Chun et al. (2019)
    part_mean_vel_x =  np.sum(bin_vel_x,axis=2)/part_count
    part_mean_vel_y =  np.sum(bin_vel_y,axis=2)/part_count
    Lx = np.zeros((len(stepNum),1))
    Lx_part = np.zeros((len(stepNum),2))
    Ly_part = np.zeros((len(stepNum),2))
    
    Lx = np.cumsum(channel_vel_mean)*rate*dt_DEM
    Lx_part = np.cumsum(part_mean_vel_x,axis=0)*rate*dt_DEM
    Ly_part = np.cumsum(part_mean_vel_y,axis=0)*rate*dt_DEM
    Lx_median[:,:,c] = np.cumsum(median_vels,axis=0)*rate*dt_DEM
    Lx_Q1[:,:,c] = np.cumsum(Q1_vels,axis=0)*rate*dt_DEM
    Lx_Q3[:,:,c] = np.cumsum(Q3_vels,axis=0)*rate*dt_DEM

    Lx_W = Lx/W
    Lx_W_part = Lx_part/W
    Ly_mean = np.mean(Ly_part,axis=1)
    Ly_Lx = Ly_part/Lx_part
    Lx_W_part_norm = Lx_W_part[-1,:]/np.max(Lx_W_part[-1,:])

    # scalar dispersion functions
    # calculate scalar dispersion function for each particle size species
    Cj_current = inner_Cj/np.transpose(count_j)
    Cj[:,:,c] = Cj_current
    # total scaler dispersion function
    C[:,c] = np.sum(Cj_current*np.transpose(phi_j),axis=1)/phi_total
    # scalar segregation function
    C_seg_inner = np.zeros((len(stepNum),2))
    for i in range(2):
        for j in range(2):
            C_seg_inner[:,i] = np.abs(Cj_current[:,i]-Cj_current[:,j]) + C_seg_inner[:,i]
    C_seg_inner /= 2
    C_seg[c,:] = 4*np.sum(C_seg_inner,axis=1)/2

bin_phi_part_mean = np.mean(bin_phi_part,axis=2)
bin_phi_total_mean = np.mean(bin_phi_total,axis=1)
vel_profile_mean = np.mean(vel_profile_norm,axis=1)

Cj_median = np.median(Cj,axis=2)
Cj_mean = np.mean(Cj,axis=2)
C_median = np.median(C,axis=1)
C_mean = np.mean(C,axis=1)
C_sigma = (1-4*C_mean)/(1-phi_total/phi_max[bi_num])
C_seg_median = np.median(C_seg,axis=0)

ax[0,0].plot(2*1000*species_sizes[:], Cj_median[-1,:], f'{markers[col_num]}{colours[col_num]}')

ax[0,1].plot(Lx_W[:-1:plot_interval], C_mean[:-1:plot_interval], f'{markers[col_num]}{colours[col_num]}')
lines = ['-','--']
for i in range(2):
    ax[0,1].plot(Lx_W[:-1:plot_interval], Cj_median[:-1:plot_interval,i], f'{lines[i]}{colours[col_num]}')

ax[1,0].plot(bin_centres_norm, vel_profile_mean, f'{markers[col_num]}{colours[col_num]}')

ax[1,1].plot(bin_centres_norm, bin_phi_total_mean, f'{markers[col_num]}{colours[col_num]}')
#ax[1,1].plot(Lx_W[:-1:plot_interval], C_seg_median[:-1:plot_interval], f'{markers[col_num]}{colours[col_num]}')

plt.show()