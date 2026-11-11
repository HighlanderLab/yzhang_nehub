import numpy as np
from scipy.interpolate import interp1d

def make_time_grid(t_min, t_max, n=200):
    return np.geomspace(t_min, t_max, n)

def interp_timepoint(times, Ne):
    """
    times, Ne: 1D arrays
    returns callable Ne(t)
    """
    order = np.argsort(times)
    f = interp1d(
        times[order],
        Ne[order],
        kind="linear", # linear function
        bounds_error=False,
        fill_value=(np.nan,np.nan)  # if out of bounds
    )
    return f

def interp_timeperiod(t_lower, t_upper, Ne):
    """
    t_lower, t_upper, Ne: 1D arrays of equal length
    returns callable Ne(t)
    """
    def f(t):
        t = np.asarray(t)
        Ne_t = np.full_like(t, np.nan, dtype=float)
        for lo, hi, val in zip(t_lower, t_upper, Ne):
            mask = (t >= lo) & (t < hi)
            Ne_t[mask] = val
        return Ne_t
    return f

def rmse_log10(Ne_est, Ne_true):
    Ne_est = np.asarray(Ne_est)
    Ne_true = np.asarray(Ne_true)

    # valid positive values only
    mask = (Ne_est > 0) & (Ne_true > 0)

    log_est = np.log10(Ne_est[mask])
    log_true = np.log10(Ne_true[mask])

    # remove nan and inf after log
    finite = np.isfinite(log_est) & np.isfinite(log_true)

    return np.sqrt(np.mean((log_est[finite] - log_true[finite])**2))

def rmse_by_window(t, Ne_est, Ne_true, t0, t1):
    mask = (t >= t0) & (t <= t1)
    return rmse_log10(Ne_est[mask], Ne_true[mask])

def cal_rmse(times,sizes,true_time,true_size,split_epoch,**kwargs):
    '''To enable fair comparison between methods that infer effective population size at discrete time points 
    and those that infer piecewise-constant trajectories, we interpolated all inferred Ne trajectories onto a common log-spaced time grid. 
    RMSE was then computed on the log10 Ne scale across the shared time domain.'''

    true_time_lower=true_time[:-1]
    true_time_upper=true_time[1:]
    Ne_true=true_size[:-1]
    f_true = interp_timeperiod(true_time_lower, true_time_upper, Ne_true)
    
    if (times[0] == 0) or (times[0] == "0"): # Time-period
        print("handle time-period")
        time_lower = times[:-1]
        time_upper = times[1:]
        Ne_size = sizes[:-1]
        f_time = interp_timeperiod(time_lower, time_upper, Ne_size)
    else: # time point
        print("handle time-point")
        f_time = interp_timepoint(np.array(times), np.array(sizes))
    t_grid = make_time_grid(split_epoch[0]+1, split_epoch[-1],**kwargs)
    t_grid = t_grid[(t_grid >= times[0]) & (t_grid <= times[-1])]
    Ne_true_g = f_true(t_grid)
    Ne_est_g = f_time(t_grid)
    rmse_epoch=[]
    for i in range(len(split_epoch)-1):
        rmse_epoch.append(rmse_by_window(t_grid, Ne_est_g, Ne_true_g, split_epoch[i], split_epoch[i+1]))
    return rmse_epoch,rmse_log10(Ne_est_g, Ne_true_g)



def cal_rmse_old(times,sizes,true_time,true_size,split_epoch=0):
    list_epoch_rmse=[]
    if times[0] != 0:
        times.insert(0,0)
        sizes.insert(0,sizes[0])
    if times[-1] < true_time[-1]:
        indics = np.searchsorted(true_time, times[-1])
        true_time = true_time[:indics]
        true_time.append(times[-1])
    if split_epoch ==0:
        combined_times = true_time
        aligned_sizes = true_size
    else:
        combined_times = np.sort(set(np.concatenate([true_time, split_epoch])))
        # Find which split_epoch values are NOT in true_times
        only_split_epoch = np.setdiff1d(split_epoch, true_time)
        # Find insertion positions of combined_times in true_times
        indices = np.searchsorted(true_time, combined_times)
        # Create a boolean mask for split_epoch-only times
        is_split_only = np.isin(combined_times, only_split_epoch)
        # Adjust indices: use i-1 for split-only times, i for others
        adjusted_indices = np.where(is_split_only, indices - 1, indices)
        # Get aligned sizes in one operation
        aligned_sizes = np.array(true_size)[adjusted_indices]

    indices = np.searchsorted(times, combined_times)
    split_times = np.split(times, indices)
    split_sizes = np.split(sizes, indices)
    total_sum_square_size=[]
    total_time = []
    for i in range(1,len(combined_times)):
        if len(split_times[i]) != 0:
            time_part = list(split_times[i])
            size_part = list(split_sizes[i])
            if time_part[0] > combined_times[i-1]:
                time_part.insert(0, combined_times[i-1])
                size_part.insert(0, last_split_size)
            if time_part[-1] < combined_times[i]:
                time_part.append(combined_times[i])
            last_split_size=split_sizes[i][-1]
        else:
            time_part = [combined_times[i-1],combined_times[i]]
            size_part = [last_split_size]
        
        sum_square_size = []
        epoch_time = []
        for epoch_index in range(len(size_part)):
            if size_part[epoch_index] != np.inf: 
                sum_square_size.append((size_part[epoch_index]-true_size[i-1])**2  * (time_part[epoch_index+1] - time_part[epoch_index]))
                epoch_time.append(time_part[epoch_index+1] - time_part[epoch_index])
        
        if len(epoch_time) !=0:  
            list_epoch_rmse.append([true_time[i-1],true_time[i], np.sqrt(sum(sum_square_size)/sum(epoch_time))])
            total_sum_square_size.extend(sum_square_size)
            total_time.extend(epoch_time)
        
    total_rmse = np.sqrt(sum(total_sum_square_size)/sum(total_time))
    return list_epoch_rmse,total_rmse 