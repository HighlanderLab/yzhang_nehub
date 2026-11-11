
import json
import numpy as np
import re
def average_time_size(filename,sampling_number_limit=np.inf):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    l_min_size = -1
    count = 0
    for pop, samplings in data.items():
        times = []
        sizes = []
        for  sampling,sampling_Ne in samplings.items():
            if type(sampling_Ne) == dict:
                for time, size in sampling_Ne.items():
                    times.append(float(time))
                    sizes.append(float(size))
            else:
                times.append(float(sampling))
                sizes.append(float(sampling_Ne))
        if l_min_size == -1:
            l_min_size = len(sizes)
        else:
            l_min_size = min(l_min_size,  len(sizes))
        if count >= sampling_number_limit: #because some msmc run >20 times, I only want take these 20
            break
        all_times.append(times)
        all_sizes.append(sizes)
        count += 1
    print(len(all_times))
    # Calculate the average size for each time point
    np_times = np.zeros((len(all_times), l_min_size))
    np_sizes = np.zeros((len(all_sizes), l_min_size))
    for i, sizes in enumerate(all_sizes):
        np_sizes[i, :] = sizes[:l_min_size]
    for i, times in enumerate(all_times):
        np_times[i, :] = times[:l_min_size]
    
    average_times = np.mean(np_times, axis=0)
    average_sizes = np.mean(np_sizes, axis=0)
    median_sizes = np.median(np_sizes, axis=0)
    return average_times,average_sizes,median_sizes,all_times,all_sizes
    
def pop_average_time_size(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    all_median_sizes = []
    all_average_sizes = []
    all_average_times = []
    pops = []
    for pop, pop_results in data.items():
        pop_times = []
        pop_sizes = []
        
        l_min_size = -1
        for individual, Ne in pop_results.items():
            times = []
            sizes = []
            for time, size in Ne.items():
                times.append(float(time))
                sizes.append(float(size))
            if len(sizes) == 1: # this individual get NAN
                pass
            else:
                if l_min_size == -1:
                    l_min_size = len(sizes)
                else:
                    l_min_size = min(l_min_size,  len(sizes))
                pop_times.append(times)
                pop_sizes.append(sizes)
        # Calculate the average size for each time point
        np_times = np.zeros((len(pop_times), l_min_size))
        np_sizes = np.zeros((len(pop_sizes), l_min_size))
        for i, every_sizes in enumerate(pop_sizes):
            np_sizes[i, :] = every_sizes[:l_min_size]
        for i, every_times in enumerate(pop_times):
            np_times[i, :] = every_times[:l_min_size]
        
        average_times = np.mean(np_times, axis=0)
        average_sizes = np.mean(np_sizes, axis=0)
        median_sizes = np.median(np_sizes, axis=0)

        all_average_times.append(average_times)
        all_average_sizes.append(average_sizes)
        all_median_sizes.append(median_sizes)
        all_times.append(pop_times)
        all_sizes.append(pop_sizes)
        pops.append(pop)

    return pops,all_average_times, all_average_sizes, all_median_sizes,all_times,all_sizes


def pop_average_sampling_time_size(filename,sampling_number_limit=np.inf):
    with open(filename, "r") as f:
        data = json.load(f)
    time_pops = {}
    size_pops = {}
    for sampling_results in data:
        for pop, pop_results in sampling_results.items():
            times = []
            sizes = []
            for time, size in pop_results.items():
                times.append(float(time))
                sizes.append(float(size))
            pop_group = re.sub(r'_\d+\.final\.txt', '', pop)
            if pop_group in time_pops:
                time_pops[pop_group].append(times)
                size_pops[pop_group].append(sizes)
            else:
                time_pops[pop_group] = [times]
                size_pops[pop_group] = [sizes]
    all_times = []
    all_sizes = []
    all_median_sizes = []
    all_average_sizes = []
    all_average_times = []
    pops = []
    for pop_group in time_pops.keys():
        pop_times = time_pops[pop_group]
        pop_sizes = size_pops[pop_group]
        if len(pop_times) > sampling_number_limit:
            pop_times[:sampling_number_limit]
            pop_sizes[:sampling_number_limit]
        l_min_size = -1
        for sizes in pop_sizes:
            if l_min_size == -1:
                l_min_size = len(sizes)
            else:
                l_min_size = min(l_min_size,  len(sizes))
        
        # Calculate the average size for each time point
        np_times = np.zeros((len(pop_times), l_min_size))
        np_sizes = np.zeros((len(pop_sizes), l_min_size))
        for i, every_sizes in enumerate(pop_sizes):
            np_sizes[i, :] = every_sizes[:l_min_size]
        for i, every_times in enumerate(pop_times):
            np_times[i, :] = every_times[:l_min_size]
        
        average_times = np.mean(np_times, axis=0)
        average_sizes = np.mean(np_sizes, axis=0)
        median_sizes = np.median(np_sizes, axis=0)

        all_average_times.append(average_times)
        all_average_sizes.append(average_sizes)
        all_median_sizes.append(median_sizes)
        all_times.append(pop_times)
        all_sizes.append(pop_sizes)
        pops.append(pop_group)

    return pops,all_average_times, all_average_sizes, all_median_sizes,all_times,all_sizes

def pop_time_size(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    l_min_size = -1
    for pop, pop_results in data.items():
        times = []
        sizes = []
        for  sampling,sampling_Ne in pop_results.items():
            if type(sampling_Ne) == dict:
                for time, size in sampling_Ne.items():
                    times.append(float(time))
                    sizes.append(float(size))
            else:
                times.append(float(sampling))
                sizes.append(float(sampling_Ne))
        pop_list.append(pop)
        all_times.append(times)
        all_sizes.append(sizes)
    return pop_list,all_times,all_sizes


def pop_time_size_tskit(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    for pop, pop_results in data.items():
        times = []
        sizes = []
        Weighted_coalrate = pop_results["Weighted_coalrate"]
        for time, size in Weighted_coalrate.items():
            times.append(float(time))
            if float(size) > 1e+8:
                sizes.append(np.inf)
            else:
                sizes.append(float(size))
        pop_list.append(pop)
        all_times.append(times)
        all_sizes.append(sizes)
    print(len(all_times))
    return pop_list,all_times,all_sizes

def pop_time_size_Relate(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    for file_name, file_results in data.items():
        for pop, pop_results in file_results.items():
            times = []
            sizes = []
            if pop.split(" - ")[0] == pop.split(" - ")[1]:
                for time, size in pop_results.items():
                    if float(size) > 1e+8:
                        break
                    else:                       
                        times.append(float(time))
                        sizes.append(float(size))
                pop_list.append(pop.replace(" - ","_"))
                all_times.append(times)
                all_sizes.append(sizes)
    return pop_list,all_times,all_sizes

def pop_time_size_stairway(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    for pop, pop_results in data.items():
        times = []
        median_sizes = []
        lower_sizes_2_5 = []
        upper_sizes_97_5 = []
        lower_sizes_12_5 = []
        upper_sizes_87_5 = []
        for  p, v in pop_results.items():
            times.append(float(p))
            median_sizes.append(float(v[0]))
            lower_sizes_2_5.append(float(v[1]))
            upper_sizes_97_5.append(float(v[2]))
            lower_sizes_12_5.append(float(v[3]))
            upper_sizes_87_5.append(float(v[4]))
        all_times.append(times)
        all_sizes.append([median_sizes,lower_sizes_2_5,upper_sizes_97_5,lower_sizes_12_5,upper_sizes_87_5])
        pop_list.append(pop)
    return pop_list,all_times,all_sizes

def pop_time_size_Phlash(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    for pop, pop_results in data.items():
        times = []
        mean_sizes = []
        median_sizes = []
        lower_sizes_2_5 = []
        upper_sizes_97_5 = []
        lower_sizes_12_5 = []
        upper_sizes_87_5 = []
        for  p, v in pop_results.items():
            times.append(float(p))
            mean_sizes.append(float(v[0]))
            median_sizes.append(float(v[1]))
            lower_sizes_2_5.append(float(v[2]))
            upper_sizes_97_5.append(float(v[3]))
            lower_sizes_12_5.append(float(v[4]))
            upper_sizes_87_5.append(float(v[5]))
        pop_list.append(pop)
        all_times.append(times)
        all_sizes.append([median_sizes,lower_sizes_2_5,upper_sizes_97_5,lower_sizes_12_5,upper_sizes_87_5])
    return pop_list,all_times,all_sizes

def pop_time_size_SMCpp(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    all_times = []
    all_sizes = []
    pop_list = []
    for file_name, file_results in data.items():
        for pop, pop_results in file_results.items():
            times = []
            sizes = []
            
            for time, size in pop_results.items():
                times.append(float(time))
                sizes.append(float(size))
            pop_list.append(file_name+pop)
            all_times.append(times)
            all_sizes.append(sizes)
    return pop_list,all_times,all_sizes