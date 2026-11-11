
import json
import numpy as np
def plot_multi_results(ax, average_times, average_sizes,median_sizes, all_times, all_sizes,method,color,color_average,shape,shape_average="--",max_time=1e+5,label=""):
    # "./con/new/GONE_results.json"
    # method='GONE'
    #color=palette[0]
    # shape='-'
    if not np.any(average_times) or not np.any(average_sizes):
        print(f"Skipping {method}, no valid data.")
        return

    # Find the cutoff index where time exceeds 1e5
    indices = [index for index, size in enumerate(average_times) if size > max_time]
    index_max_time = indices[0] if len(indices)>0 else len(average_times)
    print(f"Index max time: {index_max_time}")
    if label!="":
        label = label+"_"
    # Plot individual traces with low opacity
    for times, sizes in zip(all_times, all_sizes):
        if times[0] ==0:
            ax.step(times[:index_max_time], 
                sizes[:index_max_time], 
                color=color, linestyle=shape, alpha=0.2,where='post')
        else:
            ax.plot(times[:index_max_time], 
                sizes[:index_max_time], 
                color=color, linestyle=shape, alpha=0.2)
    if average_times[0] ==0:
        ax.step(average_times[:index_max_time], 
                average_sizes[:index_max_time], 
                color=color_average, linestyle=shape_average, label=label+'Average',where='post')
        ax.step(average_times[:index_max_time], 
            median_sizes[:index_max_time], 
            color=color_average, linestyle=shape, label=label+'Median',where='post')
    else:
        ax.plot(average_times[:index_max_time], 
            average_sizes[:index_max_time], 
            color=color_average, linestyle=shape_average, label=label+'Average')
        ax.plot(average_times[:index_max_time], 
            median_sizes[:index_max_time], 
            color=color_average, linestyle=shape, label=label+'Median')

def plot_single_results(ax, times, sizes,label,color,shape,max_time=1e+5,alpha=1):
    # "./con/new/GONE_results.json"
    # method='GONE'
    #color=palette[0]
    # shape='-'
    # Find the cutoff index where time exceeds 1e5
    indices = [index for index, size in enumerate(times) if size > max_time]
    index_max_time = indices[0] if len(indices)>0 else len(times)
    print(f"Index max time: {index_max_time}")

    # Plot individual traces with low opacity
    if times[0] ==0:
        ax.step(times[:index_max_time], 
                sizes[:index_max_time], 
                color=color, linestyle=shape, label=label, alpha=alpha, where='post')
    else:
        ax.plot(times[:index_max_time], 
                sizes[:index_max_time], 
                color=color, linestyle=shape, label=label, alpha=alpha)

def plot_single_results_with_CI(ax, times, sizes,label,color_CI,color,shape,max_time=1e+5,alpha=1):
    # "./con/new/GONE_results.json"
    # method='GONE'
    #color=palette[0]
    # shape='-'
    # Find the cutoff index where time exceeds 1e5
    indices = [index for index, size in enumerate(times) if size > max_time]
    index_max_time = indices[0] if len(indices)>0 else len(times)
    print(f"Index max time: {index_max_time}")

    # Plot individual traces with low opacity
    if len(sizes)==6:
        ax.plot(times[:index_max_time], 
                sizes[0][:index_max_time], 
                color=color, linestyle='--', label='Average')
        ax.plot(times[:index_max_time], 
                sizes[1][:index_max_time], 
                color=color, linestyle=shape, label='Median')
        ax.fill_between(times[:index_max_time],sizes[2][:index_max_time] , sizes[3][:index_max_time] ,color=color_CI, alpha=.15, linewidth=0)
        ax.fill_between(times[:index_max_time],sizes[4][:index_max_time] , sizes[5][:index_max_time], color=color_CI,alpha=.3, linewidth=0)

    else:
        ax.plot(times[:index_max_time], 
                sizes[0][:index_max_time], 
                color=color, linestyle=shape, label='Median')
        ax.fill_between(times[:index_max_time],sizes[1][:index_max_time] , sizes[2][:index_max_time] ,color=color_CI, alpha=.15, linewidth=0)
        ax.fill_between(times[:index_max_time],sizes[3][:index_max_time] , sizes[4][:index_max_time], color=color_CI,alpha=.3, linewidth=0)

def plot_single_results_without_CI(ax, times, sizes,label,color,shape,max_time=1e+5,alpha=1):
    # "./con/new/GONE_results.json"
    # method='GONE'
    #color=palette[0]
    # shape='-'
    # Find the cutoff index where time exceeds 1e5
    indices = [index for index, size in enumerate(times) if size > max_time]
    index_max_time = indices[0] if len(indices)>0 else len(times)
    print(f"Index max time: {index_max_time}")

    # Plot individual traces with low opacity
    if len(sizes)==6:
        ax.plot(times[:index_max_time], 
                sizes[1][:index_max_time], 
                color=color, linestyle=shape, label=label)
    else:
        ax.plot(times[:index_max_time], 
                sizes[0][:index_max_time], 
                color=color, linestyle=shape, label=label)

def plot_single_results_only_CI(ax, times, sizes,label,color,shape,max_time=1e+5,alpha=1):
    # "./con/new/GONE_results.json"
    # method='GONE'
    #color=palette[0]
    # shape='-'
    # Find the cutoff index where time exceeds 1e5
    indices = [index for index, size in enumerate(times) if size > max_time]
    index_max_time = indices[0] if len(indices)>0 else len(times)
    print(f"Index max time: {index_max_time}")

    # Plot individual traces with low opacity
    if len(sizes)==6:
        ax.fill_between(times[:index_max_time],sizes[2][:index_max_time] , sizes[3][:index_max_time] ,color=color, alpha=alpha/7, linewidth=0)
        ax.fill_between(times[:index_max_time],sizes[4][:index_max_time] , sizes[5][:index_max_time], color=color,alpha=alpha/3, linewidth=0)

    else:
        ax.fill_between(times[:index_max_time],sizes[1][:index_max_time] , sizes[2][:index_max_time] ,color=color,alpha=alpha/7, linewidth=0)
        ax.fill_between(times[:index_max_time],sizes[3][:index_max_time] , sizes[4][:index_max_time], color=color,alpha=alpha/3, linewidth=0)