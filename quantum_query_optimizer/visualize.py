import matplotlib.pyplot as plt
import numpy as np

def visualize(x_var, y_var, xlabel, ylabel, title, linestyle="-"):
    '''
        Parameters:
            x_var : x-axis data
            y_var : y-axis data
            xlabel : x-axis label
            ylabel : y-axis label
            title : figure title
        Returns:
            plt : plot created
    '''
    plt.rcParams['axes.facecolor'] = 'w'
    plt.plot(x_var, y_var, 'k', linestyle=linestyle)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(range(min(x_var), max(x_var)+1))
    return plt

def visualizeRuntime(solutions, title="Runtime by Input Size", filename="RuntimeByInputSize.eps"):
    '''
        Parameters:
            solutions : solutions dictionary from SDP solution
            title : plot title
            filename : output file name

        This function visualizes runtime data
    '''
    visualize(
        x_var=solutions['n_bitstring'],
        y_var=solutions['run_time'],
        xlabel='Input Size (n)',
        ylabel='Runtime (seconds)',
        title=title
    )
    plt.savefig(filename)
    plt.close() 

def visualizeComplexity(solutions, title="Query Complexity by Input Size", filename="ComplexityByInputSize.eps"):
    '''
        Parameters:
            solutions : solutions dictionary from SDP solution
            title : plot title
            filename : output file name

        This function visualizes quantum query complexity data
    '''
    visualize(
        x_var=solutions['n_bitstring'],
        y_var=solutions['query_complexity'],
        xlabel='Input Size (n)',
        ylabel='Queries',
        title=title
    )
    plt.savefig(filename)
    plt.close() 

def visualizeComplexityOR(solutions, title, filename):
    xs = [int(n) for n in solutions['n_bitstring']]
    ys = np.sqrt(xs)
    plt.plot(xs, ys, color="red")
    visualize(
        x_var=solutions['n_bitstring'],
        y_var=solutions['query_complexity'],
        xlabel='Input Size (n)',
        ylabel='Queries',
        title=title
    )
    plt.legend(['Analytical', 'Empirical'])
    plt.savefig(filename)
    plt.close()  

def visualizeRuntimeOR(all_solutions, worst_solutions, title, filename):
    visualize(
        x_var=all_solutions['n_bitstring'],
        y_var=all_solutions['run_time'],
        xlabel='Input Size (n)',
        ylabel='Runtime (seconds)',
        title=title,
        linestyle="-"
    )    
    visualize(
        x_var=worst_solutions['n_bitstring'],
        y_var=worst_solutions['run_time'],
        xlabel='Input Size (n)',
        ylabel='Runtime (seconds)',
        title=title,
        linestyle="--"
    )
    plt.legend(['All Inputs', 'Worst-case Inputs'])
    plt.savefig(filename)
    plt.close()  
