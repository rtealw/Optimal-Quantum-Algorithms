import matplotlib.pyplot as plt
import numpy as np

def visualize(x_var, y_var, xlabel, ylabel, title):
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
    plt.plot(x_var, y_var, 'k', label="Empirical", color="red")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    return plt

def visualizeRuntime(solutions, title="Runtime by Input Size", filename="RuntimeByInputSize.eps"):
    '''
        Parameters:
            solutions : solutions dictionary from SDP solution
            title : plot title
            filename : output file name

        This function visualizes runtime data
    '''
    plt = visualize(
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
    plt = visualize(
        x_var=solutions['n_bitstring'],
        y_var=solutions['query_complexity'],
        xlabel='Input Size (n)',
        ylabel='Queries',
        title=title
    )
    ns = [int(n) for n in solutions['n_bitstring']]
    x = np.arange(min(ns), max(ns)+1)
    y = np.sqrt(x)
    plt.plot(x,y, label="Theoretical", color="black")
    plt.xticks(x)
    plt.legend(['Theoretical', 'Empirical'])
    plt.savefig(filename)
    plt.close() 
