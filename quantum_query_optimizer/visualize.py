import matplotlib.pyplot as plt
import numpy as np

def setup(xvar, xlabel, ylabel, title):
    '''
        Parameters:
            xvar : integer list of x-values
            xlabel : x-axis label
            ylabel : y-axis label
            title : plot title

        This function sets up the figure
    '''
    plt.rcParams['axes.facecolor'] = 'w'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(xvar)

def visualizeRuntime(solutions, title="Runtime by Input Size", filename=None, labels=[], more_runtimes=[]):
    '''
        Parameters:
            solutions : solutions dictionary from SDP solution
            title : plot title
            filename : output file name
            labels (optional) : labels for runtimes
            more_runtimes (optional) : more runtimes to visualize on plot

        This function visualizes runtime data
    '''
    xvar = [int(n) for n in solutions['n_bitstring']]
    setup(xvar=xvar, xlabel='Runtime (seconds)', ylabel='Queries', title=title) 
    label_yes = len(labels) == 1 + len(more_runtimes)
    label_num = 0
    for runtimes in [solutions['run_time']] + more_runtimes:
        label = "" if not label_yes else labels[label_num]
        plt.plot(xvar, runtimes, label=label)
        label_num += 1
    if label_yes:
        plt.legend()
    if filename == None:
        plt.show()
    else:
        plt.savefig(filename)
    plt.close() 

def visualizeComplexity(solutions, title="Query Complexity by Input Size", filename=None, labels=[], functions=[], more_queries=[]):
    '''
        Parameters:
            solutions : solutions dictionary from SDP solution
            title : plot title
            filename : output file name
            labels (optional) : labels for sets of queries and functions
            functions (optional) : functions to plot
            more_queries (optional) : more sets of queries to plot

        This function visualizes quantum query complexity data
    '''
    xvar = [int(n) for n in solutions['n_bitstring']]
    setup(xvar=xvar, xlabel='Input Size (n)', ylabel='Queries', title=title)
    label_no = len(labels) != 1 + len(more_queries) + len(functions)
    label_num = 0
    for queries in [solutions['query_complexity']] + more_queries:
        label = "" if label_no else labels[label_num]
        plt.plot(xvar, queries, label=label)
        label_num += 1
    for function in functions:
        label = "" if label_no else labels[label_num]
        plt.plot(xvar, function(xvar), label=label)
        label_num += 1
    if not label_no:
        plt.legend()
    if filename==None:
        plt.show()
    else:
        plt.savefig(filename)
    plt.close() 
