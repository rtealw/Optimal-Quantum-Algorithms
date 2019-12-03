import matplotlib.pyplot as plt

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
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = 'w'
    plt.plot(x_var, y_var, 'k')
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
    plt.savefig(filename)
    plt.close() 
