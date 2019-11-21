import matplotlib.pyplot as plt

def visualize(x_var, y_var, xlabel, ylabel, title):
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = 'w'
    plt.plot(x_var, y_var, 'k')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

def visualizeRuntime(solutions, title="Runtime by Inpuz Size", filename="RuntimeByInputSize.eps"):
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
    visualize(
        x_var=solutions['n_bitstring'],
        y_var=solutions['query_complexity'],
        xlabel='Input Size (n)',
        ylabel='Queries',
        title=title
    )
    plt.savefig(filename)
    plt.close() 