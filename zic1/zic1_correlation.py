import sys
from scipy.stats.stats import pearsonr
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.mlab as mlab



def main(f):
    fh = open(f, 'r')
    zic1_expression = []
    other_expression = {}
    for line in fh:
        tokens = line.split('\t')
        if line.startswith('206373_at'):
            i = 0
            while i < len(tokens):
                if (i-1) % 4 == 0:
                    zic1_expression.append(tokens[i])
                i+=1
            zic1_expression = map(float, zic1_expression)
        elif line.startswith('Scan'):
            pass
        else:
            other_expression[tokens[0]] = []
            i = 0
            while i < len(tokens):
                if (i-1) % 4 == 0:
                    other_expression[tokens[0]].append(tokens[i])
                i+=1
            other_expression[tokens[0]] = map(float, other_expression[tokens[0]])
    plotting(zic1_expression, other_expression)

def plotting(zic1,comparators):
    """docstring for plotting"""
    from mapping import probe_map
    for key in comparators.keys():
        corr = pearsonr(zic1, comparators[key])
        #the string of correlation stats
        s = 'R = '+str(corr[0])+'\nP = '+str(corr[1])
        # Create a figure with size 6 x 6 inches.
        fig = Figure(figsize=(6,6))
        # Create a canvas and add the figure to it.
        canvas = FigureCanvas(fig)
        # Create a subplot.
        ax = fig.add_subplot(111)
        # Set the title.
        ax.set_title(s,fontsize=10)
        # Set the X Axis label.
        ax.set_xlabel('Samples',fontsize=8)
        # Set the Y Axis label.
        ax.set_ylabel('Normalized Expression',fontsize=8)
        # Display Grid.
        ax.grid(True,linestyle='-',color='0.75')
        # Generate the Scatter Plot.
        ax.plot(range(1,25), zic1, 'go-', label=probe_map['206373_at'])
        ax.plot(range(1,25), comparators[key], 'r^-', label=probe_map[key])
        # add the legend
        ax.legend()
        #ax.text(0.1,max(zic1),s)
        # Save the generated Scatter Plot to a PNG file.
        canvas.print_figure('correlations/'+key+'.png',dpi=500)
    

if __name__ == '__main__':
    if sys.argv[1]:
        main(sys.argv[1])
    else:
        main('processed_filtered.txt')
