import sys
from scipy.stats.stats import pearsonr
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.mlab as mlab

sample_map = {
    'Ost1-1':'A1054_01',
    'Ost1-1pr2':'A1054_07',
    'Ost2-2':'A1054_05',
    'Ost3L41':'A1054_06',
    'Ost16':'A1054_08',
    'Ost17':'A1054_09',
    'Ost19':'A1054_10',
    'Ost33':'A1054_16',
    'Ost26':'A1054_12',
    'Ost27':'A1054_13',
    'Ost34':'A1054_17',
    'Ost35':'A1054_18',
    'Ost37':'A1054_19',
    'Ost38':'A1054_04',
    'Ost39':'A1054_20',
    'Ost40':'A1054_21',
}

def main(f1, f2):
    fh = open(f1, 'r')
    zic1_expression = []
    other_expression = {}
    for line in fh:
        tokens = line.split('\t')
        if line.startswith('206373_at'):
            #zic1_expression = tokens[1:25]
            i = 0
            while i < len(tokens):
                if (i-1) % 4 == 0:
                    zic1_expression.append(tokens[i])
                i+=1
            #zic1_expression = map(float, zic1_expression)
            zic1_expression_copy = [zic1_expression[0], 
                    zic1_expression[1],
                    zic1_expression[7],
                    zic1_expression[21],
                    zic1_expression[3],
                    zic1_expression[4],
                    zic1_expression[5],
                    zic1_expression[14],
                    zic1_expression[10],
                    zic1_expression[11],
                    zic1_expression[15],
                    zic1_expression[16],
                    zic1_expression[18],
                    zic1_expression[19],
                    zic1_expression[20],
                    zic1_expression[23],
                    ]
            zic1_expression = map(float, zic1_expression_copy)
    fh.close()
    mir_expressions = {}
    fh = open(f2, 'r')
    for line in fh:
        if not line.startswith('Detector'):
            tokens = line.rstrip().split('\t')
            probe_id = tokens[0]
            expressions = [tokens[3],
                    tokens[9],
                    tokens[7],
                    tokens[8],
                    tokens[10],
                    tokens[11],
                    tokens[12],
                    tokens[18],
                    tokens[14],
                    tokens[15],
                    tokens[19],
                    tokens[20],
                    tokens[21],
                    tokens[6],
                    tokens[22],
                    tokens[23],
                ]
            mir_expressions[probe_id] = map(float, expressions)
    print len(mir_expressions)
    for key in mir_expressions.keys():
        corr = pearsonr(zic1_expression, mir_expressions[key])
        if corr[0] > 0.9 or corr[0] < -0.9:
            plotting(zic1_expression, mir_expressions[key], key, 'sig1')
        if corr[1] < 0.05:
            plotting(zic1_expression, mir_expressions[key], key, 'sig2')
        if corr[1] < (0.05/762):
            plotting(zic1_expression, mir_expressions[key], key, 'sig3')
    #plotting(zic1_expression, other_expression)
    #print zic1_expression
    #for key in other_expression.keys():
    #    print key
    #    corr = pearsonr(zic1_expression, other_expression[key])
    #    print corr
    #plotting(zic1_expression, other_expression['202310_s_at'], '206373_at', '202310_s_at', 0.84823378890913925, 1.6509007558736995e-07)

def plotting(zic1,comparator,name,sig_level):
    """docstring for plotting"""
    corr = pearsonr(zic1, comparator)
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
    ax.plot(range(1,17), zic1, 'go-', label='Zic1')
    ax.plot(range(1,17), comparator, 'r^-', label=name)
    # add the legend (bottom right for these plots)
    ax.legend(loc=4) 
    ax.set_xlim((0,17))
    # Save the generated Scatter Plot to a PNG file.
    if sig_level == 'sig1':
        folder = 'mir_r_cutoff'
    if sig_level == 'sig2':
        folder = 'mir_p_cutoff'
    if sig_level == 'sig3':
        folder = 'mir_corrp_cutoff'
    filename = name.split()[0]
    filename = filename.split('"')[1]
    print filename
    canvas.print_figure(folder+'/'+filename+'.png',dpi=500)
    

if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        main('processed_filtered.txt', 'mir_data.txt')
