from utils.StatisticAnalyser import StatisticAnalyser
from utils.PlotMerge import PlotMerge

if __name__=='__main__':
    alpha = [-0.06, -0.08, -0.1, -0.12, -0.14, -0.16, -0.2, -0.22, 
            -0.24, -0.27, -0.3, -0.33, -0.36, -0.39, -0.42, -0.45
    ]

    stats = StatisticAnalyser()
    merger = PlotMerge()

    for a in alpha:
        database = f'results/optimization/a{a}_results/results_a{a}.csv'
        
        print(f'Ploting {a}')
        stats.plot_dispersion(database, a, show=False, save=True)

    print('Making merge...')
    merger.mergeDispersion(alpha)