from utils.StatisticAnalyser import StatisticAnalyser
from utils.PlotMerge import PlotMerge

if __name__=='__main__':
    alpha_optimization = [
        -0.06, -0.08, -0.1, -0.12, -0.13, -0.133, -0.136, -0.139,
        -0.142, -0.145, -0.148, -0.151, -0.154, -0.157, -0.14, -0.16,
        -0.2, -0.22, -0.24, -0.27, -0.3, -0.303, -0.306, -0.309, -0.312, 
        -0.315, -0.318, -0.321, -0.324, -0.327, -0.33, -0.36, -0.39, -0.42,
        -0.45
    ]

    alpha_recovering = [-0.06, -0.08, -0.1, -0.12, -0.14, -0.16, -0.22, 
            -0.24, -0.27, -0.33, -0.36, -0.39, -0.42, -0.45
    ]

    percent = [
        0.9, 0.95, 0.96, 0.97, 1,
        1.03, 1.04, 1.05, 1.1
    ]

    stats = StatisticAnalyser()
    merger = PlotMerge()

    for a in alpha_optimization:
        database = f'results/optimization/a{a}_results/results_a{a}.csv'
        
        print(f'Ploting {a}')
        stats.plot_dispersion(database, a, show=False, save=True)

    print('Making dispersion merge...')
    merger.mergeDispersion(alpha_optimization)

    # for a in alpha_recovering:
    #     database = f'results/optimization/a{a}_results/results_a{a}.csv'

    #     print(f'\nPloting recover for alpha = {a}')
    #     stats.plot_recovered_err(database, a, percent, show=False, save=True)

    # print('Making recovery merge...')
    # merger.mergeRecover(alpha_recovering)