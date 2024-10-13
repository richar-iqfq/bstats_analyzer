from utils.StatisticAnalyser import StatisticAnalyser
from utils.PlotMerge import PlotMerge
import os

if __name__=='__main__':
    mode = '' # '', 'odd' or 'even'
    step = ['energy'] # dispersion, energy 

    alpha_total_optimization = [
        -0.06, -0.08, -0.1, -0.12, -0.13, -0.133, -0.136, -0.139,
        -0.14, -0.142, -0.145, -0.148, -0.151, -0.154, -0.157, -0.16,
        -0.2, -0.22, -0.24, -0.243, -0.246, -0.249, -0.252, -0.255,
        -0.258, -0.261, -0.264, -0.267, -0.27, -0.3, -0.303, -0.306,
        -0.309, -0.312, -0.315, -0.318, -0.321, -0.324, -0.327, -0.33,
        -0.36, -0.39, -0.42, -0.45
    ]

    alpha_odd_optimization = [
        -0.3, -0.303, -0.306, -0.309, 
        -0.312, -0.315, -0.318, -0.321,
        -0.324, -0.327, -0.36, -0.39,
        -0.42, -0.47, -0.5, -0.52, -0.54,
        -0.56, -0.58, -0.6, -0.63, -0.66,
        -0.69, -0.72, -0.75, -0.78, -0.81,
        -0.58, -0.33, -0.43, -0.73, -0.83, -0.55, -0.61,
        -0.47-0.07, -0.22, -0.32, -0.62, -0.72, -0.44, -0.5,
        -0.36, -0.11, -0.21, -0.51, -0.61, -0.33, -0.39,
        -0.72, -0.47, -0.57, -0.87, -0.97, -0.69, -0.75
    ]

    # -0.321
    alpha_even_optimization = [
        0.6, 0.4, 0.3, 0.25, 0.22, 0.2, 
        0.18, 0.166, 0.15, 0.133, 0.1, 
        0.066, 0.033, 0, -0.01, 
        -0.02, -0.035, -0.04, 
        -0.045, -0.05, -0.07, -0.1,
        -0.13, -0.133, -0.136, -0.145,
        -0.148, -0.151, -0.154, -0.157,
        0.18, -0.07, 0.03, 0.33, 0.43, 0.15, 0.21,
        0.066, -0.184, -0.084, 0.216, 0.316, 0.096, 0.036
    ]

    # -0.3, -0.303, -0.306,
    # -0.309, -0.312, -0.315, -0.318,
    # -0.324, -0.327

    alpha_total_recovering = [
        -0.06, -0.08, -0.1, -0.12, -0.13, -0.133, -0.136,
        -0.139, -0.14, -0.142, -0.145, -0.148, -0.151,
        -0.154, -0.157, -0.16, -0.22, -0.24, -0.243, 
        -0.246, -0.249, -0.252, -0.255, -0.258, -0.261, 
        -0.264, -0.267, -0.27, -0.3, -0.303, -0.306, 
        -0.309, -0.312, -0.315, -0.318, -0.321, -0.324,
        -0.327,-0.33, -0.36, -0.39, -0.42, -0.45
    ]

    alpha_even_recovering = [
        -0.03, -0.035, -0.04,
        -0.045, -0.05,
        0.18, -0.07, 0.03, 0.33, 0.43, 0.15, 0.21,
        0.066, -0.184, -0.084, 0.216, 0.316, 0.096, 0.036
    ]

    alpha_odd_recovering = [
        -0.36, -0.47,
        -0.58, -0.33, -0.43, -0.73, -0.83, -0.55, -0.61,
        -0.47, -0.22, -0.32, -0.62, -0.72, -0.44, -0.5,
        -0.36, -0.11, -0.21, -0.51, -0.61, -0.33, -0.39,
        -0.72, -0.47, -0.57, -0.87, -0.97, -0.69, -0.75
    ]

    percent = [
        0.9, 0.95, 0.96, 0.97, 1,
        1.03, 1.04, 1.05, 1.1
    ]

    percent = [
        0.97, 1.1
    ]

    stats = StatisticAnalyser()
    merger = PlotMerge()

    alpha_optimization = {
        '' : list(dict.fromkeys(alpha_total_optimization)),
        'odd' : list(dict.fromkeys(alpha_odd_optimization)),
        'even' : list(dict.fromkeys(alpha_even_optimization))
    }

    alpha_recovering = {
        '' : list(dict.fromkeys(alpha_total_recovering)),
        'odd' : list(dict.fromkeys(alpha_odd_recovering)),
        'even' : list(dict.fromkeys(alpha_even_recovering))
    }

    #################################################
    ############# Dispersion ########################
    #################################################

    if 'dispersion' in step:
        # Dispersion by alpha
        alpha_optimization[mode].sort()
        for a in alpha_optimization[mode]:
            if mode == 'odd':
                database = f'results/optimization/a{a}_odd_results/results_a{a}.csv'
            elif mode == 'even':
                database = f'results/optimization/a{a}_even_results/results_a{a}.csv'
            else:
                database = f'results/optimization/a{a}_results/results_a{a}.csv'
            
            print(f'Ploting {a}')
            stats.plot_dispersion(database, a, show=False, save=True, append=True, mode=mode)

        print('+'*50, end='\n')
        # Merge dispersion
        print('Making dispersion merge...')
        merger.mergeDispersion(alpha_optimization[mode], mode=mode)
        
        print('+'*50, end='\n')
        # # General dispersion stats
        print('Ploting general dispersion stats')
        stats.plot_general_dispersion_stats(mode=mode, show=False)

        print('+'*50, end='\n')
        #  General dispersion merge
        print('Making general dispersion merge')
        input_path = os.path.join(
            'results', 'optimization',
            'Statistic_Plots', 'general',
            f'by_{mode}metric'
        )

        output_path = os.path.join(
            'results', 'optimization',
            'Statistic_Plots', 'general',
            f'Dispersion_{mode}metrics.pdf'
        )

        merger.mergeFiles(input_path, output_path)

        print('+'*50, end='\n')

    #######################################################
    ####################### Recovering ####################
    #######################################################

    if 'energy' in step:
        # Energy recover by alpha
        alpha_recovering[mode].sort()
        for a in alpha_recovering[mode]:
            if mode == 'odd':
                database = f'results/optimization/a{a}_odd_results/results_a{a}.csv'
            elif mode == 'even':
                database = f'results/optimization/a{a}_even_results/results_a{a}.csv'
            else:
                database = f'results/optimization/a{a}_results/results_a{a}.csv'

            print(f'Ploting recover for alpha = {a}')
            stats.plot_recovered_err(database, a, percent, show=False, save=True, mode=mode)

        print('+'*50, end='\n')
        # Energy recover merge
        print('Making recovery merge...')
        merger.mergeRecover(alpha_recovering[mode], mode=mode)

        print('+'*50, end='\n')
        # General energy recover stats
        print('Ploting general recovering stats')
        stats.plot_general_recovering_stats(alpha_recovering[mode], percent, show=False, save=True, mode=mode)

        print('+'*50, end='\n')
        # General energy recover merge
        print('Making general recovering merge')
        input = os.path.join(
            'results', 'recovering',
            'Statistic_Plots', 'general',
            f'{mode}partial'
        )

        output = os.path.join(
            'results', 'recovering',
            'Statistic_Plots', 'general',
            f'Alpha_{mode}energy_recovering.pdf'
        )

        merger.mergeFiles(input, output)
