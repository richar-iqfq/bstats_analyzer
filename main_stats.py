from utils.StatisticAnalyser import StatisticAnalyser
import os

if __name__=='__main__':
    output_path = 'a-0.5_results'
    percent = (0.9, 0.95, 1, 1.05, 1.1)

    tolerance = {
        'Ne_calc' : 1E-4,
        'Energy_S' : 1E-7,
        'Energy_Savrg' : 1E-7,
        'Recover' : 1E-7
    }
    
    database = 'Database_AA.csv'
    energies = 'Database_Energies_AA.feather'

    sa = StatisticAnalyser(database, energies, output_path, alpha=-0.5, percent=percent)
    sa.plot_dispersion(save=True, show=False)
    sa.plot_correlation_matrix(save=True, show=False)
    #sa.plot_Err(cores='All', save=True, show=False, tolerance=tolerance)