from utils.StatisticAnalyser import StatisticAnalyser
import os

if __name__=='__main__':
    cores = 14
    calc_types=(1, 2, 3, 4, 5)
    alpha = -0.3
    
    percent = (0.9, 0.95, 0.96, 0.97, 1, 1.03, 1.04, 1.05, 1.1)

    tolerance = {
        'Ne_calc' : 1E-4,
        'Energy_S' : 1E-7,
        'Energy_Savrg' : 1E-7,
        'Recover' : 1E-7
    }

    database = 'Database_AA.csv'
    energies = 'Database_Energies_AA.feather'

    output_path = f'a{alpha}_results'
    out_file = os.path.join('results', 'optimization', output_path, f'results_a{alpha}.csv')

    st = StatisticAnalyser(out_file, energies, output_path, alpha=alpha, percent=percent)
    st.plot_Err(cores=cores, save=True, show=False, tolerance=tolerance, isdata=True)