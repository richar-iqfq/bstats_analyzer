from utils.StatisticAnalyser import StatisticAnalyser

if __name__=='__main__':
    output_path = 'a-0.27_results'
    
    database = 'data/results_a-0.27.csv'
    energies = 'data/energies.feather'

    sa = StatisticAnalyser()
    sa.plot_correlation_matrix(database, output_path, save=True, show=True)