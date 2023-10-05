from utils.StatisticAnalyser import StatisticAnalyser

if __name__=='__main__':
    output_path = 'a-0.5_results'
    
    database = 'data/data.csv'
    energies = 'data/energies.feather'

    sa = StatisticAnalyser(database, energies)
    sa.plot_correlation_matrix('general', save=True, show=True)