from utils.StatisticAnalyser import StatisticAnalyser

if __name__=='__main__':
    output_path = 'results/general'
    
    database = 'data/data.csv'
    energies = 'data/energies.feather'

    sa = StatisticAnalyser()
    sa.plot_correlation_matrix(database, output_path, save=True, show=True)