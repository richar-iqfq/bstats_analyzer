from utils.StatisticAnalyser import StatisticAnalyser

alpha = -0.06
database = 'results/optimization/a-0.06_results/results_a-0.06.csv'

stats = StatisticAnalyser()

stats.plot_dispersion(database, alpha, show=True)