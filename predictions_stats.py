from utils.PredictionAnalyzer import PredictionAnalyzer

code = 'b033cn'
alpha = -0.33
b = 1
hidden = 4

# database = 'data/FPredictions_b033cn_egqLXKXflvRSfnjwLhuqUrdFH.csv'
database = 'data/FPredictions_b022cn_BCDytstiomkCTSTUgSUUpobbp.csv'

# predictions_file = f'predictions/saved_values/Ecorrb1_b033cn_a-0.33_egqLXKXflvRSfnjwLhuqUrdFH.npy'
predictions_file = 'predictions/saved_values/Ecorrb1_b022cn_a-0.22_BCDytstiomkCTSTUgSUUpobbp.npy'
out = f'_{code}_FPredictions'
percent = 0.035

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
# mb.make_plots(save=False, show=True)
# mb.make_plotsv2(save=False, show=True)
mb.make_plotsv3(code,
                save=False, 
                show=True,
                strategy='statistical_difference',
                strategy_perc=0.07)

# mb.plot_b_sigma(save=False, strategy='statistical_difference')e
