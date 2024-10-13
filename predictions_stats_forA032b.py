from utils.PredictionAnalyzer import PredictionAnalyzer

database = 'data/Net_5Hlayerb_opt4g_A032b_FPredictions.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27_A032b.npy'
alpha = -0.27
b = 4
out = '_A032b_FPredictions'
percent = 0.035

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
# mb.make_plots(save=True, show=True)
# mb.make_plotsv2(save=False, show=True)
mb.make_plotsv3('A032b',
                save=False, 
                show=True,
                strategy='statistical_difference',
                strategy_perc=0.07)

# mb.plot_b_sigma(save=False, strategy='statistical_difference')