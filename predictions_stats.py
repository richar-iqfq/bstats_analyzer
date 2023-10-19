from utils.PredictionAnalyzer import PredictionAnalyzer

database = 'data/Net_6Hlayerb_opt4_A039_FPredictions.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27.npy'
alpha = -0.27
b = 'b4'
out = '_A039_random_state_r2val'
percent = 0.08

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
mb.make_plots(save=True, show=True)
