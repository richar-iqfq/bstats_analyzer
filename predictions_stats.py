from utils.PredictionAnalyzer import PredictionAnalyzer

database = 'Net_5Hlayer_b_opt4_Full_R7_A030_58046_FPredictions.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27.npy'
alpha = -0.27
b = 'b4'
out = '_A0252_58046'
percent = 0.08

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
mb.make_plots(save=True, show=True)
