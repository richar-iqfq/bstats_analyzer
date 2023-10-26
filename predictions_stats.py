from utils.PredictionAnalyzer import PredictionAnalyzer

database = 'data/Net_5Hlayerb_opt4_A032_outliers_count_ae_5000epochs_volume.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27_outliers_A032.npy'
alpha = -0.27
b = 'b4'
out = '_A032_5_rs112254_outlierscount_ae_volume'
percent = 0.08

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
mb.make_plots(save=True, show=True)
