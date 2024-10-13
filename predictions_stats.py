from utils.PredictionAnalyzer import PredictionAnalyzer

code = 'b002cn'
alpha = -0.22
b = 1
hidden = 6

# database = f'data/FPredictions_b013sn_uAIRZytgceDqUrQNrXtQlfdUu.csv'
#database = f'data/FPredictions_b012sn_GqNYQtKzHXbvxdGbFGiqlELUw.csv'
database = 'data/FPredictions_b002cn_imtKlNImlmLtDMpmnTqfVckUr_tuningBatch_6H.csv'
predictions_file = f'predictions/saved_values/Ecorrb{b}_a{alpha}_{code}.npy'
out = f'_{code}_FPredictions'
percent = 0.035

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out, percent)
# mb.make_plots(save=True, show=True)
# mb.make_plotsv2(save=False, show=True)
mb.make_plotsv3(code,
                save=True, 
                show=True,
                strategy='statistical_difference',#'statistical_difference',
                strategy_perc=0.07)

# mb.plot_b_sigma(save=False, strategy='statistical_difference')e
