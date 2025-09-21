import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_squared_error
import os
from matplotlib.ticker import PercentFormatter

class PredictionAnalyzer():
    def __init__(self, database, predictions_file, alpha, b, out=None, percent=0.08):
        self.alpha = alpha
        self.database = pd.read_csv(database)
        self.energies = [val for val in self.database['CIe']]
        self.predictions = np.load(predictions_file)
        self.energies_pred = self.predictions[:,0]
        self.b = b
        self.out = out
        self.percent = percent

    def plot_b_sigma(self, save=False, strategy='percentual_difference', perc=0.08):
        fig, ax = plt.subplots(2)
        fig.suptitle(f'Using percentage of {perc*100}%' if strategy=='percentual_difference' else f'Using abs(b - b_pred) > mean + 3*sigma')
        fig.set_size_inches(20, 13)
        fig.subplots_adjust(top=0.95,
                            bottom=0.05,
                            left=0.05,
                            right=0.95)

        b, b_pred = self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred']

        # boolean_outliers = np.zeros(len(b), dtype=bool)

        # Get difference
        difference = abs(b - b_pred)

        # Get stat metrics
        mean = np.mean(difference)
        sigma = np.std(difference)
        diff_max = np.max(difference)

        # Set cut value
        cut = diff_max*perc if perc else mean + 3*sigma
        
        # # Get outliers
        # for i, diff_val in enumerate(difference):
        #     if diff_val > cut:
        #         boolean_outliers[i] = True

        boolean_outliers = self._get_outliers(b, b_pred, strategy=strategy, perc=perc)

        inverse_boolean = np.invert(boolean_outliers)

        # Get limit lines
        b_upper = [val+cut for val in b]
        b_lower = [val-cut for val in b]

        ax[0].plot(b, b, '-b')
        ax[0].plot(b, b_upper, '-r')
        ax[0].plot(b, b_lower, '-y')

        ax[0].plot(b, b_pred, '.g')
        ax[0].plot(b[boolean_outliers], b_pred[boolean_outliers], '.r')

        ax[0].set_title(f'Total: {len(b)}     Outliers: {np.count_nonzero(boolean_outliers)}')

        ax[1].plot(b, b, '-b')
        ax[1].plot(b[inverse_boolean], b_pred[inverse_boolean], '.y')

        ax[1].set_title(f'Final Molecules: {len(b[inverse_boolean])}')

        for axis in ax:
            axis.set_xlabel('b')
            axis.set_ylabel('b_pred')

        plt.show()

        if save:
            path = os.path.join('predictions', 'plots', f"a{self.alpha}_{self.b}{self.out}_outliers_{'perc_'+str(perc) if perc else '3sigma'}.pdf")
            fig.savefig(path, format='pdf')
        
    def _get_outliers(self, y, y_pred, strategy='percentual_error', perc=0.08):
        '''
        Get outlier values

        Parameters
        ----------
        y `array`:
            Real values
        
        y_pred `array`:
            Predicted values

        strategy `str`:
            Algorithm to remove outliers, can be `percentual_error`, `percentual_difference` and
            `statistical_difference` 
        '''
        y = np.array(y)
        y_pred = np.array(y_pred)

        tol = self.percent*100

        # Get difference
        difference = abs(y - y_pred)

        # Get stat metrics
        mean = np.mean(difference)
        sigma = np.std(difference)
        diff_max = np.max(difference)
        
        if strategy == 'percentual_error':
            # Get percentual error
            if y.all():
                err = abs( (y - y_pred)/(y) )*100

                # Get boolean array
                boolean_outliers = err > tol
    
            else:
                boolean_outliers = np.zeros_like(y, dtype=bool)

                for i, y_value in enumerate(y):
                    if y_value != 0:
                        err = abs( (y_value - y_pred[i])/(y_value) )*100
                    else:
                        err = np.inf

                    if err > tol:
                        boolean_outliers[i] = True

        elif strategy == 'percentual_difference':
            # Set cut value by a percentage
            cut = diff_max*perc

            # Get outliers
            boolean_outliers = difference > cut
            
        elif strategy == 'statistical_difference':
            # Set cut value by statistic significance
            cut = mean + 3*sigma

            # Get outliers
            boolean_outliers = difference > cut

        else:
            raise Exception('Invalid strategy requested...')

        return boolean_outliers
    
    def __retire_outliers(self, y, y_pred, boolean_outliers):
        boolean_outliers = np.invert(boolean_outliers)

        y = y[boolean_outliers]
        y_pred = y_pred[boolean_outliers]

        boolean_outliers = np.ones(y.shape)

        return y, y_pred, boolean_outliers

    def __save_img_outliers(self, smiles_ID):
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit import RDLogger

        # Disable logs and errors
        RDLogger.DisableLog('rdApp.*')

        smiles_df = pd.read_csv('data/smiles.csv')

        img_path = os.path.join('predictions', 'plots', 'img_outliers')
        if not os.path.isdir(img_path):
            os.makedirs(img_path)

        for i, ID in enumerate(smiles_ID.values):
            row = smiles_df[smiles_df['ID'] == ID]

            try:
                ID_smile = row['smiles'].values[0]
            except:
                print(f'Missing ID: {ID}')
                continue

            Mol = Chem.MolFromSmiles(ID_smile)

            file_name = os.path.join(img_path, f'{i+1}_{ID}.jpg')

            if Mol:
                img = Draw.MolToImage(Mol)
                img.save(file_name)
            else:
                print(f'Error with: {ID}')

    def __count_species(self, IDs):
        anions = 0
        cations = 0
        radicals = 0
        neutros = 0

        for ID in IDs:
            if 'a' in ID:
                anions += 1
            elif 'c' in ID:
                cations += 1
            elif 'r' in ID:
                radicals += 1
            else:
                neutros += 1

        print('\n----Count----')
        print(f'anions: {anions}')
        print(f'cations: {cations}')
        print(f'radicals: {radicals}')
        print(f'neutros: {neutros}')

    def make_plots(self, show=False, save=True, retire_outliers=True):
        '''
        Build and save the resulting plots

        Parameters
        ----------
        show `bool`
            If true display the plots in screen, default is False

        save `bool`
            If true save plots in pdf format inside prediction/plots
        '''
        path = os.path.join('predictions', 'plots')
        if not os.path.isdir(path):
            os.makedirs(path)
        
        #------------------------------- Err ------------------------------
        enrg = np.array(self.energies)
        enrg_pred = np.array(self.energies_pred)

        initial_err = abs(((enrg_pred - enrg)/enrg)*100)

        boolean_outliers = self._get_outliers(enrg, enrg_pred)
        outliers_count = np.count_nonzero(boolean_outliers)

        outliers_ID = self.database['ID'][boolean_outliers]
        
        self.__count_species(self.database['ID'][np.invert(boolean_outliers)])

        if retire_outliers:
            energies, energies_pred, _ = self.__retire_outliers(enrg, enrg_pred, boolean_outliers)

        else:
            energies, energies_pred = enrg, enrg_pred

        MAE_e = mean_absolute_error(energies, energies_pred)
        MAE_b = mean_absolute_error(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred'])

        PE_e = np.sqrt(mean_squared_error(energies, energies_pred))
        PE_b = np.sqrt(mean_squared_error(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred']))

        Err = abs(((energies_pred - energies)/energies)*100)
        maxErr = max(Err)
        meanErr = np.mean(Err)
        
        fig, ax = plt.subplot_mosaic('AABB;CCCC')
        fig.suptitle(f'{self.b} alpha: {self.alpha}')
        fig.set_size_inches(20, 13)
        fig.subplots_adjust(wspace=0.320)
        
        # Correlation b plot
        ax['A'].scatter(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred'], s=2, color='y')
        ax['A'].plot(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}'], '-b', linewidth='0.7')
        
        r2 = np.corrcoef(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred'])

        ax['A'].set_xlabel('b_opt')
        ax['A'].set_ylabel('b Predicted')
        ax['A'].set_title(f'RMSE: {PE_b:.2f}      r$^2$: {r2[0,1]:.3f}')

        # Correlation energie plot
        ax['B'].scatter(energies, energies_pred, s=2, color='g')
        ax['B'].plot(energies, energies, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(energies, energies_pred)

        print('\n----Metrics----')
        print(f'r2: {r2[0,1]}')
        print(f'RMSE: {PE_e}')
        print((1-PE_e)*100)
        print(f'MAE: {MAE_e}\n')

        ax['B'].set_xlabel('Energies')
        ax['B'].set_ylabel('Energies Predicted')
        ax['B'].set_title(f'r$^2$: {r2[0,1]:.3f}\nRMSE: {PE_e:.2f}    Acc: {(1-PE_e)*100:.2f}%\nFinal_Molecules: {len(energies)}    Outliers_retired: {outliers_count}')

        # # B Err plot
        # b_Err = ((self.database['b_pred'] - self.database[f'B_opt{self.b}'])/self.database[f'B_opt{self.b}'])*100
        # maxBErr = max(b_Err)
        # ax['C'].hist(b_Err, 15, color='y')
        # ax['C'].set_xlabel('Percentage Error')
        # ax['C'].set_ylabel('Count')
        # ax['C'].set_title(f'Error over b_opt\nMAE = {MAE_b:.4f}      Max Error: {maxBErr:.2f}')

        # Energie Err plot
        values, bins, bars = ax['C'].hist(Err, 40, color='g', weights=np.ones(len(Err)) / len(Err))
        ax['C'].yaxis.set_major_formatter(PercentFormatter(1))
        ax['C'].set_xlabel('Absolute Percentage Error')
        
        sum = 0
        Accumulate = np.zeros(40)
        y_height = np.zeros(40)
        for i, val in enumerate(values):
            sum += val
            Accumulate[i] = sum
            y_height[i] = val

        rects = ax['C'].patches
        x_center = np.zeros(len(rects))
        for i, rect in enumerate(rects):
            x_center[i] = rect.get_x() + rect.get_width()/2

        ax['C'].plot(x_center, y_height, '-b')

        for i in range(len(x_center)):
            ax['C'].text(x_center[i], y_height[i]+0.001, f'{Accumulate[i]*100:.1f}%', fontdict={'size': 6,})

        # ax['C'].bar_label(bars, fontsize=20, color='navy')
        ax['C'].axvline(Err.mean(), color='k', linestyle='dashed', linewidth=1)
        ax['C'].set_title(f'Error over Energie    MAE = {MAE_e:.4f}\nAbsolute Mean Percentage Error: {meanErr:.2f}%    Absolute Max Error: {maxErr:.2f}%')

        if save:
            self.__save_img_outliers(outliers_ID)
            fig.savefig(os.path.join(path, f'a{self.alpha}_{self.b}{self.out}_{self.percent}.pdf'), dpi=450, format='pdf')

        # Err as a function of b
        fig_Err_b, ax_Err_b = plt.subplots(1)
        fig_Err_b.set_size_inches(20, 13)
        ax_Err_b.set_title('Error as a function of b')

        ax_Err_b.plot(self.database[f'B_opt{self.b}_pred'], initial_err, '.g')
        ax_Err_b.set_ylim(0, 1.0)
        ax_Err_b.set_xlabel('b')
        ax_Err_b.set_ylabel('Percentage error')

        if save:
            fig_Err_b.savefig(os.path.join(path, f'Error_vs_b_{self.alpha}_{self.b}{self.out}_{self.percent}.pdf'))

        if not os.path.isfile(os.path.join(path, f'Distribution_{self.alpha}_{self.b}.pdf')):
            # b and Ecorr relation with Ne
            fig_b, ax = plt.subplots(2)
            fig_b.set_size_inches(20, 13)
            fig_b.suptitle('Distribution')

            ax[0].plot(self.database['Ne'], self.database[f'B_opt{self.b}'], '.r')
            ax[0].set_xlabel('Electrons')
            ax[0].set_ylabel('b opt')

            ax[1].plot(self.database['Ne'], self.database['CIe'], '.g')
            ax[1].set_xlabel('Electrons')
            ax[1].set_ylabel('Ecorr')

            if save:
                fig_b.savefig(os.path.join(path, f'Distribution_{self.alpha}_{self.b}.pdf'), dpi=450, format='pdf')

        if not os.path.isfile(os.path.join(path, f'Ne_Distribution_{self.alpha}_{self.b}.pdf')):
            # Ne dispersion
            fig_c, ax = plt.subplots(1)
            fig_c.set_size_inches(20, 13)
            values, bins, bars = ax.hist(self.database['Ne'], max(self.database['Ne']), color='b', weights=np.ones(len(self.database['Ne'])) / len(self.database['Ne']))
            ax.yaxis.set_major_formatter(PercentFormatter(1))
            ax.set_xlabel('Ne')
            
            # sum = 0
            # Accumulate = np.zeros(max(self.database['Ne']))
            # y_height = np.zeros(max(self.database['Ne']))
            # for i, val in enumerate(values):
            #     sum += val
            #     Accumulate[i] = sum
            #     y_height[i] = val

            # rects = ax.patches
            # x_center = np.zeros(len(rects))
            # for i, rect in enumerate(rects):
            #     x_center[i] = rect.get_x() + rect.get_width()/2

            # ac_i = 1
            # for i in range(len(x_center)):
            #     if ac_i != Accumulate[i]:
            #         ax.text(x_center[i], y_height[i]+0.001, f'{Accumulate[i]*100:.2f}%', fontdict={'size': 6,})
            #         ac_i = Accumulate[i]
                
            # ax['C'].bar_label(bars, fontsize=20, color='navy')
            ax.axvline(self.database['Ne'].mean(), color='k', linestyle='dashed', linewidth=1)
            ax.set_title(f"Number of electrons\nMean = {self.database['Ne'].mean():.0f}")

            if save:
                fig_c.savefig(os.path.join(path, f'Ne_Distribution_{self.alpha}_{self.b}.pdf'), dpi=450, format='pdf')

        if show:
            plt.show()

    def make_plotsv2(self, show=False, save=True, retire_outliers=True):
        '''
        Build and save the resulting plots

        Parameters
        ----------
        show `bool`
            If true display the plots in screen, default is False

        save `bool`
            If true save plots in pdf format inside prediction/plots
        '''
        path = os.path.join('predictions', 'plots')
        if not os.path.isdir(path):
            os.makedirs(path)
        
        #------------------------------- Err ------------------------------
        enrg = np.array(self.energies)
        enrg_pred = np.array(self.energies_pred)

        initial_err = abs(((enrg_pred - enrg)/enrg)*100)

        boolean_outliers = self._get_outliers(enrg, enrg_pred)
        outliers_count = np.count_nonzero(boolean_outliers)

        outliers_ID = self.database['ID'][boolean_outliers]
        
        self.__count_species(self.database['ID'][np.invert(boolean_outliers)])

        if retire_outliers:
            energies, energies_pred, _ = self.__retire_outliers(enrg, enrg_pred, boolean_outliers)
            b_values, b_values_pred, _ = self.__retire_outliers(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}pred'], boolean_outliers)

        else:
            energies, energies_pred = enrg, enrg_pred
            b_values, b_values_pred = self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}pred']

        MAE_e = mean_absolute_error(energies, energies_pred)
        MAE_b = mean_absolute_error(b_values, b_values_pred)

        PE_e = np.sqrt(mean_squared_error(energies, energies_pred))
        PE_b = np.sqrt(mean_squared_error(b_values, b_values_pred))

        Err = abs(((energies_pred - energies)/energies)*100)
        maxErr = max(Err)
        meanErr = np.mean(Err)
        
        fig, ax = plt.subplot_mosaic('AABB;CCCC')
        fig.suptitle(f'b {self.b} alpha: {self.alpha}')
        fig.set_size_inches(20, 13)
        fig.subplots_adjust(wspace=0.320)
        
        # Correlation b plot
        ax['A'].scatter(b_values, b_values_pred, s=2, color='y')
        ax['A'].plot(b_values, b_values, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(b_values, b_values_pred)

        ax['A'].set_xlabel('b_opt')
        ax['A'].set_ylabel('b Predicted')
        ax['A'].set_title(f'RMSE: {PE_b:.2f}      r$^2$: {r2[0,1]:.3f}\nData: {len(b_values)}')

        # Correlation energie plot
        ax['B'].scatter(energies, energies_pred, s=2, color='g')
        ax['B'].plot(energies, energies, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(energies, energies_pred)

        print('\n----Metrics----')
        print(f"b_max: {max(self.database[f'B_opt{self.b}'])}")
        print(f"b_min: {min(self.database[f'B_opt{self.b}'])}")
        print(f'r2: {r2[0,1]}')
        print(f'RMSE: {PE_e}')
        print((1-PE_e)*100)
        print(f'MAE: {MAE_e}\n')

        ax['B'].set_xlabel('Energies')
        ax['B'].set_ylabel('Energies Predicted')
        ax['B'].set_title(f'r$^2$: {r2[0,1]:.3f}\nRMSE: {PE_e:.2f}    Acc: {(1-PE_e)*100:.2f}%\nFinal_Molecules: {len(energies)}    Outliers_retired: {outliers_count}')

        # Energie Err plot
        values, bins, bars = ax['C'].hist(Err, 40, color='g', weights=np.ones(len(Err)) / len(Err))
        ax['C'].yaxis.set_major_formatter(PercentFormatter(1))
        ax['C'].set_xlabel('Absolute Percentage Error')
        
        sum = 0
        Accumulate = np.zeros(40)
        y_height = np.zeros(40)
        for i, val in enumerate(values):
            sum += val
            Accumulate[i] = sum
            y_height[i] = val

        rects = ax['C'].patches
        x_center = np.zeros(len(rects))
        for i, rect in enumerate(rects):
            x_center[i] = rect.get_x() + rect.get_width()/2

        ax['C'].plot(x_center, y_height, '-b')

        for i in range(len(x_center)):
            ax['C'].text(x_center[i], y_height[i]+0.001, f'{Accumulate[i]*100:.1f}%', fontdict={'size': 6,})

        # ax['C'].bar_label(bars, fontsize=20, color='navy')
        ax['C'].axvline(Err.mean(), color='k', linestyle='dashed', linewidth=1)
        ax['C'].set_title(f'Error over Energie    MAE = {MAE_e:.4f}\nAbsolute Mean Percentage Error: {meanErr:.2f}%    Absolute Max Error: {maxErr:.2f}%')

        if save:
            fig.savefig(os.path.join(path, f'a{self.alpha}_{self.b}{self.out}_{self.percent}_v2.pdf'), dpi=450, format='pdf')

        if show:
            plt.show()

    def make_plotsv3(self, ID, show=False, save=True, retire_outliers=True, strategy='percentual_difference', strategy_perc=0.08):
        '''
        Build and save the resulting plots

        Parameters
        ----------
        show `bool`
            If true display the plots in screen, default is False

        save `bool`
            If true save plots in pdf format inside prediction/plots
        '''
        path = os.path.join('predictions', 'plots')
        if not os.path.isdir(path):
            os.makedirs(path)
        
        #Retire previous outliers
        boolean_outliers = self._get_outliers(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred'], strategy=strategy, perc=strategy_perc)

        initial_retired = np.count_nonzero(boolean_outliers)

        b, b_pred, _ = self.__retire_outliers(self.database[f'B_opt{self.b}'], self.database[f'B_opt{self.b}_pred'], boolean_outliers)
        prev_energies, prev_energies_pred, _ = self.__retire_outliers(np.array(self.energies), np.array(self.energies_pred), boolean_outliers)

        self.__count_species(self.database['ID'][np.invert(boolean_outliers)])

        #------------------------------- Err ------------------------------
        enrg = np.array(prev_energies)
        enrg_pred = np.array(prev_energies_pred)

        initial_err = abs(((enrg_pred - enrg)/enrg)*100)

        boolean_outliers = self._get_outliers(enrg, enrg_pred, strategy='percentual_error')
        outliers_count = np.count_nonzero(boolean_outliers)

        if retire_outliers:
            energies, energies_pred, _ = self.__retire_outliers(enrg, enrg_pred, boolean_outliers)
            b_values, b_values_pred, _ = self.__retire_outliers(b, b_pred, boolean_outliers)

        else:
            energies, energies_pred = enrg, enrg_pred
            b_values, b_values_pred = b, b_pred

        MAE_e = mean_absolute_error(energies, energies_pred)
        MAE_b = mean_absolute_error(b_values, b_values_pred)

        PE_e = np.sqrt(mean_squared_error(energies, energies_pred))
        PE_b = np.sqrt(mean_squared_error(b_values, b_values_pred))

        Err = abs(((energies_pred - energies)/energies)*100)
        maxErr = max(Err)
        meanErr = np.mean(Err)
        
        fig, ax = plt.subplot_mosaic('AABB;CCCC')
        fig.suptitle(f'ID: {ID}\nb: {self.b} alpha: {self.alpha}\nStrategy: {strategy}')
        fig.set_size_inches(20, 13)
        fig.subplots_adjust(wspace=0.320)
        
        # Correlation b plot
        ax['A'].scatter(b_values, b_values_pred, s=2, color='y')
        ax['A'].plot(b_values, b_values, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(b_values, b_values_pred)

        ax['A'].set_xlabel('$b_{opt}$')
        ax['A'].set_ylabel('$b_{pred}$')
        ax['A'].set_title(f'$RMSE$: {PE_b:.2f}      r$^2$: {r2[0,1]:.3f}\n'+'$b_{outliers}$: ' + f'{initial_retired}')

        # Correlation energie plot
        ax['B'].scatter(energies, energies_pred, s=2, color='g')
        ax['B'].plot(energies, energies, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(energies, energies_pred)

        print('\n----Metrics----')
        print(f"b_max: {max(b)}")
        print(f"b_min: {min(b)}")
        print(f'r2: {r2[0,1]}')
        print(f'RMSE: {PE_e}')
        print((1-PE_e)*100)
        print(f'MAE: {MAE_e}\n')

        ax['B'].set_xlabel('Energías (hartree)')
        ax['B'].set_ylabel('Energías predichas (hartree)')
        ax['B'].set_title(f'r$^2$: {r2[0,1]:.3f}\n$RMSE$: {PE_e:.2f}    Precisión: {(1-PE_e)*100:.2f}%\nMoléculas: {len(energies)}    ' + '$E_{outliers}$: ' + f'{outliers_count}')

        # Energie Err plot
        values, bins, bars = ax['C'].hist(Err, 40, color='g', weights=np.ones(len(Err)) / len(Err))
        ax['C'].yaxis.set_major_formatter(PercentFormatter(1))
        ax['C'].set_xlabel('Porcentaje de error absoluto (APE)')
        
        sum = 0
        Accumulate = np.zeros(40)
        y_height = np.zeros(40)
        for i, val in enumerate(values):
            sum += val
            Accumulate[i] = sum
            y_height[i] = val

        rects = ax['C'].patches
        x_center = np.zeros(len(rects))
        for i, rect in enumerate(rects):
            x_center[i] = rect.get_x() + rect.get_width()/2

        ax['C'].plot(x_center, y_height, '-b')

        for i in range(len(x_center)):
            ax['C'].text(x_center[i], y_height[i]+0.001, f'{Accumulate[i]*100:.1f}%', fontdict={'size': 6,})

        # ax['C'].bar_label(bars, fontsize=20, color='navy')
        ax['C'].axvline(Err.mean(), color='k', linestyle='dashed', linewidth=1)
        ax['C'].set_title(f'Error sobre energías\n$MAE$ = {MAE_e:.4f}    $MAPE$: {meanErr:.2f}%    Error absoluto máximo: {maxErr:.2f}%')

        if save:
            fig.savefig(os.path.join(path, f'a{self.alpha}_{self.b}{self.out}_{self.percent}_{strategy}.pdf'), dpi=300, format='pdf')

        if show:
            plt.show()