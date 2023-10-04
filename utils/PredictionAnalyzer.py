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

    def _get_outliers(self, y, y_pred):
        boolean_outliers = np.zeros(len(y), dtype=bool)

        for i, y_value in enumerate(y):
            tol = self.percent*100
            
            if y_value != 0:
                err = abs( (y_value - y_pred[i])/(y_value) )*100
            else:
                err = np.inf

            if err > tol:
                boolean_outliers[i] = True

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

        smiles_df = pd.read_csv('Database_SMILES.csv')

        img_path = os.path.join('predictions', 'plots', 'img_outliers')

        for i, ID in enumerate(smiles_ID.values):
            row = smiles_df[smiles_df['ID'] == ID]
            ID_smile = row['smiles'].values[0]

            Mol = Chem.MolFromSmiles(ID_smile)

            file_name = os.path.join(img_path, f'{i+1}_{ID}.jpg')

            if Mol:
                img = Draw.MolToImage(Mol)
                img.save(file_name)
            else:
                print(f'Error with: {ID}')

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
        #------------------------------- Err ------------------------------
        enrg = np.array(self.energies)
        enrg_pred = np.array(self.energies_pred)

        boolean_outliers = self._get_outliers(enrg, enrg_pred)
        outliers_count = np.count_nonzero(boolean_outliers)

        outliers_ID = self.database['ID'][boolean_outliers]
        
        if retire_outliers:
            energies, energies_pred, _ = self.__retire_outliers(enrg, enrg_pred, boolean_outliers)

        else:
            energies, energies_pred = enrg, enrg_pred
        
        MAE_e = mean_absolute_error(energies, energies_pred)
        MAE_b = mean_absolute_error(self.database['B_opt4'], self.database['b_pred'])

        PE_e = np.sqrt(mean_squared_error(energies, energies_pred))
        PE_b = np.sqrt(mean_squared_error(self.database['B_opt4'], self.database['b_pred']))

        Err = abs(((energies_pred - energies)/energies)*100)
        maxErr = max(Err)
        meanErr = np.mean(Err)
        
        fig, ax = plt.subplot_mosaic('AABB;CCCC')
        fig.suptitle(f'{self.b} alpha: {self.alpha}')
        fig.set_size_inches(20, 13)
        fig.subplots_adjust(wspace=0.320)
        
        # Correlation b plot
        ax['A'].scatter(self.database['B_opt4'], self.database['b_pred'], s=2, color='y')
        ax['A'].plot(self.database['B_opt4'], self.database['B_opt4'], '-b', linewidth='0.7')
        
        r2 = np.corrcoef(self.database['B_opt4'], self.database['b_pred'])

        ax['A'].set_xlabel('b_opt')
        ax['A'].set_ylabel('b Predicted')
        ax['A'].set_title(f'RMSE: {PE_b*100:.2f}      r$^2$: {r2[0,1]:.3f}')

        # Correlation energie plot
        ax['B'].scatter(energies, energies_pred, s=2, color='g')
        ax['B'].plot(energies, energies, '-b', linewidth='0.7')
        
        r2 = np.corrcoef(energies, energies_pred)

        ax['B'].set_xlabel('Energies')
        ax['B'].set_ylabel('Energies Predicted')
        ax['B'].set_title(f'r$^2$: {r2[0,1]:.3f}\nRMSE: {PE_e:.2f}    Acc: {(1-PE_e)*100:.2f}%\nFinal_Molecules: {len(energies)}    Outliers_retired: {outliers_count}')

        # # B Err plot
        # b_Err = ((self.database['b_pred'] - self.database['B_opt4'])/self.database['B_opt4'])*100
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
        ax['C'].set_title(f'Error over Energie    MAE = {MAE_e:.4f}\nAbsolute Mean Error: {meanErr:.2f}%    Absolute Max Error: {maxErr:.2f}%')

        # b and Ecorr relation with Ne
        fig_b, ax = plt.subplots(2)
        fig_b.set_size_inches(20, 13)
        fig_b.suptitle('Distribution')

        ax[0].plot(self.database['Ne'], self.database['B_opt4'], '.r')
        ax[0].set_xlabel('Electrons')
        ax[0].set_ylabel('b opt')

        ax[1].plot(self.database['Ne'], self.database['CIe'], '.g')
        ax[1].set_xlabel('Electrons')
        ax[1].set_ylabel('Ecorr')

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

        if show:
            plt.show()

        if save:
            path = os.path.join('predictions', 'plots')
            if not os.path.isdir(path):
                os.makedirs(path)
            
            self.__save_img_outliers(outliers_ID)
            fig.savefig(os.path.join(path, f'a{self.alpha}_{self.b}{self.out}_{self.percent}.pdf'), dpi=450, format='pdf')
            fig_b.savefig(os.path.join(path, f'Distribution_{self.alpha}_{self.b}.pdf'), dpi=450, format='pdf')
            fig_c.savefig(os.path.join(path, f'Ne_Distribution_{self.alpha}_{self.b}.pdf'), dpi=450, format='pdf')
