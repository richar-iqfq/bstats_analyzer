import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from utils.CalculatorAsyncV6 import MainCalculator
from sklearn.metrics import mean_absolute_error
from utils.Convergence import Counter

class StatisticAnalyser():
    '''
    Statistic analyzer for correlation energy recovering

    Parameters
    ----------
    database (`str`):
        Csv database with optimized values (b, theta, mu) for analisys

    energies (`str`):
        Feather database with energie values.

    output_path (`str`):
        Path to save the .npy files. root/results/recovering/saved_values/output_path
    
    alpha (`float`):
        Alpha value

    percent (`tuple`) of (`float`):
        Percent values to compute the recovering error of the form -> (0.9, 0.95, 1, 1.05, 1.1)

    Methods
    -------
    plot_dispersion(bars=False, save=False, show=True)
    Plot the b dispersion graphic

    bars (`bool`):
        if True, show label with count over each ploted bar
    
    save (`bool`):
        if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
        Default is False
    
    show(`bool`):
        if True, shows the builded plot. Default is False

    plot_correlation_matrix(save=False, show=True)
    Plot correlation matrix

    save (`bool`):
        if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
        Default is False
    
    show(`bool`):
        if True, shows the builded plot. Default is False

    plot_err(cores='All', save=False, show=False, tolerance=None)
    Plot the correlation energy recovering

    cores ('All') or (`int`):
        Number of cores/threads to use in calculation.

    save (`bool`):
        if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
        Default is False
    
    show(`bool`):
        if True, shows the builded plot. Default is False

    tolerance (`dict`):
        Dictionary with the tolerance values. Keys->'Ne_calc','Energy_S','Energy_Savrg','Recover'
    '''
    def __init__(self, database, energies, output_path, alpha, percent):
        self.database_df = pd.read_csv(database)

        self.database = database
        self.energies = energies
        self.plot_colors = ['blue', 'green', 'gray', 'orange', 'red']
        self.alpha = alpha
        self.percent = percent

        counter = Counter()
        self.count = counter.get_count(alpha)
        
        root = os.getcwd()
        self.output_path = output_path
        self.plots_path = os.path.join(root, 'results', 'recovering', 'Statistic_Plots', self.output_path)
        self.saved_path = os.path.join(root, 'results', 'recovering', 'saved_values', self.output_path)

        if not os.path.isdir(self.plots_path):
            os.makedirs(self.plots_path)

    def plot_dispersion(self, bars=False, save=False, show=False):
        '''
        Plot the b dispersion graphic

        Parameters
        ----------
        bars (`bool`):
            if True, show label with count over each ploted bar.
        
        save (`bool`):
            if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
            Default is False
        
        show(`bool`):
            if True, shows the builded plot. Defaul is False
        '''
        b_list = range(1, 5+1)

        
        # std for bi
        b_std = []
        for bi in b_list:
            b_std.append(self.database_df[f'B_opt{bi}'].std(ddof=0))

        b_var = []
        for bi in b_list:
            b_var.append(self.database_df[f'B_opt{bi}'].var(ddof=0))

        fig, axis = plt.subplot_mosaic('AABB;CCDD;.EE.')
        fig.subplots_adjust(left=0.06, bottom=0.05, right=0.97, top=0.97)
        fig.set_size_inches(20, 13)

        axes = [axis['A'], axis['B'], axis['C'], axis['D'], axis['E']]

        # Creating the plots
        for bi in (1, 2, 3, 4, 5):
            i = bi - 1
            bi = sns.histplot(self.database_df, x=f'B_opt{bi}', stat='count', kde=True, ax=axes[i], color=self.plot_colors[i])
            
            if bars:
                bi_labels = [str(v) if v else '' for v in bi.containers[0].datavalues]
                bi.bar_label(bi.containers[0], labels=bi_labels, rotation=45, fontsize=8)

        # Add the text box
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        for bi in b_list:
            i = bi - 1

            text_str = f'Dispersión\n$\sigma={b_std[i]:.4f}$\n$\sigma^{2}={b_var[i]:.4f}$'
            axes[i].text(0.05, 0.95, text_str, fontsize=12, bbox=props, transform=axes[i].transAxes, verticalalignment='top')
            axes[i].set_xlabel('')
            axes[i].set_title(f'B_opt{bi} -> {self.count[i]} Not converged', fontweight ="bold")

        if show:
            plt.show()

        if save:
            path = os.path.join(self.plots_path, f'b_dispersion_a{self.alpha}.pdf')
            fig.savefig(path, dpi=450, format='pdf')

        plt.close()

    def plot_correlation_matrix(self, save=False, show=False):
        '''
        Plot correlation matrix

        save (`bool`):
            if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
            Default is False
        
        show(`bool`):
            if True, shows the builded plot. Default is False
        '''

        df = self.database_df.drop(columns=['CIe', 'HF', 'Factor'])

        corr_m = df.corr(numeric_only=True)

        fig, ax = plt.subplots(1)
        fig.suptitle('Correlation Matrix', fontweight ="bold")
        fig.set_size_inches(9, 9)

        sns.heatmap(corr_m, linewidths=0.5, mask=(np.abs(corr_m) <= 0.3), annot=True, annot_kws={"size":5}, square=True, ax=ax)
        
        if show:
            plt.show()

        if save:
            path = os.path.join(self.plots_path, f'Correlation_matrix_a{self.alpha}.pdf')
            fig.savefig(path, dpi=450, format='pdf')

        plt.close()

    def plot_Err(self, cores='All', save=False, show=False, tolerance=None, isdata=False):
        '''
        Plot the correlation energy recovering

        Parameters:
        cores ('All') or (`int`):
            Number of cores/threads to use in calculation.

        save (`bool`):
            if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
            Default is False
        
        show(`bool`):
            if True, shows the builded plot. Default is False

        tolerance (`dict`):
            Dictionary with the tolerance values. Keys->'Ne_calc','Energy_S','Energy_Savrg','Recover'
            Default is None

        isdata (`bool`):
            if True, will search the .npy array for the values inside root/results/optimization/
            a{alpha}_results/saved_data
        '''
        b_amount=5
        Ecorr_real = self.database_df['CIe']
        Mu_real = self.database_df['Mu']
        Theta_real = self.database_df['Theta']
        
        MAE_E_results = {i+1 : {} for i in range(b_amount)}
        MAE_Mu_results = {i+1 : {} for i in range(b_amount)}
        MAE_Theta_results = {i+1 : {} for i in range(b_amount)}
        
        if isdata:
            Recover_results = {}

        print('########### Running ##########')
        for perc in self.percent:

            print(f'\nPercent: {perc}')

            if isdata:
                for i in range(b_amount):

                    file = os.path.join(self.saved_path, f'Ecorrb{i+1}_a{self.alpha}_perc{perc}.npy')
                    array_data = np.load(file)

                    Recover_results[f'Ecorr_{i+1}'] = array_data[:,0]
                    Recover_results[f'Theta_{i+1}']= array_data[:,1]
                    Recover_results[f'Mu_{i+1}'] = array_data[:,2]

            else:
                mc = MainCalculator(self.output_path, cores=cores, tolerance=tolerance)
                Recover_results = mc.run_recovering(self.database, self.energies, perc, self.alpha)

            for i in range(b_amount):

                MAE = mean_absolute_error(Ecorr_real, Recover_results[f'Ecorr_{i+1}'])
                MAE_E_results[i+1][perc] = MAE

                MAE = mean_absolute_error(Mu_real, Recover_results[f'Mu_{i+1}'])
                MAE_Mu_results[i+1][perc] = MAE

                MAE = mean_absolute_error(Theta_real, Recover_results[f'Theta_{i+1}'])
                MAE_Theta_results[i+1][perc] = MAE

        # Build and save the dataframes to csv files
        MAE_E_df = pd.DataFrame(MAE_E_results)
        path = os.path.join(self.saved_path, f'MAE_E_a{self.alpha}.csv')
        MAE_E_df.to_csv(path)

        MAE_Mu_df = pd.DataFrame(MAE_Mu_results)
        path = os.path.join(self.saved_path, f'MAE_Mu_a{self.alpha}.csv')
        MAE_Mu_df.to_csv(path)

        MAE_Theta_df = pd.DataFrame(MAE_Theta_results)
        path = os.path.join(self.saved_path, f'MAE_Theta_a{self.alpha}.csv')
        MAE_Theta_df.to_csv(path)

        # Build and save plots
        #================================= E_corr =========================================
        fig_Ecorr, axis = plt.subplot_mosaic('AABB;CCDD;.EE.')
        fig_Ecorr.suptitle('Correlation Energy Recover', fontweight='bold')
        fig_Ecorr.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Ecorr.set_size_inches(20, 13)

        axes = [axis['A'], axis['B'], axis['C'], axis['D'], axis['E']]
        
        for j in range(b_amount):
            x = np.array(self.percent)*100
            y = [val for val in MAE_E_results[j+1].values()]

            y = -np.log10(y)

            axes[j].plot(x, y, marker='+', color=self.plot_colors[j])
            axes[j].set_title(f'B_opt{j+1}', fontweight='bold')
            axes[j].set_xlabel('Percentage')
            axes[j].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            path = os.path.join(self.plots_path, f'Correlation_Energy_Recover_a{self.alpha}.pdf')
            fig_Ecorr.savefig(path, dpi=450, format='pdf')

        #================================= Mu_recover =========================================
        fig_Mu, axis = plt.subplot_mosaic('AABB;CCDD;.EE.')
        fig_Mu.suptitle('Mu Recover', fontweight='bold')
        fig_Mu.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Mu.set_size_inches(20, 13)

        axes = [axis['A'], axis['B'], axis['C'], axis['D'], axis['E']]
        
        for j in range(b_amount):
            x = np.array(self.percent)*100
            y = [val for val in MAE_Mu_results[j+1].values()]

            y = -np.log10(y)

            axes[j].plot(x, y, marker='+', color=self.plot_colors[j])
            axes[j].set_title(f'B_opt{j+1}', fontweight='bold')
            axes[j].set_xlabel('Percentage')
            axes[j].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            path = os.path.join(self.plots_path, f'Mu_Recover_a{self.alpha}.pdf')
            fig_Mu.savefig(path, dpi=450, format='pdf')

        #================================= Theta_recover =========================================
        fig_Theta, axis = plt.subplot_mosaic('AABB;CCDD;.EE.')
        fig_Theta.suptitle('Theta Recover', fontweight='bold')
        fig_Theta.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Theta.set_size_inches(20, 13)

        axes = [axis['A'], axis['B'], axis['C'], axis['D'], axis['E']]
        
        for j in range(b_amount):
            x = np.array(self.percent)*100
            y = [val for val in MAE_Theta_results[j+1].values()]

            y = -np.log10(y)

            axes[j].plot(x, y, marker='+', color=self.plot_colors[j])
            axes[j].set_title(f'B_opt{j+1}', fontweight='bold')
            axes[j].set_xlabel('Percentage')
            axes[j].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            path = os.path.join(self.plots_path, f'Theta_Recover_a{self.alpha}.pdf')
            fig_Theta.savefig(path, dpi=450, format='pdf')
        
        plt.close('all')