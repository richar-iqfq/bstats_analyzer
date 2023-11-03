import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import numpy as np
import os
from sklearn.metrics import mean_absolute_error
from utils.Convergence import Convergence

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
        Path to save the files.
    
    alpha (`float`):
        Alpha value

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
    def __init__(self):
        self.plot_colors = ['blue', 'green', 'gray', 'orange', 'red',
                            'blue', 'green', 'gray', 'orange', 'red']

    def retire_anomalies(self, data, cut_off=3):
        anomalies = np.ones(len(data), dtype=bool)

        std = np.std(data)
        mean = np.mean(data)
        anomaly_cut_off = std * cut_off

        lower_limit = mean - anomaly_cut_off
        upper_limit = mean + anomaly_cut_off

        for i, value in enumerate(data):
            if value > upper_limit or value < lower_limit:
                anomalies[i] = 0

        data = data[anomalies]
        anomalies_count = np.count_nonzero(np.invert(anomalies))

        return data, anomalies_count

    def plot_dispersion(self, database, alpha, save=False, show=False):
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
        optim_plots_path = os.path.join('results', 'optimization', f'a{alpha}_results')
        
        database_df = pd.read_csv(database)

        available_columns = database_df.columns

        conv = Convergence(database)
        count, perc = conv.get_count()
        b_list = range(1, 9)#conv.get_valid_b()

        # std for bi
        b_std = [0]*8
        b_var = [0]*8
        b_mean = [0]*8
        
        for bi in b_list:
            i = bi - 1
                
            if f'B_opt{bi}' in available_columns:
                b_std[i] = database_df[f'B_opt{bi}'].std(ddof=0)
                b_var[i] = database_df[f'B_opt{bi}'].var(ddof=0)
                b_mean[i] = database_df[f'B_opt{bi}'].mean()

        fig, axis = plt.subplot_mosaic('ABCD;EFGH')
        fig.suptitle(f'Alpha={alpha}')
        fig.subplots_adjust(left=0.04, bottom=0.075, right=0.97, top=0.94)
        fig.set_size_inches(20, 13)

        axes = [axis['A'], axis['B'], axis['C'], axis['D'], axis['E'], axis['F'], axis['G'], axis['H']]
                
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # Creating the plots
        for bi in b_list:
            i = bi - 1
            
            if f'B_opt{bi}' in available_columns:
                data = database_df[f'B_opt{bi}']

                data, anomalies_count = self.retire_anomalies(data)

                center = data.min() + (data.max() - data.min())/2

                b_std = np.std(data)
                b_mean = np.mean(data)
                b_var = np.var(data)

                size = len(data)

                values, bins, bars = axes[i].hist(data, 40, color=self.plot_colors[i], weights=np.ones(size)/size)
                axes[i].yaxis.set_major_formatter(PercentFormatter(1))
                axes[i].axvline(b_mean, color='k', linestyle='dashed', linewidth=1)            

                # Add the text box
                text_str = f'Dispersión\n$\sigma={b_std:.4f}$\n$\sigma^{2}={b_var:.4f}$\nmean={b_mean:.2f}'

                if b_mean > center:
                    x_pos = 0.05
                else:
                    x_pos = 0.65

                axes[i].text(x_pos, 0.95, text_str, fontsize=12, bbox=props, transform=axes[i].transAxes, verticalalignment='top')

                axes[i].set_xlabel(f'B_opt{bi}')
                axes[i].set_title(f"{100*perc[f'b{bi}']:.1f}% Not converged", fontweight ="bold")

        if show:
            plt.show()

        if save:
            path = os.path.join(optim_plots_path, f'b_dispersion_a{alpha}.pdf')
            fig.savefig(path, dpi=450, format='pdf')

        plt.close()

    def plot_correlation_matrix(self, database, output_folder, save=False, show=False):
        '''
        Plot correlation matrix

        save (`bool`):
            if True, saves the fig_A.
            Default is False
        
        show(`bool`):
            if True, shows the builded plot. Default is False
        '''

        database_df = pd.read_csv(database)
        df = database_df.drop(columns=['HF', 'Factor', 'Err1', 'Err2', 'Err3', 'Err4', 'Err5', 'Err6', 'Err7', 'Err8'])

        corr_m = df.corr(numeric_only=True)

        fig, ax = plt.subplots(1)
        fig.suptitle('Correlation Matrix', fontweight ="bold")
        fig.set_size_inches(9, 9)

        sns.heatmap(corr_m, linewidths=0.5, mask=(np.abs(corr_m) <= 0.3), annot=True, annot_kws={"size":6}, square=True, ax=ax)
        
        if show:
            plt.show()

        if save:
            path = os.path.join(self.optim_plots_path, output_folder)
            if not os.path.isdir(path):
                os.makedirs(path)
            
            file = os.path.join(path, 'Correlation_matrix.pdf')
            fig.savefig(file, dpi=450, format='pdf')

        plt.close()

    def plot_recovered_err(self, database, alpha, percent, save=False, show=False):
        '''
        Plot the correlation energy recovering

        Parameters:
        save (`bool`):
            if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
            Default is False
        
        show(`bool`):
            if True, shows the builded plot. Default is False
        '''
        # recover_plots_path = os.path.join('results', 'recovering', f'a{alpha}_results')
        saved_path = os.path.join('results', 'recovering')

        b_amount = 8
        database_df = pd.read_csv(database)

        convergence = Convergence(database)
        valid_b = convergence.get_valid_b()

        Ecorr_real = database_df['CIe']
        Mu_real = database_df['Mu']
        Theta_real = database_df['Theta']
        
        MAE_E_results = {i+1 : {} for i in range(b_amount)}
        MAE_Mu_results = {i+1 : {} for i in range(b_amount)}
        MAE_Theta_results = {i+1 : {} for i in range(b_amount)}
        
        Recover_results = {}

        for perc in percent:

            print(f'\nPercent: {perc}')
            
            for b in valid_b:

                file = os.path.join(saved_path, 'saved_values', f'data_alpha_{alpha}', f'Ecorrb{b}_a{alpha}_perc{perc}.npy')
                array_data = np.load(file)

                Recover_results[f'Ecorr_{b}'] = array_data[:,0]
                Recover_results[f'Theta_{b}']= array_data[:,1]
                Recover_results[f'Mu_{b}'] = array_data[:,2]

            for b in valid_b:

                MAE = mean_absolute_error(Ecorr_real, Recover_results[f'Ecorr_{b}'])
                MAE_E_results[b][perc] = MAE

                MAE = mean_absolute_error(Mu_real, Recover_results[f'Mu_{b}'])
                MAE_Mu_results[b][perc] = MAE

                MAE = mean_absolute_error(Theta_real, Recover_results[f'Theta_{b}'])
                MAE_Theta_results[b][perc] = MAE

        # Build and save the dataframes to csv files
        MAE_E_df = pd.DataFrame(MAE_E_results)
        path = os.path.join(saved_path, 'dataframes', f'MAE_E_a{alpha}.csv')
        MAE_E_df.to_csv(path)

        MAE_Mu_df = pd.DataFrame(MAE_Mu_results)
        path = os.path.join(saved_path, 'dataframes', f'MAE_Mu_a{alpha}.csv')
        MAE_Mu_df.to_csv(path)

        MAE_Theta_df = pd.DataFrame(MAE_Theta_results)
        path = os.path.join(saved_path, 'dataframes', f'MAE_Theta_a{alpha}.csv')
        MAE_Theta_df.to_csv(path)

        # Build and save plots
        #================================= E_corr =========================================
        fig_Ecorr, axis = plt.subplot_mosaic('ABCD;EFGH')
        fig_Ecorr.suptitle(f'Energy Recover\nAlpha = {alpha}', fontweight='bold')
        fig_Ecorr.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Ecorr.set_size_inches(20, 13)

        axes = [
            axis['A'], axis['B'], axis['C'], axis['D'],
            axis['E'], axis['F'], axis['G'], axis['H']
        ]
        
        for b in valid_b:
            x = np.array(percent)
            y = [val for val in MAE_E_results[b].values()]

            y = -np.log10(y)

            axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
            axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
            axes[b-1].set_xlabel('Percentage')
            axes[b-1].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', f'results_{alpha}')

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)
            
            path = os.path.join(plots_path, f'Correlation_Energy_Recover_a{alpha}.pdf')
            fig_Ecorr.savefig(path, dpi=450, format='pdf')

        #================================= Mu_recover =========================================
        fig_Mu, axis = plt.subplot_mosaic('ABCD;EFGH')
        fig_Mu.suptitle(f'Mu Recover\nAlpha = {alpha}', fontweight='bold')
        fig_Mu.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Mu.set_size_inches(20, 13)

        axes = [
            axis['A'], axis['B'], axis['C'], axis['D'],
            axis['E'], axis['F'], axis['G'], axis['H']
        ]
        
        for b in valid_b:
            x = np.array(percent)
            y = [val for val in MAE_Mu_results[b].values()]

            y = -np.log10(y)

            axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
            axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
            axes[b-1].set_xlabel('Percentage')
            axes[b-1].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', f'results_{alpha}')

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)
            
            path = os.path.join(plots_path, f'Mu_Recover_a{alpha}.pdf')
            fig_Mu.savefig(path, dpi=450, format='pdf')

        #================================= Theta_recover =========================================
        fig_Theta, axis = plt.subplot_mosaic('ABCD;EFGH')
        fig_Theta.suptitle(f'Theta Recover\nAlpha = {alpha}', fontweight='bold')
        fig_Theta.subplots_adjust(left=0.06, bottom=0.06, right=0.93, top=0.93, hspace=0.26)
        fig_Theta.set_size_inches(20, 13)

        axes = [
            axis['A'], axis['B'], axis['C'], axis['D'],
            axis['E'], axis['F'], axis['G'], axis['H']
        ]
        
        for b in valid_b:
            x = np.array(percent)
            y = [val for val in MAE_Theta_results[b].values()]

            y = -np.log10(y)

            axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
            axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
            axes[b-1].set_xlabel('Percentage')
            axes[b-1].set_ylabel('-log(MAE)')
        
        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', f'results_{alpha}')

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)
            
            path = os.path.join(plots_path, f'Theta_Recover_a{alpha}.pdf')
            fig_Theta.savefig(path, dpi=450, format='pdf')
        
        plt.close('all')