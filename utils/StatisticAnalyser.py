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

        self.__begin_alpha_dispersion()

    def __begin_alpha_dispersion(self, calc_types=8):
        self.not_converged = {
            i+1 : {} for i in range(calc_types)
        }

        self.sigma = {
            i+1 : {} for i in range(calc_types)
        }

        self.sigma_2 = {
            i+1 : {} for i in range(calc_types)
        }

        self.mean = {
            i+1 : {} for i in range(calc_types)
        }

    def retire_not_converged(self, database, type):
        database_df = database[ database[f'Err{type}'] <= 1E-4 ]

        return database_df

    def retire_anomalies(self, dataset, cut_off=3):
        data = dataset.values

        # Remove extreme values
        data = data[ data >= -200 ]
        data = data[ data <= 200 ]

        # Dispersion
        std = np.std(data)
        mean = np.mean(data)
        anomaly_cut_off = std * cut_off

        lower_limit = mean - anomaly_cut_off
        upper_limit = mean + anomaly_cut_off

        # Retire anomalies
        anomalies = np.ones(len(data), dtype=bool)
        for i, value in enumerate(data):
            if value > upper_limit or value < lower_limit:
                anomalies[i] = 0

        data = data[anomalies]
        anomalies_count = np.count_nonzero(np.invert(anomalies))

        return data, anomalies_count

    def plot_dispersion(self, database, alpha, save=False, show=False, append=False, mode=None):
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

        if mode == 'odd':
            self.optim_plots_path = os.path.join('results', 'optimization', f'a{alpha}_odd_results')
        elif mode == 'even':
            self.optim_plots_path = os.path.join('results', 'optimization', f'a{alpha}_even_results')
        else:
            self.optim_plots_path = os.path.join('results', 'optimization', f'a{alpha}_results')

        database_df = pd.read_csv(database)

        available_columns = database_df.columns

        conv = Convergence(database)
        count, perc = conv.get_count()
        b_list = range(1, 9)#conv.get_valid_b()

        # std for bi
        b_std = [0]*8
        b_var = [0]*8
        b_mean = [0]*8

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
                database_df_cleaned = self.retire_not_converged(database_df, bi)

                data = database_df_cleaned[f'B_opt{bi}']

                data, anomalies_count = self.retire_anomalies(data)

                center = data.min() + (data.max() - data.min())/2

                b_std = np.std(data)
                b_var = np.var(data)
                b_mean = np.mean(data)

                if append:
                    self.not_converged[bi][alpha] = count[f'b{bi}']
                    self.sigma[bi][alpha] = b_std
                    self.sigma_2[bi][alpha] = b_var
                    self.mean[bi][alpha] = b_mean

                size = len(data)

                values, bins, bars = axes[i].hist(data, 40, color=self.plot_colors[i], weights=np.ones(size)/size)
                axes[i].yaxis.set_major_formatter(PercentFormatter(1))
                axes[i].axvline(b_mean, color='k', linestyle='dashed', linewidth=1)

                # Add the text box
                text_str = f'Dispersion\n$\sigma={b_std:.4f}$\n$\sigma^{2}={b_var:.4f}$\nmean={b_mean:.2f}\nanomalies={anomalies_count}'

                if b_mean > center:
                    x_pos = 0.05
                else:
                    x_pos = 0.65

                axes[i].text(x_pos, 0.95, text_str, fontsize=12, bbox=props, transform=axes[i].transAxes, verticalalignment='top')

                # Add bar values
                y_height = np.zeros(40)
                for j, val in enumerate(values):
                    y_height[j] = val

                rects = axes[i].patches
                x_center = np.zeros(len(rects))
                for j, rect in enumerate(rects):
                    x_center[j] = rect.get_x() + rect.get_width()/2

                for j in range(len(x_center)):
                    if values[j] != 0:
                        axes[i].text(x_center[j], y_height[j]+0.001, '▴', fontdict={'size': 6})

                axes[i].set_xlabel(f'B_opt{bi}')
                axes[i].set_title(f"{100*perc[f'b{bi}']:.1f}% Not converged", fontweight ="bold")

            else:
                if append:
                    self.not_converged[bi][alpha] = None
                    self.sigma[bi][alpha] = None
                    self.sigma_2[bi][alpha] = None
                    self.mean[bi][alpha] = None

        if show:
            plt.show()

        if save:
            path = os.path.join(self.optim_plots_path, f'b_dispersion_a{alpha}.pdf')
            fig.savefig(path, dpi=450, format='pdf')

        plt.close()

    def plot_general_dispersion_stats(self, save=True, show=False, calc_types=8, mode=None):
        if mode == 'odd':
            output_path = os.path.join('results', 'optimization', 'Statistic_Plots', 'general', 'by_oddmetric')
        elif mode == 'even':
            output_path = os.path.join('results', 'optimization', 'Statistic_Plots', 'general', 'by_evenmetric')
        else:
            output_path = os.path.join('results', 'optimization', 'Statistic_Plots', 'general', 'by_metric')

        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        ##################################################################################
        # Convergence count
        fig_a, ax_a = plt.subplots(2, 2)
        ax_a = ax_a.flatten()
        fig_a.suptitle(f'Not-Convergence count 1/2')
        fig_a.set_size_inches(20, 13)
        fig_a.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        fig_b, ax_b = plt.subplots(2, 2)
        ax_b = ax_b.flatten()
        fig_b.suptitle(f'Not-Convergence count 2/2')
        fig_b.set_size_inches(20, 13)
        fig_b.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        for i in range(calc_types):
            bi = i+1

            # Not converged count
            x_count = self.not_converged[bi].keys()
            y_count = self.not_converged[bi].values()

            if i<4:
                ax_a[i].set_title(f'b{bi}')
                ax_a[i].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_a[i].set_xlabel('alpha')
                ax_a[i].set_ylabel('Not-Convergence Count')
            else:
                j = i-4
                ax_b[j].set_title(f'b{bi}')
                ax_b[j].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_b[j].set_xlabel('alpha')
                ax_b[j].set_ylabel('Not-Convergence Count')

        if save:
            file_a = os.path.join(output_path, 'Not-Convergence_Count_1.pdf')
            fig_a.savefig(file_a, dpi=450, format='pdf')

            file_b = os.path.join(output_path, 'Not-Convergence_Count_2.pdf')
            fig_b.savefig(file_b, dpi=450, format='pdf')

        ##################################################################################
        # Sigma
        fig_a, ax_a = plt.subplots(2, 2)
        ax_a = ax_a.flatten()
        fig_a.suptitle(f'Standard Deviation 1/2')
        fig_a.set_size_inches(20, 13)
        fig_a.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        fig_b, ax_b = plt.subplots(2, 2)
        ax_b = ax_b.flatten()
        fig_b.suptitle(f'Standard Deviation 2/2')
        fig_b.set_size_inches(20, 13)
        fig_b.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        for i in range(calc_types):
            bi = i+1

            # Sigma
            x_count = self.sigma[bi].keys()
            y_count = self.sigma[bi].values()

            if i<4:
                ax_a[i].set_title(f'b{bi}')
                ax_a[i].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_a[i].set_xlabel('alpha')
                ax_a[i].set_ylabel('$\sigma$')
            else:
                j = i-4
                ax_b[j].set_title(f'b{bi}')
                ax_b[j].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_b[j].set_xlabel('alpha')
                ax_b[j].set_ylabel('$\sigma$')

        if save:
            file_a = os.path.join(output_path, 'Standard_Deviation_1.pdf')
            fig_a.savefig(file_a, dpi=450, format='pdf')

            file_b = os.path.join(output_path, 'Standard_Deviation_2.pdf')
            fig_b.savefig(file_b, dpi=450, format='pdf')

        ##################################################################################
        # Sigma²
        fig_a, ax_a = plt.subplots(2, 2)
        ax_a = ax_a.flatten()
        fig_a.suptitle(f'Variance 1/2')
        fig_a.set_size_inches(20, 13)
        fig_a.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        fig_b, ax_b = plt.subplots(2, 2)
        ax_b = ax_b.flatten()
        fig_b.suptitle(f'Variance 2/2')
        fig_b.set_size_inches(20, 13)
        fig_b.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        for i in range(calc_types):
            bi = i+1

            # Sigma²
            x_count = self.sigma_2[bi].keys()
            y_count = self.sigma_2[bi].values()

            if i<4:
                ax_a[i].set_title(f'b{bi}')
                ax_a[i].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_a[i].set_xlabel('alpha')
                ax_a[i].set_ylabel('$\sigma^{2}$')
            else:
                j = i-4
                ax_b[j].set_title(f'b{bi}')
                ax_b[j].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_b[j].set_xlabel('alpha')
                ax_b[j].set_ylabel('$\sigma^{2}$')

        if save:
            file_a = os.path.join(output_path, 'Variance_1.pdf')
            fig_a.savefig(file_a, dpi=450, format='pdf')

            file_b = os.path.join(output_path, 'Variance_2.pdf')
            fig_b.savefig(file_b, dpi=450, format='pdf')

        ##################################################################################
        # mean
        fig_a, ax_a = plt.subplots(2, 2)
        ax_a = ax_a.flatten()
        fig_a.suptitle(f'Mean 1/2')
        fig_a.set_size_inches(20, 13)
        fig_a.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        fig_b, ax_b = plt.subplots(2, 2)
        ax_b = ax_b.flatten()
        fig_b.suptitle(f'Mean 2/2')
        fig_b.set_size_inches(20, 13)
        fig_b.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)

        for i in range(calc_types):
            bi = i+1

            # Mean
            x_count = self.mean[bi].keys()
            y_count = self.mean[bi].values()

            if i<4:
                ax_a[i].set_title(f'b{bi}')
                ax_a[i].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_a[i].set_xlabel('alpha')
                ax_a[i].set_ylabel('$\\bar{X}$')
            else:
                j = i-4
                ax_b[j].set_title(f'b{bi}')
                ax_b[j].plot(x_count, y_count, color=self.plot_colors[i], marker='.')
                ax_b[j].set_xlabel('alpha')
                ax_b[j].set_ylabel('$\\bar{X}$')

        if save:
            file_a = os.path.join(output_path, 'Mean_1.pdf')
            fig_a.savefig(file_a, dpi=450, format='pdf')

            file_b = os.path.join(output_path, 'Mean_2.pdf')
            fig_b.savefig(file_b, dpi=450, format='pdf')

        if show:
            plt.show()

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
        # df = database_df.drop(columns=['HF', 'Factor', 'Err1', 'Err2', 'Err3', 'Err4', 'Err5', 'Err6', 'Err7', 'Err8'])
        df = database_df.drop(columns=['HF', 'Factor'])

        corr_m = df.corr(numeric_only=True)

        fig, ax = plt.subplots(1)
        fig.suptitle('Correlation Matrix', fontweight ="bold")
        fig.set_size_inches(9, 9)

        sns.heatmap(corr_m, linewidths=0.5, mask=(np.abs(corr_m) <= 0.3), annot=True, annot_kws={"size":7}, square=True, ax=ax)

        if show:
            plt.show()

        if save:
            if not os.path.isdir(output_folder):
                os.makedirs(output_folder)

            file_pdf = os.path.join(output_folder, 'Correlation_matrix.pdf')
            file_png = os.path.join(output_folder, 'Correlation_matrix.png')

            fig.savefig(file_pdf, dpi=450, format='pdf')
            fig.savefig(file_png, dpi=450, format='png')

        plt.close()

    def plot_recovered_err(self, database, alpha, percent, save=False, show=False, mode=''):
        '''
        Plot the correlation energy recovering

        Parameters:
        save (`bool`):
            if True, saves the fig to root/results/recovering/Statistic_Plots/output_path/fig.pdf.
            Default is False

        show(`bool`):
            if True, shows the builded plot. Default is False
        '''
        # container per mode
        if mode == 'even':
            folder = 'dataframes_even'
            container = f'results_a{alpha}_even'
        elif mode == 'odd':
            folder = 'dataframes_odd'
            container = f'results_a{alpha}_odd'
        else:
            folder = 'dataframes'
            container = f'results_a{alpha}'

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
            for b in valid_b:
                if mode == 'even':
                    file = os.path.join(saved_path, 'saved_values', f'data_alpha_{alpha}_even', f'Ecorrb{b}_a{alpha}_perc{perc}.npy')
                elif mode == 'odd':
                    file = os.path.join(saved_path, 'saved_values', f'data_alpha_{alpha}_odd', f'Ecorrb{b}_a{alpha}_perc{perc}.npy')
                else:
                    file = os.path.join(saved_path, 'saved_values', f'data_alpha_{alpha}', f'Ecorrb{b}_a{alpha}_perc{perc}.npy')

                try:
                    array_data = np.load(file)

                    Recover_results[f'Ecorr_{b}'] = array_data[:,0]
                    Recover_results[f'Theta_{b}']= array_data[:,1]
                    Recover_results[f'Mu_{b}'] = array_data[:,2]
                except:
                    Recover_results[f'Ecorr_{b}'] = None
                    Recover_results[f'Theta_{b}']= None
                    Recover_results[f'Mu_{b}'] = None

            for b in valid_b:
                try:
                    MAE = mean_absolute_error(Ecorr_real, Recover_results[f'Ecorr_{b}'])
                    MAE_E_results[b][perc] = MAE

                    MAE = mean_absolute_error(Mu_real, Recover_results[f'Mu_{b}'])
                    MAE_Mu_results[b][perc] = MAE

                    MAE = mean_absolute_error(Theta_real, Recover_results[f'Theta_{b}'])
                    MAE_Theta_results[b][perc] = MAE
                except:
                    MAE_E_results[b][perc] = None
                    MAE_Mu_results[b][perc] = None
                    MAE_Theta_results[b][perc] = None

        # Build and save the dataframes to csv files                
        MAE_E_df = pd.DataFrame(MAE_E_results)
        path = os.path.join(saved_path, folder, f'MAE_E_a{alpha}.csv')
        MAE_E_df.to_csv(path)

        MAE_Mu_df = pd.DataFrame(MAE_Mu_results)
        path = os.path.join(saved_path, folder, f'MAE_Mu_a{alpha}.csv')
        MAE_Mu_df.to_csv(path)

        MAE_Theta_df = pd.DataFrame(MAE_Theta_results)
        path = os.path.join(saved_path, folder, f'MAE_Theta_a{alpha}.csv')
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
            try:
                x = np.array(percent)
                y = [val for val in MAE_E_results[b].values()]

                y = -np.log10(y)

                axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
                axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
                axes[b-1].set_xlabel('Percentage')
                axes[b-1].set_ylabel('-log(MAE)')
            except:
                pass

        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', container)

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)

            path = os.path.join(plots_path, f'Correlation_{mode}Energy_Recover_a{alpha}.pdf')
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
            try:
                x = np.array(percent)
                y = [val for val in MAE_Mu_results[b].values()]

                y = -np.log10(y)

                axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
                axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
                axes[b-1].set_xlabel('Percentage')
                axes[b-1].set_ylabel('-log(MAE)')
            except:
                pass

        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', container)

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)

            path = os.path.join(plots_path, f'Mu_{mode}Recover_a{alpha}.pdf')
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
            try:
                x = np.array(percent)
                y = [val for val in MAE_Theta_results[b].values()]

                y = -np.log10(y)

                axes[b-1].plot(x, y, marker='+', color=self.plot_colors[b])
                axes[b-1].set_title(f'B_opt{b}', fontweight='bold')
                axes[b-1].set_xlabel('Percentage')
                axes[b-1].set_ylabel('-log(MAE)')
            except:
                pass

        if show:
            plt.show()

        if save:
            plots_path = os.path.join(saved_path, 'Statistic_Plots', container)

            if not os.path.isdir(plots_path):
                os.makedirs(plots_path)

            path = os.path.join(plots_path, f'Theta_{mode}Recover_a{alpha}.pdf')
            fig_Theta.savefig(path, dpi=450, format='pdf')

        plt.close('all')

    def plot_general_recovering_stats(self, alpha_list, percent, save=True, show=False, mode=''):
        if mode=='even':
            path = os.path.join('results', 'recovering', 'dataframes_even')
        elif mode=='odd':
            path = os.path.join('results', 'recovering', 'dataframes_odd')
        else:
            path = os.path.join('results', 'recovering', 'dataframes')
        
        output_path = os.path.join('results', 'recovering', 'Statistic_Plots', 'general', f'{mode}partial')

        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        markers = [
            ".", "^", "v", "^",
            "<", ">", "s", "*",
            "x"
        ]

        # lines = [
        #     "solid", "dotted", "solid", "dashed",
        #     "solid", "dashdot", "solid", "dotted",
        #     "solid"
        # ]

        lines = [
            "solid", "solid", "solid", "solid",
            "solid", "solid", "solid", "solid",
            "solid"
        ]

        min_alpha = min(alpha_list)
        max_alpha = max(alpha_list)

        files = {}

        for alpha in alpha_list:
            files[alpha] = os.path.join(path, f'MAE_E_a{alpha}.csv')

        dataframes = {}

        for alpha in alpha_list:
            dataframes[alpha] = pd.read_csv(files[alpha])

        data = {
            i+1 : {perc : {} for perc in percent} for i in range(8)
        }

        for i in range(8):
            bi = i+1
            for j, perc in enumerate(percent):
                for alpha in alpha_list:
                    df = dataframes[alpha]
                    column = df[f'{bi}']
                    data[bi][perc][alpha] = column[j]

        fig_a, ax_a = plt.subplots(2,2)
        fig_a.suptitle('Energy Recover 1/2')
        fig_a.set_size_inches(20, 13)
        fig_a.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)
        ax_a = ax_a.flatten()

        fig_b, ax_b = plt.subplots(2,2)
        fig_b.suptitle('Energy Recover 2/2')
        fig_b.set_size_inches(20, 13)
        fig_b.subplots_adjust(left=0.06, bottom=0.06, right=0.96, top=0.94, hspace=0.19, wspace=0.15)
        ax_b = ax_b.flatten()

        for i in range(8):
            bi = i+1

            if i<4:
                percent_keys = data[bi].keys()

                for k, perc in enumerate(percent_keys):
                    x = data[bi][perc].keys()
                    y = data[bi][perc].values()

                    y = [-np.log10(val) for val in y]

                    ax_a[i].plot(x, y, marker=markers[k], color=self.plot_colors[i], linestyle=lines[k])
                    
                    ax_a[i].set_title(f'b{bi}')
                    ax_a[i].set_ylabel('-log(MAE)')
                    ax_a[i].set_xlabel('alpha')

                ax_a[i].legend(percent_keys)
                ax_a[i].set_xlim([min_alpha, max_alpha])
                ax_a[i].axhline(1.5, linestyle='dashed', linewidth='0.7')
            else:
                j = i-4
                percent_keys = data[bi].keys()

                for k, perc in enumerate(percent_keys):
                    x = data[bi][perc].keys()
                    y = data[bi][perc].values()

                    y = [-np.log10(val) for val in y]
                    
                    ax_b[j].plot(x, y, marker=markers[k], color=self.plot_colors[i], linestyle=lines[k])
                    ax_b[j].set_title(f'b{bi}')
                    ax_b[j].set_ylabel('-log(MAE)')
                    ax_b[j].set_xlabel('alpha')

                ax_b[j].legend(percent_keys)
                ax_b[j].set_xlim([min_alpha, max_alpha])
                ax_b[j].axhline(1.5, linestyle='dashed', linewidth='0.7')

        if show:
            plt.show()

        if save:
            fig_a.savefig(os.path.join(output_path, f'recovering{mode}_1.pdf'), dpi=450, format='pdf')
            fig_b.savefig(os.path.join(output_path, f'recovering{mode}_2.pdf'), dpi=450, format='pdf')
