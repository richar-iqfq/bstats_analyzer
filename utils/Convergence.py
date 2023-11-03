import pandas as pd
import os
import re

'''
Script for checking the convergence of each molecule calculated
'''

class Convergence():
    def __init__(self, database, calc_types=(1,2,3,4,5,6,7,8)):
        self.database = pd.read_csv(database)
        self.availabable_columns = self.database.columns
        self.calc_types = calc_types
        self.not_converged_count = {}
        self.not_converged_percent = {}

    def write_molecules(self, output_path):
        '''
        Start computation and save files.
        '''
        for calc in self.calc_types:
            not_convergence = self.database.loc[self.database[f'Err{calc}'] >= 1e-3]

            file_name = os.path.join('results', 'recovering', output_path, f'c{calc}.csv')
            if len(not_convergence) > 0:
                not_convergence.to_csv(file_name, index=False)

    def get_count(self, specific=None):
        '''
        Obtain the non convergence count for alpha value given

        Parameter
        ---------
        alpha (`float`):
            Alpha value to get count

        Returns
        -------
            A list with the b1-b5 no convergence count
        '''

        if specific:
            count = len(self.database[self.database[f'Err{specific}']>=1E-3])
        
            self.not_converged_count[f'b{specific}'] = count
            self.not_converged_percent[f'b{specific}'] = count*100/len(self.database[f'Err{type}'])

        else:
            for type in self.calc_types:
                if f'Err{type}' in self.availabable_columns:
                    count = len(self.database[self.database[f'Err{type}']>=1E-3])
                    
                    self.not_converged_count[f'b{type}'] = count
                    self.not_converged_percent[f'b{type}'] = count/len(self.database[f'Err{type}'])

                else:
                    self.not_converged_count[f'b{type}'] = 0
                    self.not_converged_percent[f'b{type}'] = 0

        return self.not_converged_count, self.not_converged_percent

    def write_count(self, alpha_list):
        '''
        Writes the count csv file

        Parameter
        ---------
        alpha_list (`list` or `tuple` of `float`)
            alpha values to write in file
        '''
        dict = {
            'alpha' : [],
            'b1' : [],
            'b2' : [],
            'b3' : [],
            'b4' : [],
            'b5' : [],
            'b6' : [],
            'b7' : [],
            'b8' : []
        }

        for alpha in alpha_list:
            count = self.get_count(alpha)
            dict['alpha'].append(alpha)
            for type in range(1, 6):
                dict[f'b{type}'].append(count[type-1])

        df = pd.DataFrame(dict)

        out = os.path.join(self.path, 'Convergence_count_results.csv')
        df.to_csv(out, index=False)

    def get_valid_b(self, max_failed=200):
        count, perc = self.get_count()

        calc_types = []

        for b in self.calc_types:
            if count[f'b{b}'] < max_failed:
                calc_types.append(b)

        return calc_types
