import pandas as pd
import os
import re

'''
Script for checking the convergence of each molecule calculated
'''

class Writter():
    '''
    Search if there are come molecules with no convergence

    Parameters
    ----------
    database (`str`):
        database path with the required info (Err1, Err2, etc)

    calc_types (`tuple`) of (`int`):
        Type of calculations to compute (1 to 5), default is (1,2,3,4,5)
    
    Methods
    -------
    run_check()
        start computation and save files.
    '''
    def __init__(self, database, calc_types=(1,2,3,4,5,6,7)):
        self.database_name = database[0:-4]
        self.database = pd.read_csv(database)
        self.calc_types = calc_types

    def run_write(self):
        '''
        Start computation and save files.
        '''
        for calc in self.calc_types:
            not_convergence = self.database.loc[self.database[f'Err{calc}'] >= 1e-3]

            file_name = self.database_name + f'c{calc}.csv'
            if len(not_convergence) > 0:
                not_convergence.to_csv(file_name, index=False)

class Analyzer():
    '''
    Check database results to find errors in convergence

    Parameters
    ----------
    re_path (`str`):
        Re type path to search recurrence. Example >> r'a-[0-9.]+_results'
    
    criteria (`float`):
        Minimum value for convergence. Default is 1E-3

    Methods
    -------
    analyze()
        Run analysis
    '''
    def __init__(self, re_path, criteria=1E-3):
        self.criteria = criteria
        self.files = []
        self.alpha = []

        root = os.getcwd()
        path = os.path.join(root, 'results', 'optimization')
        content = os.listdir(path)

        folders = [folder for folder in content if re.search(re_path, folder)]

        for folder in folders:
            alpha = re.search(r'a(-[0-9.]+)', folder).group(1)
            file = os.path.join(path, folder, f'results_a{alpha}.csv')

            if os.path.isfile(file):
                self.files.append(file)
                self.alpha.append(float(alpha))

    def analyze(self):
        '''
        Run analysis

        Returns
        -------
        List of tuples of type (alpha, type, mean)
        '''
        results = []

        for i, file in enumerate(self.files):
            df = pd.read_csv(file)
            alpha = self.alpha[i]

            for type in (1, 2, 3, 4, 5):
                count = len(df[df[f'Err{type}']>=self.criteria])
                
                if count > len(df)*0.05:
                    correct = df[df[f'Err{type}']<=self.criteria]
                    
                    b_mean = correct[f'B_opt{type}'].mean()
                    
                    results.append((alpha, type, b_mean))
        
        return results
    
class Counter():
    '''
    Count the number of molecules with no-convergence in alpha database
    
    Methods
    -------
    get_count(alpha):
        return a list with the b1-b5 no convergence count

    write_count(alpha_list):
        save a csv file with the count for all the alpha folders given
    '''
    def __init__(self, specific=None):
        root = os.getcwd()
        self.path = os.path.join(root, 'results', 'optimization')
        self.specific = specific

    def get_count(self, alpha):
        '''
        Obtain the no convergence count for alpha given

        Parameter
        ---------
        alpha (`float`):
            alpha value to check count

        Returns
        -------
            a list with the b1-b5 no convergence count
        '''
        folder = os.path.join(self.path, f'a{alpha}_results')
        
        file = os.path.join(folder, f'results_a{alpha}.csv')
        df = pd.read_csv(file)
        
        count = []

        if self.specific:
            count.append(len(df[df[f'Err{self.specific}']>=1E-3]))
        
        else:
            for type in range(1, 6):
                count.append(len(df[df[f'Err{type}']>=1E-3]))

        return count
    
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
        }

        for alpha in alpha_list:
            count = self.get_count(alpha)
            dict['alpha'].append(alpha)
            for type in range(1, 6):
                dict[f'b{type}'].append(count[type-1])

        df = pd.DataFrame(dict)

        out = os.path.join(self.path, 'Convergence_count_results.csv')
        df.to_csv(out, index=False)