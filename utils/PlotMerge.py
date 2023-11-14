from pypdf import PdfMerger
import os

class PlotMerge():
    def __init__(self):
        self.path = os.path.join('results', 'optimization', 'Statistic_Plots')

    def mergeFiles(self, input, output):
        pdfs = os.listdir(input)
        pdfs.sort()

        pdfs = [os.path.join(input, pdf) for pdf in pdfs if '.pdf' in pdf]

        merger = PdfMerger()

        for pdf in pdfs:
            merger.append(pdf)

        merger.write(output)
        merger.close()

    def mergeDispersion(self, alpha, mode=None):
        pdfs = []
        merger = PdfMerger()
        
        for a in alpha:
            if mode == 'odd':
                pdf_path = os.path.join('results', 'optimization', f'a{a}_odd_results', f'b_dispersion_a{a}.pdf')
            elif mode == 'even':
                pdf_path = os.path.join('results', 'optimization', f'a{a}_even_results', f'b_dispersion_a{a}.pdf')
            else:
                pdf_path = os.path.join('results', 'optimization', f'a{a}_results', f'b_dispersion_a{a}.pdf')
            
            pdfs.append(pdf_path)

        for pdf in pdfs:
            merger.append(pdf)

        output_path = os.path.join(self.path, 'general')
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        if mode == 'odd':
            file = os.path.join(output_path, 'b_odd_dispersion.pdf')
        elif mode == 'even':
            file = os.path.join(output_path, 'b_even_dispersion.pdf')
        else:
            file = os.path.join(output_path, 'b_dispersion.pdf')
        merger.write(file)
        merger.close()

    def mergeRecover(self, alpha):
        energy_pdfs = []
        mu_pdfs = []
        theta_pdfs = []
        
        for a in alpha:
            
            pdfs_path = os.path.join('results', 'recovering', 'Statistic_Plots', f'results_{a}')
            output_path = os.path.join('results', 'recovering', 'Statistic_Plots', 'general')

            energy_path = os.path.join(pdfs_path, f'Correlation_Energy_Recover_a{a}.pdf')
            mu_path = os.path.join(pdfs_path, f'Mu_Recover_a{a}.pdf')
            theta_path = os.path.join(pdfs_path, f'Theta_Recover_a{a}.pdf')

            energy_pdfs.append(energy_path)
            mu_pdfs.append(mu_path)
            theta_pdfs.append(theta_path)

            # Energy
            merger = PdfMerger()

            for pdf in energy_pdfs:
                merger.append(pdf)

            file = os.path.join(output_path, 'Energy_recovering.pdf')
            merger.write(file)
            merger.close()

            # Mu
            merger = PdfMerger()

            for pdf in mu_pdfs:
                merger.append(pdf)

            file = os.path.join(output_path, 'Mu_recovering.pdf')
            merger.write(file)
            merger.close()

            # Theta
            merger = PdfMerger()

            for pdf in theta_pdfs:
                merger.append(pdf)

            file = os.path.join(output_path, 'Theta_recovering.pdf')
            merger.write(file)
            merger.close()