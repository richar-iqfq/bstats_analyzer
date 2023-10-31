from pypdf import PdfMerger
import os

class PlotMerge():
    def __init__(self):
        self.path = os.path.join('results', 'optimization', 'Statistic_Plots')
        self.merger = PdfMerger()

    def mergeDispersion(self, alpha):
        pdfs = []

        for a in alpha:
            pdf_path = os.path.join('results', 'optimization', f'a{a}_results', f'b_dispersion_a{a}.pdf')
            pdfs.append(pdf_path)

        for pdf in pdfs:
            self.merger.append(pdf)

        output_path = os.path.join(self.path, 'b_dispersion.pdf')
        self.merger.write(output_path)
        self.merger.close()
