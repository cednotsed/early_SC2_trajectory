from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

aa_list = ["G", "A", "V", "L", "I",
           "T", "S", "M", "C", "P",
           "F", "Y", "W", "H", "K",
           "R", "D", "E", "N", "Q"]

dict_list = []

for aa in aa_list:
    X = ProteinAnalysis(aa)
    temp = {'amino_acid': aa, 'mw': X.molecular_weight(), 'charge': X.charge_at_pH(7),
                         'hydropathy': X.gravy()}
    dict_list.append(temp)


all_df = pd.DataFrame(dict_list)

all_df.to_csv('/mnt/c/git_repos/early_SC2_trajectory/data/metadata/aa_properties.csv', index = False)


