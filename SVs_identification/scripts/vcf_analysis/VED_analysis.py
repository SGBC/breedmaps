import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np

files = glob.glob("/Users/jj/breedmaps/SVs_identification/results/variant_effect_predictor/VEP*")
dfs = [None] * len(files)
for i in range(len(files)):
    file = files[i]
    name = file.split('/')[-1]
    df = pd.read_csv(
        file,
        sep='\t',
        dtype={'IMPACT': str}
    )
    df.rename(columns={'#Uploaded_variation': 'SV_ID'}, inplace=True)
    df['FILE'] = name
    high_impact = df[df['IMPACT'] == 'HIGH']
    dfs[i] = high_impact
    if i == 0:
        all_df = high_impact
    else:
        all_df = pd.concat([all_df, high_impact])

# new_df = dfs[2][dfs[2]['IMPACT'] == 'HIGH']
# print(new_df)

new_df1 = all_df.groupby(['Consequence']).count()
new_df2 = pd.DataFrame({'COUNT':new_df1.iloc[:,2]})
plot = new_df2.plot.pie(y='COUNT')
plt.show()
