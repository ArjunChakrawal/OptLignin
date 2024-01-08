#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.close('all')
import plotly.express as px

# In[3]:


data = pd.read_excel('est_params/exp/fit_results_table.xlsx')
name_counts = data['study_name'].value_counts()


lat = data.loc[~np.isnan(data['Lattitude_N_']), 'Lattitude_N_']
long = data.loc[~np.isnan(data['Longitude_E_']), 'Longitude_E_']


loc = pd.DataFrame({'Latitude': lat, 'Longitude': long})
fig = px.scatter_geo(loc, lat='Latitude', lon='Longitude', color_discrete_sequence=['red'])
fig.update_geos(projection_type="natural earth")  # You can choose a different projection type

fig.update_layout(title_text='Your Title',margin=dict(l=0, r=0, b=0, t=0))
fig.write_html('first_figure.html', auto_open=True)
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
export_params = dict(format='png', width=1200, height=800, scale=3)
fig.write_image("results/FigureS1.png",**export_params)

# In[] 
data_lsq = pd.read_excel('LSQ_fitting/est_params/data summary_LSQ.xlsx')

rename_dict = {'LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_': 'Lig_method',
               'LitterTypeGrass_1Leaf_2Needle_3Roots_4Wood_5Lichen_6': 'LitterType',
               'ClimateTundra_1Boreal_2CoolTemperate_3WarmTemperate_4Tropical_5': 'climate',
               'LigDecStartDay': 'tau', 'AISC0': 'ARC_0','MATC':'MAT','AIC':'BIC','AIC_Co':'BIC_Co','AIC_CT':'BIC_CT'}
data.rename(columns=rename_dict, inplace=True)

data_lsq.rename(columns=rename_dict, inplace=True)

rename_dict = {'aromaticC0': 'ARC_0'}
data.rename(columns=rename_dict, inplace=True)

for i in range(data_lsq.shape[0]):
    id = (data_lsq.iloc[i]['Lig_method'] == 1) or \
         (data_lsq.iloc[i]['Lig_method'] == 2)
    if id:
        data_lsq.at[i, 'ARC_0'] = 0.25 * data_lsq.iloc[i]['ARC_0']

data = data.loc[data['r2_co'] > 0]
data = data.loc[data['r2'] > 0]

data['model'] = 'OCP'
data_lsq = data_lsq.loc[data_lsq['r2_co'] >0]
data_lsq['model'] = 'LSC'

colmns = ['r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT',
          'rmse_CT', 'BIC','BIC_Co','BIC_CT', 'model']
newdata = pd.concat([data[colmns], data_lsq[colmns]])


# In[] Model performance figure
order = None
xlevel = 'model'
hueLevel = None

plt.figure(figsize=(8.5,4))
sns.set_style("white")
plt.clf
plt.subplot(131)
ax = sns.boxplot(data=newdata, y='rmse', x=xlevel, orient='v', palette=sns.color_palette('colorblind'), hue=hueLevel,
                 order=order, linewidth=0.5, width=0.5)
# ax.set_xticklabels([])
plt.ylim([0, 0.5])
plt.xlabel('')
plt.ylabel(r'$\it rmse$')


plt.subplot(132)
ax = sns.boxplot(data=newdata, y='r2', x=xlevel, orient='v', palette=sns.color_palette('colorblind'), hue=hueLevel,
                 order=order, linewidth=0.5, width=0.5)
plt.ylim([0.0, 1])
plt.xlabel('')
plt.ylabel(r'$\it r^2$')

plt.subplot(133)
ax = sns.boxplot(data=newdata, y='BIC', x=xlevel, orient='v', palette=sns.color_palette('colorblind'), hue=hueLevel,
                 order=order, linewidth=0.5, width=0.5)
plt.xlabel('')
plt.ylabel(r'$\it BIC $')
plt.tight_layout()
plt.savefig('results/Figure5.png', dpi=300)
plt.show()

