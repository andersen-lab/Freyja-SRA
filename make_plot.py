import dash
import pandas as pd
import pickle
import json
import requests
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from scipy import signal
from datetime import date,timedelta
import yaml

# df_meta = pd.read_csv('data/all_metadata.csv',index_col=0)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

with open( "plot_config.yml", "r" ) as f :
    plot_config = yaml.load( f, Loader=yaml.FullLoader )

#--- borrowed from SEARCH wastewater surveillance dashboard, coordinated by Scripps Research.--#
def convert_rbg_to_tuple( rgb ):
    rgb = rgb.lstrip( "#" )
    return tuple( int( rgb[i :i + 2], 16 ) for i in (0, 2, 4) )
def convert_tuple_to_rgb( r, g, b ):
    return '#%02x%02x%02x' % (int(r), int(g), int(b))
def lighten_field( value, alpha, gamma=2.2 ):
    return pow( pow(255, gamma) * (1 - alpha) + pow( value, gamma ) * alpha, 1 / gamma)
def lighten_color( r, g, b, alpha, gamma=2.2 ):
    return lighten_field(r, alpha, gamma ), lighten_field( g, alpha, gamma ), lighten_field( b, alpha, gamma )

children_dict = dict()
delta = 0.15
for key in reversed( list( plot_config.keys() ) ):
    for value in ["name", "members"]:
        assert value in plot_config[key], f"YAML entry {key} is not complete. Does not contain '{value}' entry."
    if "color" not in plot_config[key]:
        assert "parent" in plot_config[key], f"YAML entry {key} is incomplete. Must specify either a 'color' or 'parent' entry."
        if plot_config[key]["parent"] in children_dict:
            children_dict[plot_config[key]["parent"]] += 1
        else:
            children_dict[plot_config[key]["parent"]] = 1
        child_idx = children_dict[plot_config[key]["parent"]]
        parent_color = plot_config[plot_config[key]["parent"]]["color"]
        parent_color = convert_rbg_to_tuple( parent_color )
        plot_config[key]["color"] = convert_tuple_to_rgb( *lighten_color( *parent_color, alpha=1.0-(delta*child_idx) ) )

#----#

with open('lineages.yml', 'r') as f:
        try:
            lineages_yml = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError('Error in lineages.yml file: ' + str(exc))

lineage_info = {}
for lineage in lineages_yml:
    lineage_info[lineage['name']] = {'children': lineage['children']}

# config = checkConfig(config)

for i, i_config in plot_config.items():
    for lin in i_config['members']:
        if '.X' in lin:
            plot_config[i]['members'].extend(lineage_info
                                        [lin.replace('.X', '')]
                                        ['children']
                                        )
# aggregate all the members into list
config_members = [lin for val in plot_config.values() for lin in val[
                  'members']]

# # method for if lineage crumbs present. 
# def group_lineages(lineage_dict):
#     curated_lineages = ['BA.1',
#         'BA.2',
#         'BA.2.12',
#         'BA.4',
#         'BA.5',
#         'XBB',
#         'XBB.1.16',
#         'XBB.1.5',
#         'XBB.1.9',
#         'XBB.1.9.1',
#         'XBB.1.9.2',
#         'EG.5',
#         'BF.7',
#         'BQ.1',
#         'FL.1.5.1']
    
#     lineage_crumbs = lineage_dict['crumbs'].split(';')[::-1]
#     for lin in lineage_crumbs:
#         if lin in curated_lineages:
#             lineage_dict['name'] = lin
#             return lineage_dict
# ​
#     lineage_dict['name'] = 'Other'
#     return lineage_dict



#Load json file
data = []
with open('outputs/aggregate/aggregate_demix.json') as f:
    for line in f:
        data.append(json.loads(line))

df = pd.DataFrame(data)
print(len(df))
#df = df[df['coverage'] > 50]
print(len(df))
df = df.explode('lineages')

df.index = range(len(df))

#############
#map lineages 

def get_name(val, dict0):
    values = []
    for key, value in dict0.items():
        if val in value['members']:
            values.append(dict0[key]['name'])
    return values[0]

import copy
from tqdm import tqdm
config_members = set(config_members)
df['abundance'] = df['lineages'].apply(lambda x:x['abundance'])
df['lineage'] = df['lineages'].apply(lambda x:get_name(x['name'],plot_config) if x['name'] in config_members else "Recombinants" if x['name'][0]=='X' else "Other")
############


# from freyja-sra
# df['lineages'] = df['lineages'].apply(group_lineages)
# df_lineages = pd.json_normalize(df['lineages'])
# df = pd.concat([df.drop(['lineages'], axis=1,), df_lineages], axis=1)
df = df.drop(['lineages'], axis=1,)
###

# Scale by population
df['abundance'] = df['abundance'] * df['ww_population'].astype(float)

# Bin collection dates by week
df['collection_date'] = pd.to_datetime(df['collection_date'], yearfirst=True)
# df['collection_date'] = df['collection_date'].dt.to_period('M').dt.start_time

df = df[['collection_date','lineage','abundance']]

# Aggregate by collection date
df = df.groupby(['collection_date', 'lineage']).agg({'abundance': 'mean'}).reset_index()

# # Make the abundance for each collection date sum to 1
# df['abundance'] = df.groupby('name')['abundance']
df['abundance'] = df['abundance'] / df.groupby('collection_date')['abundance'].transform('sum')
###

df = df.pivot(index='collection_date',columns='lineage').fillna(0)

ordering = [v['name'] for v in plot_config.values()] + ['Recombinants','Other']
colors0 = [v['color'] for v in plot_config.values()] + ["#6A6C6E","#DDDDDD"]
colors0 = colors0[::-1]
df = df.abundance[ordering[::-1]]
### group by week. 
df_ = df.groupby(pd.Grouper(freq='MS')).mean()

fig,ax = plt.subplots(figsize=(10.5,5))

for i in range(0, df_.shape[1]):
    ax.bar(df_.index, df_.iloc[:, i],
           width=20, bottom=df_.iloc[:, 0:i].sum(axis=1),
           label=df_.columns[i], color=colors0[i])
    # ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
ax.set_xlim(df_.index.min()-timedelta(days=15),df_.index.max()+timedelta(days=15))
ax.set_ylim(0,1)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()

plt.savefig('USA_stackplot_Monthly.pdf')



### group by month. 
df_ = df.groupby(pd.Grouper(freq='W')).mean()

fig,ax = plt.subplots(figsize=(10.5,5))

for i in range(0, df_.shape[1]):
    ax.bar(df_.index, df_.iloc[:, i],
           width=7, bottom=df_.iloc[:, 0:i].sum(axis=1),
           label=df_.columns[i], color=colors0[i])
# ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
ax.set_xlim(df_.index.min()-timedelta(days=15),df_.index.max()+timedelta(days=15))
ax.set_ylim(0,1)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()
plt.savefig('USA_stackplot_weekly.pdf')


### average across days, rolling 
df_ = df.groupby(pd.Grouper(freq='D')).mean()
windowSize=28
df_ = df_.rolling(windowSize, center=True,min_periods=0).mean()

fig,ax = plt.subplots(figsize=(10.5,5))
ax.stackplot(df_.index,df_.T,labels=df_.columns,colors=colors0)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlim(df_.index.min(),df_.index.max())
ax.set_ylim(0,1)
ax.set_ylabel('Lineage prevalence')
fig.tight_layout()

plt.savefig('USA_stackplot_daily.pdf')
plt.close('all')
    # print(df)

    # fig = px.bar(df, x='collection_date', y='abundance', color='name', title='Lineage Prevalences in USA Wastewater Samples', barmode='stack')
    # # fig.update_layout(legend=dict(
    # #     orientation="h",
    # #     yanchor="bottom",
    # #     x=0.5,
    # #     y=-1.0,
    # #     xanchor="center"))

    # return fig
