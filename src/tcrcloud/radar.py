from sklearn import preprocessing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skbio
import airr

plt.rcParams["font.family"] = "serif"

def convergence(df,length):
    counts=df['counts']
    total=counts.sum()
    conv_value=0
    for i in df['counts']:
        if i > 1:
            conv_value += (i/total)
    return conv_value

def Dfifty(df,length):
    counts=df['counts']
    total=counts.sum()
    counter=0
    
    for i in counts.cumsum():
        counter += 1
        if i > total/2:
            break
    if length >= 10000:
        return (counter*100)/10000
    elif length < 10000:
        return (counter*100)/length

def format_data(args):
    validate = airr.validate_rearrangement(args.rearrangements) 
    reader = airr.read_rearrangement(args.rearrangements)
    empty_list=[]

    # keep only the Junction, Vgene and Repertoire ID columns
    keys = ['junction_aa', 'v_call', 'junction','repertoire_id']

    for row in reader:
        empty_list.append({x:row[x] for x in keys})

    df = pd.DataFrame(empty_list)

    # replace cells without junction with Nan
    df['junction_aa'].replace('', np.nan, inplace=True)
    
    # delete lines with Nan
    df.dropna(subset=['junction_aa'], inplace=True)
    
    # delete lines with an X on the junction
    df = df[~df.junction_aa.str.contains("X")]
    
    # Series to select only one Vgene when there are multiple in the column
    temp_v_call  = df.v_call.str.split(",", n = 0, expand = True)
    
    # replace the v_call column with the new one
    df['v_call'] = temp_v_call[0]

    # format the df to aggregate by junction and count them
    aggregate=df.pivot_table(index=['junction_aa','junction','v_call'], 
                                aggfunc='size').reset_index()

    # replace the old information in the df
    df['junction_aa']= aggregate['junction_aa']
    df['junction']= aggregate['junction']
    df['v_call']= aggregate['v_call']
    df['counts']= aggregate[0]
    
    # the size of the df changed so removing extra lines
    df.dropna(subset=['junction_aa'], inplace=True)
    
    # sort the df by highest numbers of counts
    df=df.sort_values(by='counts', ascending=False)
    
    # remove allele information from v_call and keep only the gene information
    df['v_call'] = df.apply(lambda x: x['v_call'][:-3], axis = 1)
    df['chain'] = df.apply(lambda x: x['v_call'][2], axis = 1)

    return df

def radar(args):
    categories = [
                    'D50\nIndex',
                    'Chao1\nIndex',
                    'Convergence', 
                    'Shannon\nIndex', 
                    'Simpson\nIndex', 
                    'Gini\nindex',
                    'Distinct\nCDR3',]
    categories = [*categories, categories[0]]

    minmax_scale = preprocessing.MinMaxScaler(feature_range=(0, 1))

    samples_df=format_data(args)
    samples = samples_df.groupby(['chain','repertoire_id'])
    keys = [key for key, _ in samples]
    patients=[]

    for j in keys:
        df=samples.get_group(j)
        length=len(df.pivot_table(index=['junction'], 
                                aggfunc='size').reset_index())
        counts = df['counts'].tolist()
        distinct= np.array([0,length,10000]).reshape(-1, 1)
        convergencelist= np.array([0,convergence(df,length),1]).reshape(-1, 1) 
        shannon= np.array([0,skbio.diversity.alpha_diversity('shannon',counts)[0],15]).reshape(-1, 1)
        simpson= np.array([0,skbio.diversity.alpha_diversity('simpson',counts)[0],1]).reshape(-1, 1)
        chao= np.array([0,skbio.diversity.alpha_diversity('chao1',counts)[0],120000]).reshape(-1, 1)
        gini= np.array([0,skbio.diversity.alpha_diversity('gini_index',counts)[0],1]).reshape(-1, 1)
        dfiftylist= np.array([0,Dfifty(df,length),1]).reshape(-1, 1)
        metrics=[]
        metrics.append(j[1]+' '+j[0])
        metrics.append(minmax_scale.fit_transform(dfiftylist)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(chao)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(convergencelist)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(shannon)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(simpson)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(gini)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(distinct)[1].astype(np.float))
        patients.append(metrics)
    
    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(patients[0]))

    plt.figure(figsize=(10,12))
    plt.subplot(polar=True)
    len_label=0
    radar_colours=['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                    '#882255', '#44AA99', '#999933', '#AA4499']
    for i in patients:        
        patient = i[1:]
        patient = [*patient, patient[0]]
        try:
            thecolour=radar_colours.pop()
        except IndexError:
            thecolour='#BBBBBB'
        
        plt.plot(label_loc, patient, label=i[0], linewidth=4.0,color=thecolour)
        plt.fill(label_loc, patient, linewidth=4.0, alpha=0.6,color=thecolour)
        if len(i[0]) > len_label:
            len_label = len(i[0])
    plt.title('Repertoire profile\ncomparison', size=20, y=1.10)

    plt.text(label_loc[0], 0.06,'5',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[0], 0.26,'15',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[0], 0.46,'25',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[0], 0.66,'35',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[0], 0.86,'45',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[0], 1.03,'50',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    
    plt.text(label_loc[1], 0.06,'11000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[1], 0.26,'33000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[1], 0.46,'55000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[1], 0.66,'77000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[1], 0.86,'99000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[1], 1.03,'110000',horizontalalignment='left',verticalalignment='center',fontsize=12,fontweight='bold')
    
    plt.text(label_loc[2], 0.06,'0.1',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[2], 0.26,'0.3',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[2], 0.46,'0.5',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[2], 0.66,'0.7',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[2], 0.86,'0.9',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[2], 1.03,'1',horizontalalignment='center',verticalalignment='bottom',fontsize=12,fontweight='bold')
    
    plt.text(label_loc[3], 0.06,'1.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[3], 0.26,'4.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[3], 0.46,'7.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[3], 0.66,'10.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[3], 0.86,'13.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[3], 1.03,'15',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    
    plt.text(label_loc[4], 0.06,'0.1',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[4], 0.26,'0.3',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[4], 0.46,'0.5',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[4], 0.66,'0.7',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[4], 0.86,'0.9',horizontalalignment='right',verticalalignment='center',fontsize=12,fontweight='bold')
    plt.text(label_loc[4], 1.03,'1',horizontalalignment='right',verticalalignment='top',fontsize=12,fontweight='bold')
    
    plt.text(label_loc[5], 0.06,'0.1',horizontalalignment='center',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[5], 0.26,'0.3',horizontalalignment='center',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[5], 0.46,'0.5',horizontalalignment='center',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[5], 0.66,'0.7',horizontalalignment='center',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[5], 0.86,'0.9',horizontalalignment='center',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[5], 1.03,'1',horizontalalignment='center',verticalalignment='top', fontsize=12,fontweight='bold')
    
    plt.text(label_loc[6], 0.06,'1000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[6], 0.26,'3000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[6], 0.46,'5000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[6], 0.66,'7000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[6], 0.86,'9000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    plt.text(label_loc[6], 1.03,'10000',horizontalalignment='left',verticalalignment='center', fontsize=12,fontweight='bold')
    
    plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9],[])
    plt.tick_params(pad=30,labelsize=16)
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    outputname = args.rearrangements[:-4]+"_radar"+".png"
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2),fontsize=16)
    plt.tight_layout()
    plt.savefig(outputname,dpi=300)