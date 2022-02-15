from sklearn import preprocessing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skbio
import airr

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
                    'Distinct\nCDR3',
                    'Convergence', 
                    'Shannon\nIndex', 
                    'Simpson\nIndex', 
                    'Chao1\nIndex',
                    'Gini\nindex',
                    'D50\nIndex']
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
        distinct= np.array([0,length,50000]).reshape(-1, 1)
        convergencelist= np.array([0,convergence(df,length),1]).reshape(-1, 1) 
        shannon= np.array([0,skbio.diversity.alpha_diversity('shannon',counts)[0],15]).reshape(-1, 1)
        simpson= np.array([0,skbio.diversity.alpha_diversity('simpson',counts)[0],1]).reshape(-1, 1)
        chao= np.array([0,skbio.diversity.alpha_diversity('chao1',counts)[0],70000]).reshape(-1, 1)
        gini= np.array([0,skbio.diversity.alpha_diversity('gini_index',counts)[0],1]).reshape(-1, 1)
        dfiftylist= np.array([0,Dfifty(df,length),1]).reshape(-1, 1)
        metrics=[]
        metrics.append(j[1]+' '+j[0])
        metrics.append(minmax_scale.fit_transform(distinct)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(convergencelist)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(shannon)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(simpson)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(chao)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(gini)[1].astype(np.float))
        metrics.append(minmax_scale.fit_transform(dfiftylist)[1].astype(np.float))
        patients.append(metrics)
    
    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(patients[0]))

    plt.figure(figsize=(16,8))
    plt.subplot(polar=True)
    for i in patients:        
        patient_1 = i[1:]
        patient_1 = [*patient_1, patient_1[0]]
        plt.plot(label_loc, patient_1, label=i[0], linewidth=4.0)
        plt.fill(label_loc, patient_1, label=i[0], linewidth=4.0, alpha=0.6)
    plt.title('Repertoire profile\ncomparison', size=20, y=1.05)
    plt.yticks([])
    plt.tick_params(pad=22,labelsize=16)
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    
    outputname = args.rearrangements[:-4]+"_radar"+".png"
    plt.legend(loc='upper right',bbox_to_anchor=(1.7, 1.0))
    plt.savefig(outputname,dpi=300)