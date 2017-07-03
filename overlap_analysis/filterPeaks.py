__author__ = 'peeyush'
import alignment.commons as paths
import os
Path = paths.path()
basepath = Path.basepath


def filterpeaks(peak_data, name, filtering=True):
    #filtered_peak_data = {}
    #for k, v in peak_data.iteritems():
    #print k, v.shape
    name = name
    df = peak_data
    sample = name.split(' ')
    colnames = df.columns.values.tolist()
    if filtering:
        print(colnames)
        indices1 = [i for i, s in enumerate(colnames) if sample[0] in s]
        indices2 = [i for i, s in enumerate(colnames) if sample[2] in s]
        #indices3 = [i for i, s in enumerate(colnames) if "Input" in s]
        for i in indices1:
            if "RA" not in colnames[i] and "norm" not in colnames[i]:
                condition = colnames[i]
                break
            elif "RA" in colnames[i] and "norm" not in colnames[i]:
                condition = colnames[i]
                break
        for i in indices2:
            if "RA" not in colnames[i] and "norm" not in colnames[i]:
                control = colnames[i]
                break
            elif "RA" in colnames[i] and "norm" not in colnames[i]:
                control = colnames[i]
                break
        else:
            raise ValueError("Filtering Sample name differs from column name.")
        print('Sample lane:'+condition)
        print('Control lane:'+control)
        #inputcol = colnames[indices3[0]]
        #print inputcol
        exclude_from_filtering = ['H3K36me3', 'H3K27me3', 'H3K4me1']

        ## condition for simple filtering of preaks
        if any(s in condition for s in exclude_from_filtering):
            df1 = df[df[condition] >= 2*df[control]]
            final = df1
        else:
            df1 = df[df[condition] >= 2*df[control]]
            df2 = df1[((df1['stop']-df1['start'])/df1[condition]) <= 15]
            final = df2
        print(final.shape)
    else:
        final = df
        print(final.shape)
    with open(basepath + '/further_analysis/filtered/filteredPeaksCount.txt', 'a') as file:
        file.write(name+'\t'+str(len(df))+'\t'+str(len(final))+'\n')
    #filtered_peak_data[name] = final
    dir_path = basepath + '/further_analysis/filtered/'+name
    paths.ensure_path(dir_path)
    samPath = os.path.join(dir_path, name+'.txt')
    final.to_csv(samPath, sep="\t", header=True)
    final.index = range(len(final))
    return final, dir_path


def peaks_in_allsamples(peaksdirpath=None):
    '''
    Create a txt file for peaks count in all the files in csv folder
    :param peaksdirpath:
    :return:
    '''
    import pandas as pd
    if peaksdirpath is None:
        peaksdirpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/'
    #print(os.listdir(peaksdirpath))
    files = os.listdir(peaksdirpath)
    file = open(basepath + '/further_analysis/RawPeakCount.txt', 'w')
    file.write('Sample_name\tcontrol_sample\traw_peaks\n')
    if len(files) > 0:
        for peakfile in files:
            sample = peakfile.split('vs')
            if (os.path.isfile(os.path.join(peaksdirpath, peakfile))) and (peakfile.endswith('.csv')):
                #print('here')
                df = pd.read_csv(os.path.join(peaksdirpath, peakfile), sep='\t', header=0)
                if len(sample) > 1:
                    file.write(sample[0].strip()+'\t'+sample[1].strip()+'\t'+str(len(df))+'\n')
                else:
                    file.write(sample[0].strip()+'\t--\t'+str(len(df))+'\n')
    print('Done!!')
    file.close()



def extract_data(dataset, list_of_values, column):
    '''
    This will filter data based on values in column provided
    :param dataset: dataframe
    :param idlist: list
    :return: dataframe
    '''
    newdataset = dataset[dataset[column].isin(list_of_values)]
    return newdataset