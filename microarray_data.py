__author__ = 'peeyush'

from pandas import read_csv
import pandas as pd


df1 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/HW_diff_exp/Exp1.csv', sep=',')
df2 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/HW_diff_exp/Exp2.csv', sep=',')
df3 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/HW_diff_exp/Exp3.csv', sep=',')
mean_df = pd.DataFrame(df1)
for row in range(0, df1.shape[0]):
    for col in range(1, df1.shape[1]):
        mean = df1.iloc[row, col] + df2.iloc[row, col] + df3.iloc[row, col]
        #print mean
        if not isinstance(mean, basestring):
            mean_df.iloc[row, col] = mean/3
            #print 'div', mean_df.iloc[row, col]
        else:
            mean_df.iloc[row, col] = mean

for row in range(0, mean_df.shape[0]):
    for col in range(1, mean_df.shape[1]):
        mean = mean_df.iloc[row, col]
        #print mean
        if isinstance(mean, basestring):
            if mean.count('NC') > 1:
                mean_df.iloc[row, col+1] = 0
                #print 'div', mean_df.iloc[row, col]

mean_df.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/HW_diff_exp/mean_FC_HW.csv', sep='\t', header=True, index=False)


expressiondf = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/until72_HW.csv', sep=',')
peakdf = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/nearest_gene_diffPeaks_P6.csv', sep=',')

def expression4peaks(peakdf, expressiondf):
    '''
    This function will take a dataframe consisting genenames and the expression dataframe.
    :return:
    '''
    expressiondf = expressiondf.set_index(expressiondf['SYMBOL'].str.lower())
    newDF = pd.DataFrame(peakdf[['chr', 'start', 'stop', 'Next transcript gene name', 'log2FoldChange_P6', 'next5genes']])
    newDF['next5genes_expr'] = 0
    newDF['expression'] = 0
    for index, cols in newDF.iterrows():
        neargene = cols['next5genes'].split(',')
        for gene in neargene:
            print gene
            if gene.lower() in expressiondf.index:
                expr = expressiondf.loc[gene.lower(), 'FC_72'].tolist()
                if type(expr) == list:
                    expr = sum(expr)/2
                if newDF.loc[index, 'next5genes_expr'] == 0:
                    newDF.loc[index, 'next5genes_expr'] = gene+':'+str(expr)
                    newDF.loc[index, 'expression'] = expr
                else:
                    newDF.loc[index, 'next5genes_expr'] = newDF.loc[index, 'next5genes_expr'] + ';' + gene + ':' + str(expr)
                    #newDF.loc[index, 'expression'] = str(newDF.loc[index, 'expression']) + ';' + str(expr)
    newDF.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/peaks2expression.csv', sep=',', header=True, index_label=False, index=False)
    return newDF
