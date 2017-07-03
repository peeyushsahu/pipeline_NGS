__author__ = 'peeyush'
import math


def peakLengthDistribution(peak_data):
    '''
    Compute distribution of peak length.
    :param peak_data:
    :return:
    '''
    length = []
    hist_list = [0]*30
    bins = [0]*30
    nOfBins = 30
    for count, row in peak_data.iterrows():
        length.append(row['stop']-row['start'])
    maxL = max(length)
    binSpan = math.floor(maxL/30)
    print ("Length of one bin is:", binSpan)
    for i in length:
        binNo = math.floor(i/binSpan)
        #print binNo
        if binNo <= 0:
            hist_list[0] += 1
        else:
            hist_list[int(binNo)-1] += 1
    count = 0
    for i in range(0, 30):
        bins[i] = int(count)
        count = count + binSpan
    return length, hist_list, bins

#plt.hist(length, bins=30, histtype='stepfilled', normed=True, color='b', label='H3K4')
#plt.hist(lengt1, bins=30, histtype='stepfilled', normed=True, color='r', alpha=0.5, label='PRMT6')
#plt.title("H3K4/PRMT6 peak length distribution")
#plt.xlabel("Peak Length")
#plt.ylabel("Occurrence")
#plt.legend()
#plt.show()
#plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/plots/' + self.name + '.pdf')
#plt.clf()