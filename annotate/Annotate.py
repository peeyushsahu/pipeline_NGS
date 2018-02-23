__author__ = 'peeyush'
import pandas as pd
import alignment.commons as paths
import re
import os
import timeit
import mygene
import numpy as np
import traceback
#Path = paths.path()
#basepath = Path.basepath


def geneid_converter(listofids, input_identifier = None, output_identifier = None):
    '''
    Method will take gene identifier and convert into desired one. Identifiers availble:
    [u'accession', u'alias', u'biocarta', u'chr', u'end', u'ensemblgene', u'ensemblprotein', u'ensembltranscript',
    u'entrezgene', u'exons', u'flybase', u'generif', u'go', u'hgnc', u'homologene', u'hprd', u'humancyc', u'interpro',
    u'ipi', u'kegg', u'mgi', u'mim', u'mirbase', u'mousecyc', u'name', u'netpath', u'pdb', u'pfam', u'pharmgkb', u'pid',
    u'pir', u'prosite', u'ratmap', u'reactome', u'reagent', u'refseq', u'reporter', u'retired', u'rgd', u'smpdb',
    u'start', u'strand', u'summary', u'symbol', u'tair', u'taxid', u'type_of_gene', u'unigene', u'uniprot',
    u'wikipathways', u'wormbase', u'xenbase', u'yeastcyc', u'zfin']
    :param input_identifier: input identifier eg. entrezgene
    :param output_identifier: list of output identifier eg. ["symbol", "ensembl.gene"]
    :param listofids: list of ids to be mapped eg. ['1', '10', '10001']
    :return: DataFrame of mapped ids
    '''
    if input_identifier is None:
        input_identifier = "entrezgene"
    if output_identifier is None:
        output_identifier = ["symbol", "ensembl.gene"]
    mygene_object = mygene.MyGeneInfo()
    mapped_dataframe = mygene_object.querymany(listofids, scopes=input_identifier, fields=output_identifier, species="human", as_dataframe=True)
    return mapped_dataframe


def parse_gtf(path):
    '''
    Reads gtf file and parse it for further use with internal handlers.
    '''
    from pandas import read_csv
    print('Parsing gtf file...')
    colnames = ['chr', 'transcript_annotation', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'further_information']
    gtf_file = read_csv(path, sep='\t')
    gtf_file.columns = colnames
    gtf_file['chr'] = gtf_file['chr'].astype(str)
    return gtf_file


def gtf_binary_search(List, key, imin, imax):
    '''
    Binary search function uses recursive method to search the position/index of peak-summit in the GTF database.
    :param List: DF of sorted GTF file only with 'start'.
    :param key: Position to be searched in the data
    :param imin: First index
    :param imax: Last index
    :return:
    '''
    if key < List[imin]:
        return imin, imin+1
    if key > List[imax]:
        return imax-1, imax
    if imax < imin:
        return 'KEY NOT FOUND'
    elif (imax - imin) == 1 and key > List[imin] and key < List[imax]:
        return imin, imax
    else:
        imid = imin + int((imax - imin)/2)
        if List[imid] > key:
            #print('upper', List[imid], key, imin, imid)
            # key is in upper subset
            return gtf_binary_search(List, key, imin, imid)
        elif List[imid] < key:
            #print('lower', List[imid], key)
            # key is in lower subset
            return gtf_binary_search(List, key, imid, imax)
        elif List[imid] == key:
            #print('found', List[imid], key)
            # key is in lower subset
            return imid, imid+1


class GenerateAnnotationDB:

    def __init__(self, gtf_path, outpath, tss_distance_up=1500, tss_distance_down=500, upstream=5000):
        self.gtf_path = gtf_path
        self.outpath = outpath
        self.tssdistance_up = tss_distance_up
        self.tssdistance_down = tss_distance_down
        self.upstream = upstream
        self.gtf_file = parse_gtf(self.gtf_path)
        self.chr_lengths = self.get_chr_size()

    def call_me(self):
        '''
        Create transcript db and chr vectors
        :return:
        '''
        self.make_gtf_transcript_db()
        self.annotate_chrvectors()

    def get_chr_size(self):
        print('Reading chr sizes...')
        chr_size_path = os.path.join(self.outpath, 'chrom_sizes')
        chr_lengths = pd.read_csv(chr_size_path, header=None, sep='\t')
        chr_lengths[0] = chr_lengths[0].astype(str)
        return chr_lengths

    def make_gtf_transcript_db(self):
        print('Creating transcript db from GTF...')
        gtf_file = self.gtf_file
        gtf_file = gtf_file[gtf_file['feature'] == 'transcript']
        f = open(os.path.join(self.outpath, 'gtf_transcript4vector.db'), mode='wb')
        a = open(os.path.join(self.outpath, 'gtf_transcript4annotation.db'), mode='wb')
        f.write(bytes('chr\tstart\tstop\tstrand\tgene_name\tgene_id\ttranscript_name\ttranscript_id\tgene_biotype\n', 'UTF-8'))
        a.write(bytes('chr\tstart\tstop\tstrand\tgene_name\tgene_id\ttranscript_name\ttranscript_id\tgene_biotype\n', 'UTF-8'))
        for ind, row in gtf_file.iterrows():
            strand = row['strand']
            chr = row['chr']
            start = str(row['start'])
            stop = str(row['stop'])
            gtf_anno = re.split(';| |"', row['further_information'])
            gene_name = gtf_anno[gtf_anno.index('gene_name') + 2]
            tnx_name = gtf_anno[gtf_anno.index('transcript_name') + 2]
            gene_id = gtf_anno[gtf_anno.index('gene_id') + 2]
            tnx_id = gtf_anno[gtf_anno.index('transcript_id') + 2]
            biotype = gtf_anno[gtf_anno.index('gene_biotype') + 2]
            f.write(bytes(chr+'\t'+start+'\t'+stop+'\t'+strand+'\t'+gene_name+'\t'+gene_id+'\t'+tnx_name+'\t'+tnx_id+'\t'+biotype+'\n', 'UTF-8'))
            if strand == '-':
                a.write(bytes(chr+'\t'+stop+'\t'+start+'\t'+strand+'\t'+gene_name+'\t'+gene_id+'\t'+tnx_name+'\t'+tnx_id+'\t'+biotype+'\n', 'UTF-8'))
            else:
                a.write(bytes(chr+'\t'+start+'\t'+stop+'\t'+strand+'\t'+gene_name+'\t'+gene_id+'\t'+tnx_name+'\t'+tnx_id+'\t'+biotype+'\n', 'UTF-8'))
        f.close()

    def get_transcript_db(self):
        tr_db_path = os.path.join(self.outpath, 'gtf_transcript4vector.db')
        transcript_db = pd.read_csv(tr_db_path, header=0, sep='\t')
        transcript_db['chr'] = transcript_db['chr'].astype(str)
        transcript_db = transcript_db[(transcript_db['chr'].str.len() < 4) | (transcript_db['chr'].str.contains('GL'))]
        return transcript_db

    def chromosome_vactors(self):
        '''
        create chromosome zero vector for region annotaion.
        '''
        chr_vector = {}
        #print(self.chr_lengths[0])
        for k, v in self.chr_lengths.iterrows():
            chr_vector[str(v[0])] = np.zeros(v[1]+1, dtype=np.uint32)
        return chr_vector

    def save_chr_vector(self, chr_vector):
        for chro, vec in chr_vector.items():
            if not os.path.exists(os.path.join(self.outpath, 'chr_vector')):
                os.mkdir(os.path.join(self.outpath, 'chr_vector'))
            f = open(os.path.join(self.outpath, 'chr_vector', str(chro)),  "wb")
            np.save(f, vec)
            f.close()

    def annotate_chrvectors(self):
        '''
        fill chromosome vector with respective region annotation
        '''
        print('Annotating chromosome vectors...')
        chr_vector = self.chromosome_vactors()
        # we will annotate tss, exon, upstream
        gene_gtf = self.gtf_file
        gene_gtf = gene_gtf[(gene_gtf['chr'].str.len() < 4) | (gene_gtf['chr'].str.contains('GL'))]

        for feature, sym in zip(['exon', 'transcript'],[1, 3]):
            gtffile = gene_gtf[gene_gtf['feature'] == feature]
            gtffile = gtffile.sort_values(['chr', 'start'], ascending=True)
            #print(gtffile)
            gtf_chr_group = gtffile.groupby('chr')
            # Annotating exons
            if feature == 'exon':
                for chr, df in gtf_chr_group:
                    vector = chr_vector[str(chr)]
                    for ind, row in df.iterrows():
                        start = row['start'] - 1  # converting into o based system
                        stop = row['stop']
                        vector[start: stop] = sym
            # Annotating transcript
            if feature == 'transcript':
                for chr, df in gtf_chr_group:
                    vector = chr_vector[str(chr)]
                    for ind, row in df.iterrows():
                        start = row['start'] - 1
                        stop = row['stop']
                        if row['strand'] == '-':
                            vector[stop-self.tssdistance_down: stop+self.tssdistance_up] = sym
                            tr_vec = np.unique(vector[stop+self.tssdistance_up: stop+self.upstream])
                            if len(tr_vec) == 1 and tr_vec[0] == 0:
                                vector[stop+self.tssdistance_up: stop+self.upstream] = 4
                        else:
                            vector[start-self.tssdistance_up: start+self.tssdistance_down] = sym
                            tr_vec = np.unique(vector[start-self.upstream: start-self.tssdistance_up])
                            if len(tr_vec) == 1 and tr_vec[0] == 0:
                                vector[start-self.upstream: start-self.tssdistance_up] = 4
        # now we will annotate intron
        transcripts = self.get_transcript_db()
        transcript_group = transcripts.groupby('chr')

        for key, df in transcript_group:
            #print(key, list(df['chr'])[0])
            vector = chr_vector[str(key)]
            df = df.sort_values(['chr', 'start'], ascending=[True, True])
            for ind, row in df.iterrows():
                tr_vec = vector[row['start']-1: row['stop']]
                tr_vec[tr_vec == 0] = 2
                tr_vec[tr_vec == 4] = 2
                vector[row['start']-1: row['stop']] = tr_vec
            chr_vector[str(key)] = vector
        self.save_chr_vector(chr_vector) # save chr vector in bin files
        #return chr_vector


class AnnotateNearestGenes:
    '''
    This is will create an object for next gene annotation for peaks.
    '''
    def __init__(self, dataframe, path4annotations='/ps/imt/f/Genomes/geneAnotations', maxdist=1000, maxgenes=2):
        self.dataframe = dataframe
        self.path4annotations = path4annotations
        self.maxdist = maxdist
        self.maxgenes = maxgenes

    def get_gtf4annotation(self):
        tr_db_path = os.path.join(self.path4annotations, 'gtf_transcript4annotation.db')
        transcript_db = pd.read_csv(tr_db_path, header=0, sep='\t')
        transcript_db['chr'] = transcript_db['chr'].astype(str)
        transcript_db = transcript_db[(transcript_db['chr'].str.len() < 4) | (transcript_db['chr'].str.contains('GL'))]
        transcript_db = transcript_db.sort_values(['chr', 'start'], ascending=[True, True])
        transcript_db.index = range(0, len(transcript_db))
        return transcript_db

    def next_genes_annotator(self):
        '''
        This method will take a datatframe and compute 5 nearest genes from the summit of the peak.
        :param dataframe:
        :param path: To GTF file
        :return:
        '''
        print('Process: Reading GTF file')
        start = timeit.default_timer()
        gtffile = self.get_gtf4annotation()
        chr_group = gtffile.groupby('chr')['start']

        dataframe = self.dataframe.sort_values(by='chr', ascending=True)
        nearGene_df = pd.DataFrame()

        print('Process: Annotating peaks with genes')
        for ind, row in dataframe.iterrows():
            try:
                List = chr_group.get_group(row['chr'])
                key = row['start'] + row['summit']
                loc = gtf_binary_search(List, int(key), min(List.index), max(List.index))

                if loc != 'KEY OUT OF BOUND':
                    outdf = self.next_genes(gtffile, loc, row)
                    nearGene_df = nearGene_df.append(outdf, ignore_index=True)
                else:
                    print('Key not found:', key)
            except:
                print('Problem in chr:', row['chr'])
                traceback.print_exc()
                pass

        stop = timeit.default_timer()
        print('\nTime elapsed:', stop-start,' sec')
        return nearGene_df

    def next_genes(self, gtffile, position, row):
        '''
        This method actually search for five nearest genes in GTF file.
        :param gtffile:
        :param position:
        :param maxdist: maximum distance from peak in kbs
        :return:
        '''
        maxdist = self.maxdist * 1000  # converting kb into basepairs
        geneList = []
        minusindex = min(position)
        plusindex = max(position)
        neargene = pd.DataFrame(columns=['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name',
                                         'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'gene_name', 'transcript_name',
                                         'distance2tss', 'strand', 'summit'])
        while len(geneList) < self.maxgenes:
            # print geneList, upstreamindex, downstreamindex
            dist_list = []
            for index in [plusindex, minusindex]:
                if index in gtffile.index:
                    gtf_anno = gtffile.iloc[index]
                    dist2tss = int((row['start'] + row['summit']) - gtf_anno['start'])
                    dist_list.append(abs(dist2tss))
                    gene_name = gtf_anno['gene_name']

                    if (gene_name not in geneList) and (abs(dist2tss) < abs(maxdist)):
                        strand = gtf_anno['strand']
                        tnx_name = gtf_anno['transcript_name']

                        if strand == '-':
                            dist_tss = int((row['start'] + row['summit']) - gtf_anno['start']) * -1
                        else:
                            dist_tss = int((row['start'] + row['summit']) - gtf_anno['start'])

                        newrow = {'chr': row['chr'],
                                'start': int(row['start']),
                                'stop': int(row['stop']),
                                'Next transcript strand': row['Next transcript strand'],
                                'Next transcript gene name': row['Next transcript gene name'],
                                'GenomicPosition TSS=1250 bp, upstream=5000 bp': row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                                'gene_name': gene_name,
                                'transcript_name': tnx_name,
                                'strand': strand,
                                'summit': row['summit'],
                                'distance2tss': dist_tss}

                        neargene = neargene.append(pd.Series(newrow), ignore_index=True)
                        geneList.append(gene_name)
                        if len(geneList) >= self.maxgenes:
                            break
                plusindex += 1
                minusindex -= 1
            try:
                if (dist_list[0] > maxdist) and (dist_list[1] > maxdist):
                    break
            except:
                print(row, dist_list, position)
                traceback.print_exc()
                pass
        return neargene


class AnnotatePeaks:

    def __init__(self, dataframe, path4annotations=None):
        self.dataframe = dataframe
        self.path4annotations = os.path.join(path4annotations, 'geneAnotations')

    def call_me(self):
        '''
        Annotate dataframe with peaks.
        '''
        df = self.annotate_region()
        return self.annotate_distance(df)

    def get_chr_vector(self, chro):
        '''
        Returns chromosome vector for selected chromosome.
        '''
        return np.load(os.path.join(self.path4annotations, 'chr_vector', chro))

    @staticmethod
    def map_annotation(indices):
        '''
        Maps vector annotation with appropriate genomic region string
        :param indices:
        :return:
        '''
        if indices == 1:
            return 'exon'
        if indices == 2:
            return 'intron'
        if indices == 3:
            return 'tss'
        if indices == 4:
            return 'upstream'
        if indices == 0:
            return 'intergenic'

    def annotate_region(self):
        '''
        This will annotate peaks with genomic regions (tss, intron, exon, upstream, intergenic)
        '''
        print('Annotating with genomic regions...')
        peak_df = self.dataframe
        peak_df.index = range(0, len(peak_df))
        peak_df['chr'] = peak_df['chr'].astype(str)
        peak_df_group = peak_df.groupby('chr')

        annotated_df = peak_df#[['chr', 'start', 'stop', 'summit']]
        annotated_df.loc[:, 'genomic_annotation'] = np.array(['NA'] * len(annotated_df))
        ind_ge_annotation = annotated_df.columns.get_loc('genomic_annotation')

        for chr, df in peak_df_group:
            #print('Chr:', chr)
            try:
                vector = self.get_chr_vector(chr)
                for ind, row in df.iterrows():
                    peak_summit = row['start'] + row['summit']
                    indices = vector[peak_summit]
                    #print(ind, ind_ge_annotation)
                    annotated_df.iloc[ind, ind_ge_annotation] = self.map_annotation(indices)

            except:
                print('chromosome not found:', chr)
                traceback.print_exc()
                pass
        return annotated_df

    def get_gtf4annotation(self):
        tr_db_path = os.path.join(self.path4annotations, 'gtf_transcript4annotation.db')
        transcript_db = pd.read_csv(tr_db_path, header=0, sep='\t')
        transcript_db['chr'] = transcript_db['chr'].astype(str)
        transcript_db = transcript_db[(transcript_db['chr'].str.len() < 4) | (transcript_db['chr'].str.contains('GL'))]
        transcript_db = transcript_db.sort_values(['chr', 'start'], ascending=[True, True])
        transcript_db.index = range(0, len(transcript_db))
        return transcript_db

    def annotate_distance(self, annotated_df):
        '''
        This will annotate a peak with nearest transcript and distance to nearest tss.
        :return:
        '''
        print('Annotating with gene and distance...')
        import timeit
        start = timeit.default_timer()

        gtfFile = self.get_gtf4annotation()
        chr_group = gtfFile.groupby('chr')['start']

        dataframe = annotated_df
        dataframe = dataframe.sort_values(['chr'], ascending=True)
        dataframe.loc[:, 'nearestGene'] = np.array(['-'] * len(dataframe))
        dataframe.loc[:, 'nearestGene_id'] = np.array(['-'] * len(dataframe))
        dataframe.loc[:, 'dist2tss'] = np.array([0] * len(dataframe))
        dataframe.loc[:, 'strand'] = np.array(['-'] * len(dataframe))
        dataframe.loc[:, 'nearestTranscript'] = np.array(['-'] * len(dataframe))
        dataframe.loc[:, 'nearestTranscript_id'] = np.array(['-'] * len(dataframe))
        dataframe.loc[:, 'gene_biotype'] = np.array(['-'] * len(dataframe))
        ind_nearestTss = dataframe.columns.get_loc('nearestGene')
        ind_nearestTss_id = dataframe.columns.get_loc('nearestGene_id')
        ind_dist2tss = dataframe.columns.get_loc('dist2tss')
        ind_strand = dataframe.columns.get_loc('strand')
        ind_tnx = dataframe.columns.get_loc('nearestTranscript')
        ind_tnx_id = dataframe.columns.get_loc('nearestTranscript_id')
        ind_biotype = dataframe.columns.get_loc('gene_biotype')

        for ind, row in dataframe.iterrows():
            try:
                List = chr_group.get_group(row['chr'])
                key = row['start'] + row['summit']
                #print(List.shape, key, min(List.index), max(List.index))
                loc = gtf_binary_search(List, int(key), min(List.index), max(List.index))

                if not loc == 'KEY NOT FOUND':
                    nearest_gene = self.next_genes(gtfFile, loc, row)
                    #print(nearest_gene)
                    dataframe.iloc[ind, ind_nearestTss] = list(nearest_gene.keys())[0]
                    dataframe.iloc[ind, ind_dist2tss] = list(nearest_gene.values())[0][0]
                    dataframe.iloc[ind, ind_strand] = list(nearest_gene.values())[0][1]
                    dataframe.iloc[ind, ind_tnx] = list(nearest_gene.values())[0][2]
                    dataframe.iloc[ind, ind_tnx_id] = list(nearest_gene.values())[0][3]
                    dataframe.iloc[ind, ind_biotype] = list(nearest_gene.values())[0][4]
                    dataframe.iloc[ind, ind_nearestTss_id] = list(nearest_gene.values())[0][5]
                else:
                    print('Key not found:', key)
            except:
                print('Problem in chr:', row['chr'])
                traceback.print_exc()
                pass
        stop = timeit.default_timer()
        print('\nTime elapsed:', stop-start,' sec')
        return dataframe

    @staticmethod
    def next_genes(gtffile, position, row):
        '''
        This method actually search for nearest genes in GTF file.
        :param gtffile:
        :param position:
        :param maxdist: maximum distance from peak in kbs
        :return:
        '''
        nearest_gene = {}
        plusindex = position[0]
        minusindex = position[1]
        dist = float('inf')
        strand = '-'
        gene_name = '-'
        tnx_name = '-'
        tnx_id = '-'
        biotype = '-'
        gene_id = '-'
        dist_tss = 0
        # print geneList, upstreamindex, downstreamindex
        for index in [plusindex, plusindex+1, minusindex, minusindex-1]:
            if index in gtffile.index:
                gtf_anno = gtffile.iloc[index]
                dist2tss = int((row['start'] + row['summit']) - gtf_anno['start'])
                if abs(dist2tss) < abs(dist):
                    strand = gtf_anno['strand']
                    gene_name = gtf_anno['gene_name']
                    tnx_name = gtf_anno['transcript_name']
                    tnx_id = gtf_anno['transcript_id']
                    biotype = gtf_anno['gene_biotype']
                    gene_id = gtf_anno['gene_id']

                    if strand == '-':
                        dist_tss = int((row['start'] + row['summit']) - gtf_anno['start']) * -1
                    else:
                        dist_tss = int((row['start'] + row['summit']) - gtf_anno['start'])

                    dist = dist2tss
        nearest_gene[gene_name] = [dist_tss, strand, tnx_name, tnx_id, biotype, gene_id]
        return nearest_gene












