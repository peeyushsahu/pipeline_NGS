## converting mouse gene names into human gene names


library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#genes = c("Zfp286", "Tmx2")


# readd dataframe for genenames
DE_gene_df = read.csv('/ps/imt/e/HL60_Christene/GSE40155_RUNX-KO/results_RunxOE_CT_probeInt_FC2_pval0.05.tsv', sep = '\t')
print(head(DE_gene_df))

reg.gene = DE_gene_df$gene.symbols.x

reg.gene = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = reg.gene ,mart = mouse, attributesL = c("hgnc_symbol","chromosome_name", "start_position"), martL = human, uniqueRows=T)

reg.gene.1 = merge(reg.gene, DE_gene_df, by.x = 'MGI.symbol', by.y = 'gene.symbols.x', all.x = TRUE)
print(head(reg.gene.1))

write.table(x = reg.gene.1, file = '/ps/imt/e/HL60_Christene/GSE40155_RUNX-KO/mouse2human_results_RunxOE_CT_probeInt_FC2_pval0.05.tsv', sep = '\t', row.names = FALSE)
