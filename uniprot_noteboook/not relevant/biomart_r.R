require(biomaRt)



mart<-useMart(biomart = “ensembl”, dataset =  “mmusculus_gene_ensembl”)mart <- useDataset(dataset=”mmusculus_gene_ensembl”, mart=mart)mapping <-  getBM(attributes=c(“mgi_symbol”,”ensembl_gene_id”,”entrezgene_id”), filters = “mgi_symbol”, mart=mart, values=data, uniqueRows=TRUE, bmHeader = T)