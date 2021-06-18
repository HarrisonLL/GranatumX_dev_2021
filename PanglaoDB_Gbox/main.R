source('./granatum_sdk.R')

# get arguments
# arguments can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleArgument"
#num_clusters <- gn_get_arg("exampleArgument")

# get imports
# imports can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleImport"
assay_file <- gn_get_uploaded_file_path("assayFile")
tmp <- load(assay_file)
genemat <- get(tmp)

# our gn_get_import() returned an object of the "assay" class
# this is a named list with 3 keywords: "matrix" contains the
# gene expression values, "sampleIds" contain the cell IDs
# and "geneIds" contain the gene IDs

# INSERT THE BODY OF YOUR SCRIPT HERE
# for example, you may transpose the gene expression data as
# assay$matrix <- assay$matrix.T
#print(genemat@Dim[1])
#print(genemat@Dim[2])
datamat <- rep(0,genemat@Dim[1]*genemat@Dim[2])
#dgcmat <- as.matrix(summary(genemat))

#print(ncol(dgcmat))
#print(nrow(dgcmat))

for(j in 1:length(genemat@i)){
	datamat[genemat@i[j] + 1] = genemat@x[j]
}

Pangassay <- list(matrix = datamat,
		  sampleIds = dimnames(genemat)[2],
		  geneIds = dimnames(genemat)[1])
# export results
# here, we use keywords specified in the "extractFrom"
# field in the exports section of the package.yaml file - in 
# our case, "exampleExport"
gn_export_statically(Pangassay,"PanglaoDB assay")
gn_add_result("Successfully download PanglaoDB data!",data_type="markdown")
gn_commit()
