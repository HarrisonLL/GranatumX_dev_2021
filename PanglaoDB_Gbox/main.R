source('./granatum_sdk.R')

# get arguments
# arguments can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleArgument"
num_clusters <- gn_get_arg("exampleArgument")

# get imports
# imports can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleImport"
assay <- gn_get_import("exampleImport")

# our gn_get_import() returned an object of the "assay" class
# this is a named list with 3 keywords: "matrix" contains the
# gene expression values, "sampleIds" contain the cell IDs
# and "geneIds" contain the gene IDs

# INSERT THE BODY OF YOUR SCRIPT HERE
# for example, you may transpose the gene expression data as
assay$matrix <- assay$matrix.T

# export results
# here, we use keywords specified in the "extractFrom"
# field in the exports section of the package.yaml file - in 
# our case, "exampleExport"
gn_export_statically(assay,"exampleExport")
gn_add_result("Example concluding message",data_type="markdown")
gn_commit()
