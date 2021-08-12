source('./granatum_sdk.R')
library(RJSONIO)
library(jsonlite)
library(reticulate)
reticulate::install_miniconda()

# get arguments
# arguments can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleArgument"
#num_clusters <- gn_get_arg("exampleArgument")


urlsafebase64encode <- function(x, ...){
	gsub("+", "-", gsub("/", "_", base64_enc(x), fixed = TRUE), fixed = TRUE)
}
# get imports
# imports can be accessed using keywords specified 
# in the "injectInto" field of the arguments section
# of the package.yaml file - in our case, "exampleImport"
assay_file <- gn_get_uploaded_file_path("assayFile")

#sra_num <- gn_get_arg("SRA_num")
#srs_num <- gn_get_arg("SRS_num")
#url_head <- "https://panglaodb.se/data_dl.php?"
#url_tail <- "filetype=R&datatype=readcounts"
#url_head <- paste(url_head, sra_num, sep = "sra=")
#url_head <- paste(url_head, srs_num, sep = "&srs=")
#url_link = paste(url_head, url_tail, sep = "&")

#file_path = "./"
#file_name = paste(sra_num, srs_num, sep = "-")
#file_name = paste(file_name, "Rdata", sep = ".")
#dl_target = paste(file_path, file_name, sep = "")
#download.file(url_link, dl_target, mode = "wb")
#system("wget https://panglaodb.se/data_dl.php?sra=SRA553822&srs=SRS2119548&filetype=R&datatype=readcounts")
#dir()

#tmp <- load(file_name)

#tmp <- load("SRA553822_SRS2119548.sparse.RData")
tmp <- load(assay_file)
genemat <- get(tmp)
timestart<-Sys.time()
# our gn_get_import() returned an object of the "assay" class
# this is a named list with 3 keywords: "matrix" contains the
# gene expression values, "sampleIds" contain the cell IDs
# and "geneIds" contain the gene IDs

# INSERT THE BODY OF YOUR SCRIPT HERE
# for example, you may transpose the gene expression data as
# assay$matrix <- assay$matrix.T
#print(genemat@Dim[1])
#print(genemat@Dim[2])

start <- 1
step <- 1000
gene_num <- genemat@Dim[1]
cell_num <- genemat@Dim[2]

print(gene_num, flush=TRUE)
print(cell_num, flush=TRUE)

#dgcmat <- as.matrix(summary(genemat))

#print(ncol(dgcmat))
#print(nrow(dgcmat))
row_pos <- genemat@i+1
col_pos <- findInterval(seq(genemat@x)-1,genemat@p[-1])+1
val <- genemat@x
#col_num <- 1
#prev_idx <- 0
#index_ptr <- 1
output = list()
output[["origin data size"]] = c(gene_num, cell_num)
output[["current chunk"]] = c("test", "col")
output[["suggested chunk"]]= list()
output[["suggested chunk"]][["test"]] = c("col",1000)
output[["suggested chunk"]][["deepimpute"]] = c("col", 3456)
output[["suggested chunk"]][["log-transform"]] = c("row", 3456)
#c("test"=c("col",1000), "deepimpute"=c("col", 3456), "log-transform"=c("row", 3456))

count = 1

while (start <= cell_num){

	end <- start + step - 1
	if (end > cell_num){
		end <- cell_num
	}
	datamat <- rep(0,gene_num * (end - start + 1))
	datamat <- matrix(datamat, ncol = (end - start + 1)) 
	for (i in seq_along(val)){
		if (col_pos[i] <= end && col_pos[i] >= start){
			if(col_pos[i] %% step == 0){
				datamat[row_pos[i],step] <- val[i]
			}
			else{
				datamat[row_pos[i],col_pos[i]%%step] <- val[i]
			}
		}
	}
	exported_assay <- list(matrix = datamat,
			       sampleIds = genemat@Dimnames[[2]][start:end],
			       geneIds = genemat@Dimnames[[1]])
	#compressed_assay <- memCompress(toJSON(exported_assay), type="gzip")
	#print("encoding..", flush=TRUE)
	#base64_assay <- urlsafebase64encode(compressed_assay)
	#print(base64_assay[1:10], flush = TRUE)
	
	# reticulate sourcing python script
	reticulate::source_python("compression.py")
	base64_assay <- compress_chunk(exported_assay)
	
	assay_name <- paste("chunk", count, sep = "")
	output[[assay_name]] <- base64_assay
	print("finishing one chunk", flush=TRUE)
	count <- count + 1
	rm('exported_assay')
	rm('datamat')
	rm('compressed_assay')
	rm('base64_assay')
	start <- end + 1
	print(Sys.time()-timestart, flush = TRUE)
}

#print("dimension transform")
#datamat = matrix(datamat, ncol = cell_nums)
#dim(datamat) <- c(genemat@Dim[1],cell_nums)
#print(ncol(datamat))
#print(nrow(datamat))

###############################################
#Random select columns
#export_data = as.data.frame(datamat)
#random_cols = sample(ncol(export_data), 3000)
#output = as.matrix(export_data[,random_cols])
###############################################
#Pangassay <- list(matrix = output,
#		  sampleIds = genemat@Dimnames[[2]][random_cols],
#		  geneIds = genemat@Dimnames[[1]])
# export results
# here, we use keywords specified in the "extractFrom"
# field in the exports section of the package.yaml file - in 
# our case, "exampleExport"
print("start to export..", flush = TRUE)

gn_export_statically(output,"PanglaoDB assay")
gn_add_result("Successfully download PanglaoDB data!",data_type="markdown")
gn_commit()
