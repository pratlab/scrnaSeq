#' Preprocess 10X
#'
#' 
#' This function takes a data-frame (genes x cells), 
#' creates a Seurat object with it and filters the object for default tags such as  
#' Min and max nFeature_RNA and % of MT
#' 
#'
#' @param data data-frame
#' @param percent_mt integer [0-100]
#' @param max_features integer [0-Inf]
#' @param min_features integer [0-Inf]
#' @return Preprocessed Seurat object
#' @export
preprocess_10X<-function(data,percent_mt=20, max_features=5000, min_features=200){

	seurat.obj = CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
	seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = '^mt-|^MT-')
	seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt<percent_mt)
	return(seurat.obj)
}
#' Load 10X
#'
#' 
#' This function takes the path to a 10X output folder and instanciates the Seurat object
#'
#' @param file string (path to file)
#' @param percent_mt integer [0-100]
#' @param max_features integer [0-Inf]
#' @param min_features integer [0-Inf]
#' @return Preprocessed Seurat object
#' @export
load_10X<-function(file,name='10X-project',percent_mt=20, max_features=5000, min_features=200){
	data = Read10X(data.dir = file)
	seurat.obj = preprocess_10X(data,percent_mt=percent_mt,max_features=max_features,min_features=min_features)
	return(seurat.obj)
}
#' Reduce Dimensions
#'
#' 
#' This function reduces the dimensions of the Normalized Seurat object
#' Runs PCA and then UMAP and then performs clustering 
#'
#' @param seurat.obj S4 instance
#' @param ndims integer [3-100]
#' @param res double [0-3]
#' @return Seurat object with reduction embeddings
#' @export
reduce_dim<-function(seurat.obj,ndims=15,res=.1){
	seurat.obj <- RunPCA(seurat.obj, verbose = TRUE)
	seurat.obj <- RunUMAP(seurat.obj, dims = 1:ndims)
	seurat.obj <- FindNeighbors(seurat.obj, dims = 1:ndims)
	seurat.obj <- FindClusters(seurat.obj, resolution = res)
	return(seurat.obj)	
}
#' Integrate batches
#'
#' 
#' Use the SCT approach of Seurat to integrate batches together 
#' 
#'
#' @param obj.list list of the S$ objects to be merged
#' @return Integrated Seurat object
#' @export
integrate_batches<-function(obj.list){

	features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 5000)
	obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = TRUE)
	anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
	    anchor.features = features, verbose = TRUE)
	integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
	    verbose = TRUE)
	integrated<-reduce_dim(integrated)
	return(integrated)
	
}
#' Merge samples 
#'
#' 
#' This function allows to instantiate a merged Seurat object of several distinct samples
#' 
#'
#' @param files character vector of the files to merge 
#' @return Merged Seurat object 
#' @export
merge_samples<-function(files){

	to.merge=list()
	for(f in names(files)){
		to.merge[[f]]=load_10X(files[[f]])
	}	
	to.merge = merge(to.merge[[1]],to.merge[2:length(to.merge)])
	to.merge = SCTransform(to.merge)
	to.merge = reduce_dim(to.merge)
	return(to.merge)
}
#' Recluster
#'
#' 
#' Takes a subset of a Seurat object and performs reclustering by instanciating a new object  
#' 
#'
#' @param sub.obj Seurat object 
#' @return Processed reclustering
#' @export
recluster<-function(sub.obj,name='10X-subcluster',batch=NA){

	if(is.na(batch)){
		new = preprocess_10X(sub.obj[['RNA']]@data)
		new = SCTransform(new)
		new = reduce_dim(new)
	}else{
		obj.list=list()
		for(name in names(sub.obj)){
			obj.list[[name]]=preprocess_10X(sub.obj[['RNA']]@data)
			obj.list[[name]]=SCTransform(obj.list[[name]])
		}
		new = integrate_batches(obj.list)
	}
	return(new)
	
}