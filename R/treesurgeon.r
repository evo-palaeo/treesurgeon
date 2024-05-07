###################### treesurgeon v.0.2 ###################### 
################### J. N. Keating, 06/04/24 ###################

################### Data ###################

#' Early vertebrate skeleton data
#'
#' Phylogeny and associated morphological character data modified from Keating & Donoghue (2016), O'shea et al (2019), de Buffrénil et al (2016) and Witzmann (2009). 
#'
#' @docType data
#'
#' @usage data(vert_data)
#'
#' @format time tree of class 'phylo' representing the relationships of the 71 vertebrate taxa. Associated morphological character data ('morph') is included as an object of class 'data.frame'.
#'
#' @keywords datasets
#'
#' @references Keating, J.N. and Donoghue, P.C., 2016. Histology and affinity of anaspids, and the early evolution of the vertebrate dermal skeleton. Proceedings of the Royal Society B: Biological Sciences, 283(1826), p.20152917. \cr 
#' O'Shea, J., Keating, J.N. and Donoghue, P.C., 2019. The dermal skeleton of the jawless vertebrate Tremataspis mammillata (Osteostraci, stem‐Gnathostomata). Journal of morphology, 280(7), pp.999-1025. \cr 
#' de Buffrénil, V., Clarac, F., Canoville, A. and Laurin, M., 2016. Comparative data on the differentiation and growth of bone ornamentation in gnathostomes (Chordata: Vertebrata). Journal of Morphology, 277(5), pp.634-670. \cr 
#' Witzmann, F, 2009. Comparative histology of sculptured dermal bones in basal tetrapods, and the implications for the soft tissue dermis. Palaeodiversity, 2(233), p.e270.
#' 
#' @examples
#' data(vert_data)
#' plot(vert_data, cex = 0.5, no.margin = T)
#' head(vert_data$morph)
"vert_data"


#' Dinosaur-bird integument data
#'
#' Phylogeny and associated integument character data modified from Cockx et al (2024, in prep). 
#'
#' @docType data
#'
#' @usage data(dino_data)
#'
#' @format time tree of class 'phylo' representing the relationships of the 96 archosaur taxa (predominantly theropods). Integument character data ('morph') is included as an object of class 'data.frame'. Character state decriptions ('state_labels') for each of the three integument characters are included as a list.
#'
#' @keywords datasets
#'
#' @references Cockx et al (2024, in prep)
#' 
#' @examples
#' data(dino_data)
#' plot(dino_data, cex = 0.5, no.margin = T)
#' head(dino_data$morph)
"dino_data"


#' Total-evidence dataset from Lavoue (2016).
#'
#' A list comprising morphological and molecular data partitions which make up the total-evidence phylogenetic dataset of Osteoglossiformes from Lavoue (2016). 
#'
#' @docType data
#'
#' @usage data(Lavoue2016)
#'
#' @format A 'list' object
#'
#' @keywords datasets
#'
#' @references Lavoué, S., 2016. Was Gondwanan breakup the cause of the intercontinental distribution of Osteoglossiformes? A time-calibrated phylogenetic test combining molecular, morphological, and paleontological evidence. Molecular Phylogenetics and Evolution, 99, pp.34-43.
#' @examples
#' data(Lavoue2016)
#' combined_data <- cat_data(Lavoue2016$standard, Lavoue2016$dna, use.part.info = F) 
#' ## Visualise 
#' if(!require("BiocManager", quietly = TRUE)){
#' 	install.packages("BiocManager")
#' }
#' BiocManager::install("ComplexHeatmap")
#' library(ComplexHeatmap)
#' df <- t(as.data.frame(combined_data))
#' df <- tolower(df)
#' df[df == "?"] <- NA
#' df[df == "-"] <- NA
#' df[df == "n"] <- NA
#' states <- as.character(na.omit(unique(as.character(df))))
#'cols <- rep("x", length(states))
#'names(cols) <- states
#'morph_states <- suppressWarnings(which(is.na(as.numeric(names(cols))) == F))
#'mol_states <- which(names(cols) %in% c("a", "c", "g", "t"))
#'cols[morph_states] <- hcl.colors(n = length(morph_states), palette = "Hawaii")
#'cols[mol_states] <- hcl.colors(n = length(mol_states), palette = "zissou")
#'cols[-c(mol_states, morph_states)] <- hcl.colors(n = length(cols[-c(mol_states, morph_states)]), palette = "Cividis")
#'Heatmap(df[,1:1500], row_names_side = "left", col = cols, na_col = "white", name = "states")
#' 
"Lavoue2016"
################### Functions ###################

#' Read nexus data file. 
#'
#' This function is a modification of the ape function read.nexus.data(). The function has been modified so that it can read MrBayes files with partitions comprising different data types (e.g. datatype=mixed). 
#' @param file a file name specified by either a variable of mode character, or a double-quoted string.
#' @param use.part.info logical. If true, matrices of different data types will be stored separately. If true, all data types will be concatenated. 
#' @return a list object of class nexdat or multi_nexdat.  
#' @details The read.nexus.data() function from ape has been modified so that it can read MrBayes mixed data files, either as a single matrix with data types concatenated or as a list of matrices organised by datatype. In addition to this change, the function has been modified so that it can handle spaces in polymorphic or uncertain characters (e.g. (0 1), {0 1}). Due to this change, the function no longer checks if there are spaces in taxon names. Please ensure that spaces in taxon names are removed prior to using this function as they will cause the data to be read incorrectly.   
#' @examples
#' ## Use read.nexdat to read a data file in NEXUS format into object x
#' ## Not run: x <- read_nexdat("file.nex")
#' @export

read_nexdat <- function (file, use.part.info = F)
{
    "replace.single.quotes" <- function(x, start = 1L){
        z <- unlist(gregexpr("'", x))
        if (length(z) %% 2) {
            #warning("wrong number of single quotes around labels")
            warning(paste0("odd number of single quotes (", length(z), "): label(s) unchanged"))
            return(x)
        }
        i <- 1
        while (i < length(z)) {
            tmp <- substr(x, z[i], z[i + 1])
            substr(x, z[i], z[i + 1]) <- gsub("\\s+", "_", tmp)
            i <- i + 2
        }
        gsub("'", "", x)
    }

    "find.ntax" <- function(x) {
        for (i in 1:NROW(x)) {
            if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        ntax
    }
    "find.nchar" <- function(x) {
        for (i in 1:NROW(x)) {
            if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        nchar
    }
    "find.matrix.line" <- function(x) {
        for (i in 1:NROW(x)) {
            if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
                matrix.line <- as.numeric(i)
                break
            }
        }
        matrix.line
    }

    "find.datatype.line" <- function(x) {
      for (i in 1:NROW(x)) {
        if (any(f <- grep("format datatype", x[i], ignore.case = TRUE))) {
          matrix.line <- as.numeric(i)
          break
        }
      }
      matrix.line
    }


    "get.partition.info" <- function(x){
      part.info <- data.frame()
      x <- strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]), ",")[[1]]
      for(i in 1:length(x)){
        str_info <- stringr::str_extract_all(x[[i]], "\\b\\w+\\b")[[1]]
        newrow <- data.frame("partition" = tolower(str_info[[1]]), "start" = as.numeric(str_info[[2]]), "end" = as.numeric(str_info[[3]]))
        part.info <- rbind(part.info, newrow)
      }
      return(part.info)
    }

    "trim.whitespace" <- function(x) {
        gsub("\\s+", "", x)
    }
    "trim.semicolon" <- function(x) {
        gsub(";", "", x)
    }
    "is.poly.start" <- function(x) {
        return("(" == x || "{" == x)
    }
    "is.poly.end" <- function(x) {
        return(")" == x || "}" == x)
    }
    "get.polymorphism" <- function(x) {
        position <- 1
        while (position <= length(x)) {
            if (is.poly.start(x[position])) {
                poly_end <- position + 1
                while (!is.poly.end(x[poly_end])) {
                  poly_end <- poly_end + 1
                  if (is.poly.start(x[poly_end]) || poly_end >
                    length(x)) {
                    stop("missing closing bracket for a polymorphism at position ",
                      position)
                  }
                }
                x[position] <- paste0(x[(position + 1):(poly_end -
                  1)], collapse = "/")
                x <- x[-c((position + 1):poly_end)]
            }
            position <- position + 1
        }
        return(x)
    }
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
        comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0
    single.quotes <- grepl("'", X)
    if (any(single.quotes)) {
        to.replace <- which(single.quotes)
        for (j in to.replace) {
            X[[j]] <- replace.single.quotes(X[[j]])
        }
    }
    for (j in start.reading:NROW(X)) {
        Xj <- trim.semicolon(X[j])
        if (Xj == "") {
            break
        }
        if (any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE,
            ignore.case = TRUE))) {
            break
        }
        Xj <- gsub("\t", " ", Xj)
        Xj <- paste(sub(" .*", "", Xj), "  ",gsub(" ", "", sub("^\\S+\\s+", '', Xj), fixed = TRUE) )
        ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
        if (length(ts) > 2) {
            stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
        }
        if (length(ts) != 2) {
            stop("nexus parser failed to read the sequences (ts!=2)")
        }
        Seq <- trim.whitespace(ts[2])
        Name <- trim.whitespace(ts[1])
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
        if (any(l <- grep(nAME, names(Obj)))) {
            tsp <- strsplit(Seq, NULL)[[1]]
            if (any("(" %in% tsp || "{" %in% tsp)) {
                tsp <- get.polymorphism(tsp)
            }
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[l]][p] <- tsp[k]
                chars.done <- k
            }
        } else {
            names(Obj)[i] <- Name
            tsp <- strsplit(Seq, NULL)[[1]]
            if (any("(" %in% tsp || "{" %in% tsp)) {
                tsp <- get.polymorphism(tsp)
            }
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[i]][p] <- tsp[k]
                chars.done <- k
            }
        }
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax) {
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            if (tot.nchar == nchar * ntax) {
                print("ntot was more than nchar*ntax")
                break
            }
            pos <- tot.nchar
        } else {
            i <- i + 1
        }
    }
    if (tot.ntax != 0) {
        cat("ntax:", ntax, "differ from actual number of taxa in file?\n")
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            cat(names(Obj[i]), "has", length(Obj[[i]]), "characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
    Obj <- lapply(Obj, tolower)
    if(use.part.info == T){
      part_line <- find.datatype.line(X)
      p.info <- get.partition.info(X[[part_line]])
      Obj.info <- list()
      for(i in 1: nrow(p.info)){
        Obj.info[[i]] <- lapply(Obj, function(x) x[p.info[i,2]:p.info[i,3]])
        class(Obj.info[[i]]) <- c("nexdat", "list")
      }
      names(Obj.info) <- p.info[,1]
      Obj <- Obj.info
      class(Obj) <- c("multi_nexdat", "list")
    } else {
      class(Obj) <- c("nexdat", "list")
    }
    Obj
}


#' Convert phyDat_to_nexdat.
#'
#' This function coverts phylogenetic data in phyDat format to list format, equivalent to the output of read.nexus.data() or read.nexdat(). 
#' @param x an object of class 'phyDat'.
#' @return an object of class 'nexdat'. 
#' @details Note that phyDat objects can be converted to matrix objects using the as.character() function.  
#' @examples
#' install.packages("phangorn")
#' library(phangorn)
#' data("Laurasiatherian")
#' x <- phyDat_to_nexdat(Laurasiatherian)
#' @export

phyDat_to_nexdat<- function(x){
	if(class(x) != "phyDat"){
		stop("data is not phyDat object!")
	}
	x.mat <-  as.character(x)
	x.list <- lapply(1:nrow(x.mat), function(i) x.mat[i,])
	names(x.list) <- rownames(x.mat)
    class(x.list) <- c("nexdat", "list")
	return(x.list)
}


#' Get ultrametric edge lengths
#'
#' Function to generate ultrametric edge lengths for tree plots.
#' @importFrom ape rtree Nedge plot.phylo node.depth.edgelength
#' @importFrom phangorn Descendants
#' @param tree a tree object of class 'phylo'.
#' @return a tree object of class 'phylo' with ultrametric edge lengths.
#' @details This function is equivalent to plot.tree(node.depth = 2) except that a vector of edge.lengths are included in the output 'phylo' object.
#' @examples
#' t1 <- rtree(30);
#' plot(t1);
#' t2 <- u_edge_lengths(t1);
#' plot(t2);
#' @export

u_edge_lengths <- function(tree){
	if(class(tree) != "phylo"){
		stop("tree is not 'phylo' object")
	}
	tree$edge.length <- rep(1, Nedge(tree))
	nheights <- node.depth.edgelength(tree)
	theight <- max(nheights)
	for(i in 1:Ntip(tree)){
		edge_id <- which(tree$edge[,2] == i)
		tree$edge.length[edge_id] <- tree$edge.length[edge_id] + (theight - nheights[i])
	}
	for(i in 1:Nnode(tree)){
   		kids <- unlist(Descendants(tree, i+Ntip(tree), type = "tips"))
    	clade_tip_L <- unlist(lapply(kids, function(x) tree$edge.length[[which(tree$edge[,2] == x)]]))
    	mintipL <- min(clade_tip_L)
    	if(mintipL > 1){
    		Ldiff <- mintipL -1
    		for(j in kids){
    			edge_j <- which(tree$edge[,2] == j)
    			tree$edge.length[edge_j] <- tree$edge.length[edge_j] - Ldiff
    		}
    		edge_i <- which(tree$edge[,2] == i+Ntip(tree))
    		tree$edge.length[edge_i] <- tree$edge.length[edge_i] + Ldiff
    	}
    }
	return(tree)
}


#' Calculate support values for a specific tree given a sample of trees
#'
#' Function to calculate support values for a specific tree (e.g. a consensus tree, maximum clade credibility tree, maximum likelihood tree etc.) given a sample of trees (e.g. posterior trees, bootstrap replicate trees, most parsimonious trees).
#' @param trees a tree object of class 'multiPhylo'.
#' @param con a tree object of class 'phylo'.
#' @param digits integer indicating the number of decimal places to round to.
#' @return A tree object of class 'phylo' with node labels equivalent to clade support.
#' @details This function is a wrapper for the ape function prop.clades.
#' @examples
#' best_tree <- rtree(20)
#' sample_trees <- list()
#' for(i in 1:100){ sample_trees[[i]] <- rSPR(best_tree, 1)}
#' class(sample_trees) <- "multiPhylo"
#' best_tree <- clade_support(sample_trees, best_tree)
#' plot(best_tree)
#' nodelabels(text = best_tree$node.label)
#' @export

clade_support <- function(trees, con, digits = 2){
	if(is.rooted(con) == F){
		stop("Consensus tree is not rooted!")
	}
	con$node.label <- round(prop.clades(con, trees)/length(trees), digits)
 	con$node.label[is.na(con$node.label)] <- 0
	return(con)
}

#' Root a tree below a user-specified node.
#'
#' Function to root a tree at the midpoint of the branch ancestral to the user-specified node. 
#' @usage root_node(tree, node)
#' @param tree a tree object of class 'phylo'.
#' @param node number of the node or tip descending from the target branch in tree$edge. Alternatively, the user can specify a numeric vector of node/tip numbers, in which case the function will attempt to root the tree on the edge descending from the most recent common ancestor of the specified nodes/tips.
#' @return A tree object of class 'phylo'.
#' @details This function is a wrapper for the phytools function reroot(), which is more flexible as you can specify a position on an edge at which to root the tree. However, 99 times out of 100, I want to root at the midpoint of an edge. This function does this automatically. It also allows for rooting a tree with no edge.lengths, which returns an error when attempted with reroot(). See also the castor function root_in_edge(), which is similar to reroot().
#' @examples
#' ## Example with edge lengths
#' t <- rtree(10)
#' plot(t, label.offset = 0.1)
#' tiplabels()
#' t2 <- root_node(t, 2)
#' plot(t2, label.offset = 0.1)
#' tiplabels()
#' ## Example without edge lengths
#' t3 <- t
#' t3$edge.length <- NULL
#' plot(t3, label.offset = 0.1)
#' tiplabels()
#' t4 <- root_node(t3, 2)
#' plot(t4, label.offset = 0.1)
#' tiplabels()
#' @export

root_node <- function(tree, node){
	if(length(node) > 1){
		if(is.rooted(tree) == F){
			stop("tree is unrooted!")
		}
		node <- mrca.phylo(tree, node = c(node))
	}
	edge.lengths <- TRUE
	if(length(tree$edge.length) == 0){
		edge.lengths <- FALSE
		tree$edge.length <- rep(1, nrow(tree$edge))
	}
	n.edge <- which(tree$edge[,2] == node)
	tree <- reroot(tree, node, position = tree$edge.length[[n.edge]]/2)
	if(edge.lengths == FALSE){
		tree$edge.length <- NULL
	}
	return(tree)
}


#' Root multiple trees below a user-specified node.
#'
#' Function to root all trees in a multiPhylo object at the midpoint of a branch ancestral to a user-specified node. 
#' @usage root_node(trees, node)
#' @param trees a tree object of class 'multiPhylo'.
#' @param node number of the node or tip descending from the target branch in tree$edge. Alternatively, the user can specify a numeric vector of node/tip numbers, in which case the function will attempt to root the tree on the edge descending from the most recent common ancestor of the specified nodes/tips.
#' @return a tree object of class 'phylo'.
#' @details see help for root_node(). This function assumes that tip labels are identical in each tree of the multiPhylo object. Internal node labels will differ between trees of different topologies, thus it is recommended that multiPhylo objects are rooted by specifying a tip or multiple tips.  
#' @examples
#' trees <- rmtree(100, 10)
#' trees2 <- root_multi(trees, 2)
#' @export

root_multi <- function(trees, node){
	trees <- .compressTipLabel(trees)
	if(length(node) > 1){
		if(any(is.rooted(trees) == F)){
			stop("some trees are unrooted!")
		}
	}
	trees <- lapply(trees, FUN = root_node, node = 1)
	class(trees) <- "multiPhylo"
	return(trees)
}


#' Collapse nodes below a threshold value using node labels.
#'
#' Function to collapse all nodes of a phylo object with numeric node labels that are less than user-specified threshold value. 
#' @param tree a tree object of class 'phylo'.
#' @param threshold a numeric value. All nodes with a numeric node.label less than the threshold will be collapsed. 
#' @return a tree object of class 'phylo'.
#' @details This function requires the phylo object to have numeric node labels specified (e.g. tree$node.label). 
#' @examples
#' ## simulate two groups of trees
#' set.seed(8)
#' t1 <- rtree(10)
#' trees1 <- t1
#' for(i in 1:49){ 
#'   trees1 <- c(trees1, rSPR(t1, 1))
#' }
#' t2 <- rSPR(t1, 3)
#' t2 <- multi2di(t2)
#' trees2 <- t2
#' for(i in 1:49){ 
#'   trees2 <- c(trees2, rSPR(t2, 1))
#' }
#' trees <- c(trees1, trees2)
#' ## get maximum clade credibility tree
#' mcc <- maxCladeCred(trees)
#' plot(mcc)
#' ## calculate clade_support values
#' mcc <- clade_support(trees, mcc, digits = 2)
#' nodelabels(text = mcc$node.label)
#' ## collapse nodes with less than 50% support 
#' mcc_collapsed <- collapse_thresh(mcc, 0.5)
#' plot(mcc_collapsed)
#' nodelabels(text = mcc_collapsed$node.label)
#' @export

collapse_thresh <- function(tree, threshold){
	if(is.numeric(tree$node.label) == F){
		stop("tree requires numeric node labels!")
	}
	if(length(tree$node.label) != Nnode(tree)){
		stop("node.label length does not equal the number of tree nodes!")
	}
	nl <- tree$node.label[which(tree$node.label >= threshold)]
	tree <- TreeTools::CollapseNode(tree, nodes =  which(tree$node.label < threshold) + Ntip(tree))
	tree$node.label <- nl
	return(tree)
}


#' Make perfectly symmetrical tree.
#'
#' Function to generate a perfectly symmetrical (balanced) tree with equal edge lengths and a height of 1. 
#' @param n number of tips to generate. Must be in the geometric sequence: a1=2, r=2 (e.g. 2, 4, 8, 16, 32...).
#' @return a tree object of class 'phylo'.
#' @examples
#' t <- make_sym_tree(64)
#' plot(t, show.tip.label = F)
#' @export

make_sym_tree <- function(n = 32){
	if(((log2(n))%%1==0) == FALSE){
		stop("n is not in the geometric sequence: a1= 2, r = 2")
	}
	taxa <- paste("t", 1:n, sep = "")
	taxa <- as.list(taxa)
	repeat{
		taxa2 <- list()
		j <- 0
		for(i in 1:(length(taxa)/2)){
			taxa2[[i]] <- paste("(", taxa[[(i + j)]], ",", taxa[[(i + j + 1)]], ")", sep = "")
			j <- j + 1
		}
		taxa <- taxa2
		if(length(taxa) == 1){
			break
		}
	}
	tree <- paste(taxa[[1]], ";", sep = "")
	tree <- read.newick(text = tree)
	tree$edge.length <- rep(1/ (log2(n)), nrow(tree$edge))
	return(tree)
}


#' Make perfectly asymmetrical tree.
#'
#' Function to generate a perfectly asymmetrical (unbalanced) tree with equal edge lengths and a height of 1. 
#' @param n number of tips to generate.
#' @return a tree object of class 'phylo'.
#' @examples
#' t <- make_asym_tree(64)
#' plot(t, show.tip.label = F)
#' @export

make_asym_tree <- function(n = 32){
	taxa <- paste("t", 1:n, sep = "")
	taxa <- as.list(taxa)
	bl <- 1/(n - 1)
	tree <- paste(taxa[[1]], ":", bl, sep = "")
	for(i in 2:n){
		tree <- paste("(", tree, ",", taxa[[i]], ":", (bl * (i-1)), ")", sep = "")
		if(i != n){
			tree <- paste(tree, ":", bl, sep = "")
		}
	}
tree <- paste(tree, ";", sep = "")
tree <- read.newick(text = tree)
return(tree)
}


#' Get tip priors.
#'
#' Function to extract tip priors for ancestral state estimation from a character vector, list or data.frame. 
#' @param morph a character vector, list or data.frame coded in 'standard' format (e.g. morphological data).
#' @param treat.as.observed logical. If true, uncertainty (e.g. missing data, state ambiguities) are treated as observed and assigned a probability of 1. If false, uncertainty is treated as not observed directly. Probability is split equally between ambiguous states. 
#' @param extar_state logical. If TRUE, inapplicable codings ('-') will be treated as an additional state.
#' @return A 'list' of numeric matrices of length n, where n is equal to the number of characters included in the input. Each matrix is size i * j, where i is the number of taxa and j is the number of states included in the nth character. Each cell represents the likelihood of observing state j in taxon i.  
#' @details This function converts categorical characters, coded in 'standard' format (e.g. '0', '1', '2', '-', '?', '0/1') into a 2D tip prior matrix, or list of 2D tip prior matrices. These matrices can be used instead of character vectors in some ancestral state estimation functions (e.g. castor: asr_mk_model(), phytools: fitMk()). The function assigns equal prior probability to observing any state if data are coded as missing ('?') or inapplicable ('-'). If data are coded as polymorphic or uncertain between specified states (e.g. '0/1'. '(0 1)'), prior probability is divided equally between possible states. If extra_state is set to TRUE, inapplicable characters are considered an additional state. This may be useful when dealing with datasets that employ an atomised coding strategy in which absence / presence of a trait is coded seperately from the condition of the trait (e.g. red, blue). In these cases, the coding of 'inapplicable' in the condition character is logically equivalent to absent.    
#' @examples
#' ## Load data
#' data(vert_data)
#' head(vert_data$morph)
#' ## Convert secondary osteon character to tip priors, interpreting inapplicable ('-') as an additional state (i.e. absent).
#' tp <- get_tip_priors(vert_data$morph[,2], extra_state = T)
#' colnames(tp[[1]]) <- c("bone absent", "secondary osteons absent", "secondary osteons present")
#' ## conduct ancestral state analysis using phytools.
#' fitER <- fitMk(tree = vert_data, x = tp[[1]], model = "ER")
#' plot(ancr(fitER), legend = "topleft")
#' @export

get_tip_priors <- function(morph, treat.as.observed = TRUE, extra_state = F){
	if(is.atomic(morph) == T){
		morph_df <- as.data.frame(morph)
	} else if(class(morph) == "list"){
		morph_df <- t(as.data.frame(morph))
	} else {
		morph_df <- morph
	}
	results <- list()
	for(n in 1:ncol(morph_df)){
		states <- list()
		inap <- numeric()
		uncertain <- numeric()
		for(i in 1: nrow(morph_df)){
			state_i <- morph_df[i,n]
			if(is.na(suppressWarnings(as.numeric(state_i)))){
				if(extra_state == T){
					if(state_i == "-"){
						inap <- c(inap, i)
					}
					if(state_i == "?"){
						uncertain <- c(uncertain, i)
					}
				}
				if(extra_state == F){
					if(state_i == "?"|| state_i == "-"){
						uncertain <- c(uncertain, i)
					} 
				}
				if(grepl("/", state_i, fixed = TRUE)){
					states[[i]] <- as.numeric(unlist(strsplit(state_i, "/"))) + 1
				} else if(grepl("(", state_i, fixed = TRUE)){
					state_i <- gsub("\\(", "",
					gsub("\\)", "",
					gsub(" ", "",
    				state_i)))
    				states[[i]] <- as.numeric(unlist(strsplit(state_i, ""))) + 1
    			} else if(grepl("{", state_i, fixed = TRUE)){
					state_i <- gsub("\\{", "",
					gsub("\\}", "",
					gsub(" ", "",
    				state_i)))
    				states[[i]] <- as.numeric(unlist(strsplit(state_i, ""))) + 1
    			}
			} else {
				states[[i]] <- as.numeric(state_i) + 1
			}
		}
		n_states <- max(na.omit(unlist(states)))
		if(extra_state == T){
			if(length(inap) > 0){
				n_states <- n_states + 1
				states <- lapply(states, function(x){ x + 1})
				for(j in inap){
					states[j] <- 1
				}
			}
		}
		if(n_states == -Inf){
			tip_priors <- NA
		} else {
			for(j in uncertain){
				states[[j]] <- 1:n_states
			}
			tip_priors <- matrix(0, nrow(morph_df), n_states)
			for(k in 1: nrow(morph_df)){
				for(l in states[[k]]){
					if(treat.as.observed == T){
						tip_priors[k, l] <- 1
					} else {
						tip_priors[k, l] <- 1 / length(states[[k]])
					}
				}
			}
			rownames(tip_priors) <- row.names(morph_df)
			colnames(tip_priors) <- 1:ncol(tip_priors)
			results[[n]] <- tip_priors
		}
	}
	return(results)
}

########## FIX ABOVE! Does not work for nexdat! ########## 

#' Amalgamate tip priors.
#'
#' Function to amalgamate tip priors for ancestral state estimation. 
#' @param tp_1 an object of class 'data.frame' compising the tip priors of the controlling character.
#' @param tp_2 an object of class 'data.frame' compising the tip priors of the dependent character.
#' @return An object of class 'data.frame' compising the amalgamated tip priors of the controlling and dependent character. 
#' @details This function amalgamates two seperate sets of tip priors, corresponding to the controlling and dependent character respectively, into a single set of tip priors corresponding to a Structured Markov Model (SMM) equipped with hidden states (see Tarasov 2019).
#' @references Tarasov, S., 2019. Integration of anatomy ontologies and evo-devo using structured Markov models suggests a new framework for modeling discrete phenotypic traits. Systematic biology, 68(5), pp.698-716.
#' @examples
#' ## Load data
#' data(vert_data)
#' head(vert_data$morph)
#' tp <- get_tip_priors(vert_data$morph, extra_state = F)
#' tp_SMM <- amal_tip_priors(tp[[1]], tp[[2]])
#' colnames(tp_SMM) <- c("aA", "aP", "pA", "pP")
#' ## Get SMM using rphenoscate
#' devtools::install_github("uyedaj/rphenoscate")
#' library(rphenoscate)
#' c1 <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("a", "p"), c("a", "p")) )
#' c2 <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("A", "P"), c("A", "P")) )
#' SMM_ind <- amaSMM_2Q(c1, c2, controlling.state = NULL, diag.as = 0)
#' ## conduct ancestral state analysis using phytools and plot.
#' fitSMM <- fitMk(tree = vert_data, x = tp_SMM, model = SMM_ind)
#' ancSMM <- ancr(fitSMM)
#' plot(ancSMM, legend = "topleft")
#' ## Collapse hidden states to observed states and plot
#' ancSMM$ace <- as.matrix(data.frame("bone absent" = ancSMM$ace[,1] + ancSMM$ace[,2], "secondary osteons absent" = ancSMM$ace[,3], "secondary osteons present" = ancSMM$ace[,4]))
#' attributes(ancSMM)$data <- as.matrix(data.frame("bone absent" = tp_SMM[,1] + tp_SMM[,2], "secondary osteons absent" = tp_SMM[,3], "secondary osteons present" = tp_SMM[,4]))
#' plot(ancSMM, legend = "topleft")
#' @export

amal_tip_priors <- function(tp_1, tp_2){
	if(all(rownames(tp_1) != rownames(tp_2))){
    	stop("rownames for tp_1 and tp_2 differ!")
	}
	if(is.null(colnames(tp_1))){
		colnames(tp_1) <- 1:ncol(tp_1)
	}
	if(is.null(colnames(tp_2))){
		colnames(tp_2) <- 1:ncol(tp_2)
	}
    result <- matrix(0, nrow(tp_1), ncol(tp_1) * ncol(tp_2))
    col_n <- 1
    c_names <- character()
    for(i in 1:ncol(tp_1)){
        result[, col_n:(col_n -1 + ncol(tp_2))] <- tp_1[,i] * tp_2
        c_names <- c(c_names, paste(colnames(tp_1)[[i]], colnames(tp_2), sep = ""))
        col_n <- col_n + ncol(tp_2)
    }
    colnames(result) <- c_names
	rownames(result) <- rownames(tp_1)
    return(result)
}

#' Get a contrast matrix.
#'
#' Function to automatically extract a contrast matrix from a character vector, list or data.frame of standard categorical data (e.g. morphological data).
#' @param morph a character vector, list or data.frame coded in 'standard' format (e.g. morphological data).
#' @return A contrast matrix used by the phangorn function phyDat() to interpret standard categorical data. 
#' @details This function extracts a contrast matrix from categorical characters, coded in 'standard' format (e.g. '0', '1', '2', '-', '?', '0/1'). The contrast matrix is used by the phangorn phyDat() function in order to interpret symbols denoting character state ambiguity (e.g. polymorphism, missing data etc.). Data in phyDat format can be used to estimate phylogeny in R (e.g. phangorn::parsimony(), TreeSearch::MaximizeParsimony()).  
#' @examples
#' ## Load data
#' install.packages("TreeSearch")
#' library(TreeSearch)
#' data("inapplicable.datasets")
#' morph <- inapplicable.datasets$Vinther2008
#' ## get_contrast matrix
#' cm <- get_contrast(morph)
#' ## Convert data to phyDat format
#' pdat <- phyDat(data = morph, type = "USER", contrast = cm)
#' ## Run parsimony analysis
#' trees <- MaximizeParsimony(dataset = pdat)
#' ## Plot consensus tree and support values
#' con_tree <- consensus(trees, p = 0.5)
#' con_tree <- root_node(con_tree, 1)
#' con_tree$node.label[-1] <- round(as.numeric(con_tree$node.label[-1]), 2)
#' plot(ladderize(con_tree), cex = 0.5, no.margin = T)
#' nodelabels(text = con_tree$node.label)
#' @export

get_contrast <- function(morph){
	if(is.atomic(morph) == T){
		morph_df <- as.data.frame(morph)
	} else if(class(morph) == "list"){
		morph_df <- t(as.data.frame(morph))
	} else {
		morph_df <- morph
	}
	u.states <- unique(unlist(as.vector(morph_df)))
	states <- list()
	uncertain <- numeric()
	for(i in 1:length(u.states)){
		state_i <- u.states[[i]]
		if(is.na(suppressWarnings(as.numeric(state_i)))){
			if(state_i == "?"| state_i == "-"){
				uncertain <- c(uncertain, i)
			} else if(grepl("/", state_i, fixed = TRUE)){
				states[[i]] <- as.numeric(unlist(strsplit(state_i, "/"))) + 1
			} else if(grepl("(", state_i, fixed = TRUE)){
				state_i <- gsub("\\(", "",
				gsub("\\)", "",
				gsub(" ", "",
    			state_i)))
    			states[[i]] <- as.numeric(unlist(strsplit(state_i, ""))) + 1
    		} else if(grepl("{", state_i, fixed = TRUE)){
				state_i <- gsub("\\{", "",
				gsub("\\}", "",
				gsub(" ", "",
    			state_i)))
    			states[[i]] <- as.numeric(unlist(strsplit(state_i, ""))) + 1
    		}
		} else {
			states[[i]] <- as.numeric(state_i) + 1
		}
	}
	n_states <- max(na.omit(unlist(states)))
	if(n_states == -Inf){
		contrast <- NA
	} else {
		for(j in uncertain){
			states[[j]] <- 1:n_states
		}
		contrast <-  matrix(0, length(states), n_states)
		rownames(contrast) <- u.states
		colnames(contrast) <- 1:n_states - 1
		for(k in 1:length(states)){
			contrast[k, states[[k]]] <- 1
		}
	}
	return(contrast)
}

#' Generate an exponential node calibration
#'
#' Function to generate an exponential node calibration for Bayesian molecular clock analysis. 
#' @usage exp_calib(age_min, age_max, position, xlim = NULL, ylim = NULL)
#' @param age_min hard minimum fossil calibration.
#' @param age_max soft maximum fossil calibration.
#' @param position numeric value between 0 and 1, which determines the percentile of the soft maximum age. 
#' @param xlim numeric vector defining the dimensions of the x axis. If NULL, xlim is determined automatically.
#' @param ylim numeric vector defining the dimensions of the y axis. If NULL, ylim is determined automatically.    
#' @return a histogram of the prior distribution with labelled soft maximum and hard minimum fossil calibrations. The function also prints the mrbayes command for defining the calibration. 
#' @references Benton, M.J., Donoghue, P.C., Asher, R.J., Friedman, M., Near, T.J. and Vinther, J., 2015. Constraints on the timescale of animal evolutionary history.
#' @examples
#' ##vertebrate calibration from Benton et al (2015): hard min = 457.5, soft max = 636.1.
#' par(mar = c(2.2,2.2,1.5,1.5))
#' par(mfrow=c(2,1))
#' ##plot calibration assuming a very good fossil record (i.e. low probability the age of the node is older than the soft maximum)
#' exp_calib(457.5, 636.1, 0.99, xlim = c(0, 1000), ylim = c(0, 0.025))
#' ##plot calibration assuming a poor fossil record (i.e. high probability the age of the node is older than the soft maximum)
#' exp_calib(457.5, 636.1, 0.75, xlim = c(0, 1000), ylim = c(0, 0.025))
#' @export

exp_calib <- function(age_min, age_max, position, xlim = NULL, ylim = NULL){
	if(position > 1 | position <= 0){
		stop("position must be between 0 and 1")
	}
	if(age_min >= age_max){
		stop("age_min must be < age_max")
	}
	diff <- age_max - age_min
	scale_param <- (log((position-1)*-1) /diff) *-1
	mean_calib <- 1/scale_param
	x <- rexp(100000, rate = scale_param) + age_min
	age <- seq(from = age_min, to = max(x), length.out = 100)
	d <- dexp(x = seq(from = 0, to = max(x) - age_min, length.out = 100), rate = scale_param)
	df <- data.frame(age = age, d = d)
	df <- rbind(data.frame(age = c(0, age_min - 0.00001), d = c(0, 0)), df)
	df <- rbind(df, data.frame(age = max(x) * 100, d = 0))
	if(is.null(xlim)){
		xlim <- c(0, max(x))
	}
	if(is.null(ylim)){
		ylim <- c(0, max(d))
	}
	plot(df, type = "l",  xlim=xlim, ylim=ylim, main = paste("Hard minimum = ", age_min, ", soft maximum = ", age_max, ", percentile = ", position * 100, sep = ""),  xlab= "Age", cex.main=0.9, ylab = "Density")
	abline(v=age_max, col="blue")
	abline(v=age_min, col="red")
	return(cat("offsetexp(", age_min, ", ", age_min + mean_calib, ");", "\n", sep = ""))
}

#' Generate a gamma-shaped node calibration
#'
#' Function to generate a gamma-shaped node calibration for Bayesian molecular clock analysis. 
#' @usage gamma_calib(age_min, age_max, position, xlim = NULL, ylim = NULL)
#' @param age_min hard minimum fossil calibration.
#' @param age_max soft maximum fossil calibration.
#' @param shape shape parameter for the gamma distribution.
#' @param position numeric value between 0 and 1, which determines the percentile of the soft maximum age. 
#' @param xlim numeric vector defining the dimensions of the x axis. If NULL, xlim is determined automatically.
#' @param ylim numeric vector defining the dimensions of the y axis. If NULL, ylim is determined automatically.    
#' @return a histogram of the prior distribution with labelled soft maximum and hard minimum fossil calibrations. The function also prints the mrbayes command for defining the calibration. 
#' @references Benton, M.J., Donoghue, P.C., Asher, R.J., Friedman, M., Near, T.J. and Vinther, J., 2015. Constraints on the timescale of animal evolutionary history.
#' @examples
#' ##vertebrate calibration from Benton et al (2015): hard min = 457.5, soft max = 636.1.
#' par(mar = c(2.2,2.2,1.5,1.5))
#' par(mfrow=c(3,1))
#' ##plot calibration assuming a low probability the age of the node is older than the soft maximum age and a low probability the age of the node is older than the hard-minimum age.
#' gamma_calib(457.5, 636.1, 1, 0.99, xlim = c(0, 1000), ylim = c(0, 0.025))
#' ##plot calibration assuming a low probability the age of the node is older than the soft maximum age and a high probability the age of the node is older than the hard-minimum age.
#' gamma_calib(457.5, 636.1, 3, 0.99, xlim = c(0, 1000), ylim = c(0, 0.025))
#' ##plot calibration assuming a high probability the age of the node is older than the soft maximum age and a high probability the age of the node is older than the hard-minimum age.
#' gamma_calib(457.5, 636.1, 3, 0.75, xlim = c(0, 1000), ylim = c(0, 0.025))
#' @export

gamma_calib <- function(age_min, age_max, shape, position, xlim = NULL, ylim = NULL){
	if(position > 1 | position <= 0){
		stop("position must be between 0 and 1")
	}
	if(age_min >= age_max){
		stop("age_min must be < age_max")
	}
	diff <- age_max - age_min
	x <- c(runif(100, 0.000001, 0.00001) ,runif(100, 0.00001, 0.0001), runif(100, 0.0001, 0.001), runif(100, 0.001, 0.01), runif(100, 0.01, 0.1), runif(100, 0.1, 1), runif(100, 1, 10), runif(100, 10, 100), runif(100, 100, 1000), runif(100, 1000, 10000))
	tests <- unlist(lapply(x, function(y) pgamma(q = diff, shape = shape, rate = y)))
	tests <- data.frame(x = x, tests = tests)
	tests <-tests[order(tests$tests),]

	l_bound <- tests$x[[max(which(as.numeric(tests$tests) <= position))]]
	u_bound <- tests$x[[min(which(as.numeric(tests$tests) >= position))]]

	for(i in 1:1000){
		x <- runif(10000, l_bound, u_bound)
		tests <- unlist(lapply(x, function(y) pgamma(q = diff, shape = shape, rate = y)))
		tests <- data.frame(x = x, tests = tests)
		tests <-tests[order(tests$tests),]

		l_bound <- tests$x[[max(which(as.numeric(tests$tests) <= position))]]
		u_bound <- tests$x[[min(which(as.numeric(tests$tests) >= position))]]
		if(round(l_bound, 8) == round(u_bound, 8)){
			break
		}
	}
	rate_param <- mean(c(l_bound, u_bound))
	mean_calib <- shape/rate_param
	x <- rgamma(100000, shape = shape, rate = rate_param) + age_min
	if(is.null(xlim)){
		xlim <- c(0, max(x))
	}
	age <- seq(from = age_min, to = max(x), length.out = 100)
	d <- dgamma(x = seq(from = 0, to = max(x) - age_min, length.out = 100), shape = shape, rate = rate_param)
	df <- data.frame(age = age, d = d)
	if(is.null(ylim)){
		ylim <- c(0, max(d))
	}
	df <- rbind(data.frame(age = c(0, age_min - 0.00001), d = c(0, 0)), df)
	df <- rbind(df, data.frame(age = max(x) * 100, d = 0))
	plot(df, type = "l",  xlim=xlim, ylim = ylim, main = paste("Hard minimum = ", age_min, ", soft maximum = ", age_max, ", shape = ", shape, ", percentile = ", position * 100, sep = ""),  xlab= "Age", cex.main=0.9, ylab = "Density")
	abline(v=age_max, col="blue")
	abline(v=age_min, col="red")
	sd_dist <- sqrt(shape)/rate_param
	return(cat("offsetgamma(", age_min, ", ", age_min + mean_calib, ", ", sd_dist, ");", "\n", sep = ""))
}


#' remove partition information
#'
#' Function to remove partition information from a 'multi_nexdat' object.
#' @param x an object of class 'multi_nexdat'.
#' @return an object of class 'nexdat'.
#' @details This function removes partition information from a 'multi_nexdat' object, thereby converting it into a 'nexdat' object. Sequences from seperate partitions are conctenated.  
#' @examples
#' data(Lavoue2016)
#' partitioned_data <- cat_data(Lavoue2016$standard, Lavoue2016$dna, use.part.info = T)
#' concatenated_data <- remove_part_info(partitioned_data)  
#' 
#' @export 


remove_part_info <- function(x){
    if(any(class(x) == "multi_nexdat") == F){
        stop("Object is not a multi_nexdat!")
    }
    x2 <- list()
    taxa <- sort(unique(unlist(lapply(x, names))))
    for(j in taxa){
        char <- character()
        for(i in names(x)){
            char <- c(char, x[[i]][[j]])
        }
        x2[[j]] <- char
    }
    class(x2) <- c("nexdat", "list")
    return(x2) 
}


#' Combine phylogenetic data
#'
#' Function to combine lists of phylogenetic data (e.g. those output by the ape function read.nexus.data()).  
#' @param ... two or more objects of class 'list' comprising named vectors of phylogenetic characters. 
#' @param use.part.info logical. If TRUE, the function resturns a seperate list for each input list. 
#' @param part.names a character vector of partition names e.g. 'dna', 'standard', 'protein'. If not supplied, the function will attempt to guess the correct partition type. 
#' @return A list, or multiple lists, of phylogenetic data of length T, where T equals the sum of unique taxa in the input lists.
#' @details This function can be used to combine seperate lists of phylogenetic characters with non-overlapping taxa. If a taxon does not occur in all input lists, it will be coded as missing '?' for those input list(s). Can be used, for example, to combine molecular and morphological datasets with non-overlapping taxa. Note that taxa must be named identically in each input for the data to be combined succesfully.    
#' @examples
#' data(Lavoue2016)
#' combined_data <- cat_data(Lavoue2016$standard, Lavoue2016$dna, use.part.info = F) 
#' ## Visualise 
#' if(!require("BiocManager", quietly = TRUE)){
#' 	install.packages("BiocManager")
#' }
#' BiocManager::install("ComplexHeatmap")
#' library(ComplexHeatmap)
#' df <- t(as.data.frame(combined_data))
#' df <- tolower(df)
#' df[df == "?"] <- NA
#' df[df == "-"] <- NA
#' df[df == "n"] <- NA
#' states <- as.character(na.omit(unique(as.character(df))))
#'cols <- rep("x", length(states))
#'names(cols) <- states
#'morph_states <- suppressWarnings(which(is.na(as.numeric(names(cols))) == F))
#'mol_states <- which(names(cols) %in% c("a", "c", "g", "t"))
#'cols[morph_states] <- hcl.colors(n = length(morph_states), palette = "Hawaii")
#'cols[mol_states] <- hcl.colors(n = length(mol_states), palette = "zissou")
#'cols[-c(mol_states, morph_states)] <- hcl.colors(n = length(cols[-c(mol_states, morph_states)]), palette = "Cividis")
#'Heatmap(df[,1:1500], row_names_side = "left", col = cols, na_col = "white", name = "states")
#' 
#' @export

cat_data <- function(..., use.part.info = F, part.names = NULL){
	partitions <- list(...)
	taxa <- sort(unique(unlist(lapply(partitions, names))))
	new_data <- list()
	for(i in 1:length(partitions)){
        partitions[[i]] <- lapply(partitions[[i]], tolower)
		new_data[[i]] <- list()
        # attempt to guess partition names

        if(is.null(part.names)){
            px <- unlist(partitions[[i]])
            # remove missing data
            px <- px[-which(px == "?"| px == "n" | px == "x" | px == "-")]
            px_l <- length(px)
            if(length(suppressWarnings(na.omit(as.numeric(px))))/px_l > 0.9){
                names(new_data)[[i]] <- "standard"
            } else if(length(px[which(px == "a" | px == "c" | px == "g" | px == "t" | px == "u")])/px_l > 0.9){
                names(new_data)[[i]] <- "dna"
                if(length(px[which(px == "u")]) > 0){
                    names(new_data)[[i]] <- "rna"
                }
            } else {
                names(new_data)[[i]] <- "protein"
            }
        }
        names(partitions) <- tolower(part.names)
		for(j in taxa){
			if(is.null(partitions[[i]][[j]])){
				new_data[[i]][[j]] <- rep("?", length(partitions[[i]][[1]]))
			} else { 
				new_data[[i]][[j]] <- partitions[[i]][[j]]
			}
		}
	}
    class(new_data) <- c("multi_nexdat", "list")
    if(use.part.info == F){
        new_data <- remove_part_info(new_data)  
    }
    return(new_data)
}


#' Calculate shared bipartitions between trees.
#'
#' Function to calculate the number of shared bipartitions between a reference tree and a single or multiple comparsion tree(s). 
#' @param tree1 an object of class 'phylo'.
#' @param tree2 an object of class 'phylo' or 'multiPhylo'.
#' @details This is a symmetric measure and is closely related to the Robinson-Foulds distance. 
#' @return A numeric value equal to the number of shared bipartitions between the two trees. 
#' @references Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.
#' @examples
#' t1 <- rtree(20)
#' t2 <- rSPR(t1, 4)
#' SB(t1, t2)
#' @export

SB <- function(tree1, tree2){
	tree1 <- unroot(tree1)
	tree2 <- unroot(tree2)
	RF <- RF.dist(tree1, tree2)
	part_all <- (Nnode(tree1) -1) + (Nnode(tree2) -1)
	return((part_all - RF)/2)
}


#' Calculate unique bipartitions in a tree.
#'
#' Function to calculate the number of bipartitions in a reference tree, that are not found in a comparsion tree.  
#' @param tree1 an object of class 'phylo'.
#' @param tree2 an object of class 'phylo'.
#' @details This is an asymmetric measure and is closely related to the Robinson-Foulds distance. UB(tree1, tree2) + UB(tree2, tree1) = RF.dist. 
#' @return A numeric value equal to the number of  bipartitions between the two trees. 
#' @references Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.
#' @examples
#' t1 <- rtree(20)
#' t2 <- rSPR(t1, 4)
#' UB(t1, t2)
#' @export

UB <- function(tree1, tree2){
	tree1 <- unroot(tree1)
	tree2 <- unroot(tree2)
	RF <- RF.dist(tree1, tree2)
	part1 <- (Nnode(tree1) -1)
	part2 <- (Nnode(tree2) -1)
	part_all <- (Nnode(tree1) -1) + (Nnode(tree2) -1)
	SB <- ((part_all - RF)/2)
	return(part1 - SB)
}


#' Calculate the likelihood of the 'total_garbage' model
#'
#' Function to calculate the likelihood of the 'total_garbage' model from Harmon (2019).  
#' @importFrom ape rtree 
#' @importFrom phytools sim.history fitMk
#' @param tip_states an object of class 'character'.
#' @details The 'total_garbage' model from Harmon (2019) is designed to test if your data provide no informatio about historical patterns of character change. It is equivalent to drawing states at random from a hat. Consider a set of tip_states and an evolutionary model with transition matrix Q. If likeihood(Q) given tip_states is similar to likelihood(total_garbage) given tip_states, then tip_states contains little historical information from which to infer state changes. However if likeihood(Q) given tip_states is greater than likelihood(total_garbage) given tip_states, tip_states does contain some information about historical character changes.    
#' @return A numeric value equal to the likelihood of the 'total_garbage' model given tip_states. 
#' @references Harmon, L., 2019. Phylogenetic comparative methods: learning from trees. https://lukejharmon.github.io/pcm/
#' @examples
#' ## simulate some data with slow evolutionary rates
#' tree <- rtree(50)
#' Q <- matrix(c(-0.2, 0.2, 0.2, -0.2), nrow = 2, ncol = 2)
#' true_history <- sim.history(tree, Q)
#' ## Estimate Q.
#' fitMK <- fitMk(tree, x = true_history$states, model = "ER")
#' ## Compare likelihoods. 
#' fitMK$logLik
#' total_garbage(true_history$states)
#' fitMK$logLik > total_garbage(true_history$states)
#' ## simulate some data with very fast evolutionary rates
#' tree <- rtree(50)
#' Q <- matrix(c(-100, 100, 100, -100), nrow = 2, ncol = 2)
#' true_history <- sim.history(tree, Q)
#' ## Estimate Q.
#' fitMK <- fitMk(tree, x = true_history$states, model = "ER")
#' ## Compare likelihoods. 
#' fitMK$logLik
#' total_garbage(true_history$states)
#' fitMK$logLik > total_garbage(true_history$states)
#' @export

total_garbage <- function(tip_states){
	n <- length(tip_states)
	state_count <- table(tip_states)
	res <- numeric()
	for(i in 1:length(state_count)){
		res <- c(res, as.numeric(state_count[[i]] %*% log(state_count[[i]]/n)))
	}
	return(sum(res))
}


#' Calculate AIC, AICc and BIC
#' 
#' Function to calculate the AIC, AICc and BIC of a model. 
#' @param loglik numeric. The log likelihood of the model. 
#' @param k numeric. The number of free parameters of the model. 
#' @param n numeric. The sample size (i.e. the number of observations). 
#' @return A numeric vector comprising the AIC, AICc and BIC. 
#' @references Akaike, H., 1998. Information theory and an extension of the maximum likelihood principle. In Selected papers of hirotugu akaike (pp. 199-213). New York, NY: Springer New York. \cr 
#' Sugiura, N., 1978. Further analysis of the data by akaike's information criterion and the finite corrections: Further analysis of the data by akaike's. Communications in Statistics-theory and Methods, 7(1), pp.13-26. \cr 
#' Hurvich, C.M. and Tsai, C.L., 1991. Bias of the corrected AIC criterion for underfitted regression and time series models. Biometrika, 78(3), pp.499-509.  \cr 
#' Schwarz, G., 1978. Estimating the dimension of a model. The annals of statistics, pp.461-464.
#' @examples
#' ## simulate some binary data with the ARD model.
#' tree <- rtree(50)
#' Q <- structure(c(-0.2, 0.8, 0.2, -0.8), dim = c(2L, 2L))
#' true_history <- sim.history(tree, Q)
#' ## Estimate Q using ER and ARD models.
#' fitER <- fitMk(tree, x = true_history$states, model = "ER")
#' fitARD <- fitMk(tree, x = true_history$states, model = "ARD")
#' ## Calculate information criteria. 
#' get_info(fitER$logLik, 1, 50)
#' get_info(fitARD$logLik, 2, 50)
#' @export 

get_info <- function(loglik, k, n){
	AIC <- 2*k-2*loglik
	AICc <- AIC + (2*k^2 + 2*k)/ (n - k - 1)
	BIC <- k*log(n)-2*loglik
	x <- c("AIC" = AIC, "AICc" = AICc, "BIC" = BIC)
	return(x)
}


#' Compare models using model weights
#'
#' Function to compare models based on AIC, AICc or BIC weights. 
#' @param loglik numeric. A numeric vector comprising log likelihoods of competing model. 
#' @param k a numeric vector comprising the number of parameters of competing model. 
#' @param n a numeric vector comprising the number of samples of competing model.
#' @param method either "AIC", "AICc" or "BIC".  
#' @param m.names an optional character vector of model names.  
#' @return A data.frame comprising the log likelihood, number of parameters, number of samples, information criterion and model weight of each model.     
#' @references Akaike, H., 1998. Information theory and an extension of the maximum likelihood principle. In Selected papers of hirotugu akaike (pp. 199-213). New York, NY: Springer New York. \cr 
#' Sugiura, N., 1978. Further analysis of the data by akaike's information criterion and the finite corrections: Further analysis of the data by akaike's. Communications in Statistics-theory and Methods, 7(1), pp.13-26. \cr 
#' Hurvich, C.M. and Tsai, C.L., 1991. Bias of the corrected AIC criterion for underfitted regression and time series models. Biometrika, 78(3), pp.499-509.  \cr 
#' Schwarz, G., 1978. Estimating the dimension of a model. The annals of statistics, pp.461-464.
#' @examples
#' ## simulate some binary data with the ARD model.
#' tree <- rtree(50)
#' Q <- structure(c(-0.2, 0.8, 0.2, -0.8), dim = c(2L, 2L))
#' true_history <- sim.history(tree, Q)
#' ## Estimate Q using ER and ARD models.
#' fitER <- fitMk(tree, x = true_history$states, model = "ER")
#' fitARD <- fitMk(tree, x = true_history$states, model = "ARD")
#' ## Calculate model weights. 
#' logL <- c(fitER$logLik, fitARD$logLik)
#' m.wts <- comp_models(loglik = logL, k = c(1, 2), n = 50, method = "BIC", m.names = c("ER", "ARD")) 
#' m.wts 
#' 
#' @export 

comp_models <- function(loglik, k, n, method = "AIC", m.names = NULL){
  x <- which(c("AIC", "AICc", "BIC") == method)
  l <- length(loglik)
  if(is.null(m.names)){
  	m.names <- 1:l
  }
  info <- get_info(loglik, k, n)[(l*x - (l-1)):(l*x)]
  weight <- exp(-0.5 * (info - min(info))) / sum(exp(-0.5 * (info - min(info))))
  results <- data.frame("model" = m.names, "loglik" = loglik, "k" = k, "n" = n, "info" = info, "info_wt" = weight)
  colnames(results)[5:6] <- c(method, paste(method, "_wt", sep = ""))
  rownames(results) <- NULL
  return(results)
}


#' Simulate a tree under Gause's law
#'
#' Function to simulate a tree under the assumption of competitive exclusion (Gause's law). 
#' @param b a numeric value between 0 and 1. The probability of speciation per extant lineage per generation.
#' @param n numeric. The desired number of species to simulate. 
#' @param t numeric. The desired number of generations to simulate. 
#' @param ext logical. If T, simulate a mass extinction event.  
#' @param ext_t numeric. The desired generation in which to simulate a mass extinction.
#' @param ext_t numeric. The desired generation in which to simulate a mass extinction.
#' @param ext_s numeric. The severity of the mass extinction between 0 and 1.
#' @return An object of class 'phylo'. 
#' @details This function simulates an evolutionary process in discrete time under the assumption of competitive exclusion (Gause's law). The progenitor species starts with an integer value = 0, which represents the species' niche. Each generation, a species can speciate (produce a daughter species). Speciation is assumed to be budding. Daughter species inherit the parents trait +/- 1. This assumes daughter species occupy a niche close to it's parent. Each generation, the order in which species attempt to speciate is randomised. If a species attempts to speciate, yet the daughter lineage would occupy an already occupied niche, speciation does not occur. If the niche of a species becomes occupied before the species can attempt to speciate, that species instead becomes exitinct. Users can optionally specify a mass extinction using the ext logical argument. The mass extinction will occur at the generation specied by the user using the ext_t argument. The number of lineages that become extinct at the mass extinction interval is equal to the number of extant species * the mass extinction severity (ext_s), rounded to the nearest integer. Lineages that become extinct are determined randomly. The sim_g_tree() function effectively simulates a process whereby the extinction rate is proportional to the number of extant species. Becuase of this, it produces tree topologies quite unlike those of birth-death trees. Trees produced by sim_g_tree() tend to be more asymmetrical (unbalanced) and have a higher proportion of extinct lineages. Diversification occurs linearly, rather than exponentially.  
#' @examples
#' ## Simulate a tree of approx. 200 taxa
#' t1 <- sim_g_tree(b = 0.2, n = 200, t = 100000, ext = F)
#' plot(ladderize(t1), show.tip.label = F, direction = "upwards", no.margin = T)
#' ## Simulate a tree of 200 generations, with no extinction.
#' set.seed(1)
#' t2 <- sim_g_tree(b = 0.2, n = 10000, t = 200, ext = F)
#' plot(ladderize(t2), show.tip.label = F, direction = "upwards", no.margin = T)
#' ## Simulate a tree of 200 generations, with a mass extinction at gen 100.
#' set.seed(1)
#' t3 <- sim_g_tree(b = 0.2, n = 10000, t = 200, ext = T, ext_t = 100, ext_s = 0.9)
#' plot(ladderize(t3), show.tip.label = F, direction = "upwards", no.margin = T)
#' 
#' @export 

sim_g_tree <- function(b, n = 100, t = 1000, ext = F, ext_t = NULL, ext_s = NULL){
	if(b <= 0 | b >1 ){
		stop("Birth rate (b) must be between 0 and 1")
	}
	if(ext == T){
		if(ext_s <= 0 | b >1 ){
			stop("Mass extinction severity (ext_s) must be between 0 and 1")
		}
	}
	#numeric vector to store trait info for taxa that were alive in the previous gen.
	trait_old <- c('t1' = 0)
	#character vector to store all species (living + extinct)
	species <- "t1"
	tree <- "t1x"
	#data frame to store birth and death info.
	bd_df <- data.frame(taxon = 't1', birth = 0, death = NA)
	#start evolution
	gen <- 0
	repeat{
		gen <- gen + 1
		#numeric vector to store trait info for taxa that are alive in the current gen.
		trait_new <- numeric()
		#Check gen == ext_t. If so, simulate a mass extinction of severity = ext_s. 
		if(ext == T){
			if(gen == ext_t){
				victims <- sample(names(trait_old), round(length(trait_old) * ext_s))
				trait_old <- trait_old[-which(names(trait_old) %in% victims)]
				for(i in victims){
					tree <- gsub(paste(i, "x", sep = ""), paste(i, ":", gen - bd_df[which(bd_df$taxon == i),2], sep = ""), tree)
					bd_df[which(bd_df$taxon == i),3] <- gen
				}
			}
		}
		#check if everything is extinct, if so break
		if(length(trait_old) == 0){
			break
		}
		#loop through each taxon that is alive.
		# note - bodge with as.character is so you can sample a vector of varying lengths including 1. 
		alive <- as.character(which(is.na(bd_df$death)))
		for(i in bd_df[as.numeric(sample(alive)), 1]){
			trait_x <- trait_old[[i]]
			names(trait_x) <- i
			#determine if species i goes extinct. 
			if(trait_x %in% trait_new == F){
				trait_new <- c(trait_new, trait_x)
				#determine if each alive taxon speciates (produces a daughter lineage). 
				if(runif(1) <= b){
					#new taxa inherits parents trait + or - 1   
					trait_y <- trait_x + sample(c(-1, 1), 1)
					if(trait_y %in% trait_new == F){
						n_tax <- paste("t", (length(species)+1), sep = "")
						species <- c(species, n_tax)
						names(trait_y) <- n_tax
						trait_new <- c(trait_new, trait_y)
						anc_blen <- gen - bd_df[which(bd_df$taxon == i),2]
						bd_df[which(bd_df$taxon == i),2] <- gen
						bd_df <- rbind(bd_df, data.frame(taxon = n_tax, birth = gen, death = NA))
						tree <- gsub(paste(i, "x", sep = ""), paste("(", i, "x, ", n_tax, "x):", anc_blen, sep = ""), tree)
					}
				}
			} else { 
				tree <- gsub(paste(i, "x", sep = ""), paste(i, ":", gen - bd_df[which(bd_df$taxon == i),2], sep = ""), tree)
				bd_df[which(bd_df$taxon == i),3] <- gen
			}
		}
		trait_old <- trait_new
		#break if species limit is reached or surpassed 
		if(length(species) >= n){
			break
		}
		#break if gen == t
		if(gen == t){
			break
		}
	}
	bd_df[which(is.na(bd_df[,3])),3] <- gen 
	for(i in names(trait_old)){
		x <- which(bd_df$taxon == i)
		tree <- gsub(paste(i, "x", sep = ""), paste(i, ":", bd_df[x, 3] - bd_df[x, 2], sep = ""), tree)
	}
	t1 <- read.tree(text = paste(tree, ";", sep = ""))
	t1$root.time <- gen
	return(t1)
}


#' Calculate tree to tree distances in parallel
#'
#' Function to calculate tree to tree distances in parallel using the foreach package. 
#' @param trees an object of class 'multiPhylo'.
#' @param method character specifying the distance metric. Can be "RF", "quartet", "CID", "MSID" or "SPR".
#' @param slices numeric. Number of sub matrices to compute
#' @param normalise logical. If T, distance metric is normalised.  
#' @return A square matrix of tree to tree distances. 
#' @details This function uses the foreach package to calculate tree to tree distances in parallel. Foreach allows for parallel computation of matrices, however tree to tree distance matrices are symmetrical, thus simply using foreach to calculate each cell of the matrix would result in 50 percent redundancy. To lower the amount of redundancy, the function splits the square matrix into a user specifed number of rectangular matrices (using the argument 'slices'). These are dealt with in serial, but each bvalue of the rectangular matrices is computed in parallel. For example, if length(trees) == 99, and slices == 3, the function will split the calculation into three submatrices x[1:33, 1:99], x[34:66, 34:99] and x[67:99, 67:99]. This lowers the redundancy from 50 percent to less than 25 percent. In theory, tree distance matrices can be computed without any redundancy, for example using the future.apply package. However, in practice, I have found this to be nowhere near as fast. Whether or not the dist_m() function is faster than computing the distances in serial depends on a number of factors including the number of trees, the number of tips, the complexity of the distance metric and the number of availible cores. In general, this function is recommended for computing distances between > 1000 trees. The ideal number of slices is difficult to determine. A greater number of slices decreases redundancy but increases the number of serial operations. In general, I recommend setting slices to 3. Unfortunately, Windows does not support forking and thus this function will only work for Linux or Mac users.
#' @examples
#' ## simulate two groups of trees
#' t1 <- rtree(40)
#' trees1 <- t1
#' for(i in 1:199){ 
#'   trees1 <- c(trees1, rSPR(t1, 1))
#' }
#' t2 <- rSPR(t1, 3)
#' t2 <- multi2di(t2)
#' trees2 <- t2
#' for(i in 1:199){ 
#'   trees2 <- c(trees2, rSPR(t2, 1))
#' }
#' trees <- c(trees1, trees2)
#' 
#' ## Register parallel backend for Windows with doParallel
#' cl <- makeCluster(2)
#' doParallel::registerDoParallel(cl)
#' 
#' ## Register parallel backend for Linux or MacOS with doParallel
#' nc <- parallel::detectCores()
#' doParallel::registerDoParallel(cores=nc)
#'
#' ## Calculate tree distances
#' tree_dists <- dist_m(trees, method = "RF", slices = 3, normalise = T)
#' ## Conduct PCA analysis on tree distances
#' pca_res <- prcomp(tree_dists)
#' ## Check proportion of variances of principal components
#' sum_res <- summary(pca_res)
#' # install.packages("ggplot2")
#' # install.packages("factoextra")
#' # library("ggplot2")
#' # library("factoextra")
#' fviz_eig(pca_res, addlabels=TRUE, hjust = -0.3)
#' ## plot PCA results
#' ggplot(pca_res$x, aes(PC1, PC2)) + geom_point(alpha = 0.4, size = 2) +  xlab(paste("PC1 = ", sum_res$importance[2,1], sep = "")) + ylab(paste("PC2 = ", sum_res$importance[2,2], sep = "")) + theme_minimal()
#' 
#' @export 

dist_m <- function(trees, method, slices = 3, normalise = F){
	if(normalise == T & method == "SPR"){
		warning("SPR distance cannot be normalised.")
	}
	col_length <- as.numeric(table(cut(rep(1:length(trees)), slices)))
	tree_sets <- list()
	trees2 <- trees
	for(i in 1:slices){
		tree_sets[[i]] <- head(trees2, col_length[[i]])
		tree_sets[[i]] <- head(trees2, col_length[[i]])
		trees2 <- tail(trees2, (length(trees2) - col_length[[i]]))
	}
	trees2 <- tail(trees, (length(trees) - 1))
	dist_list <- list()
	for(j in 1:slices){
		dist_list[[j]] <- foreach(x = iter(tree_sets[[j]]), .combine = 'cbind')%dopar% {
			x1 <- vapply( iter(trees2), function(y) {
				if(method == "RF"){
					return(RF.dist(x, y, normalize = normalise))
				}
				else if(method == "quartet"){
    				Qstat <- Quartet::QuartetStatus(x, y)
    				Qdist <- 2*Qstat[[4]]
    				if(normalise == T){
    					Qdist <- Qdist/Qstat[[1]]
    				}
    				return(Qdist)
				}
        		else if(method == "CID"){
    				return(TreeDist::ClusteringInfoDist(x, y, normalize = normalise))
				}
        		else if(method == "MSID"){
          			return(TreeDist::MatchingSplitInfoDistance(x, y, normalize = normalise))
        		}
				else if(method == "SPR"){
					return(SPR.dist(x, y))
				}
    	}, numeric(1) )
	}
	if(j < slices){
		trees2 <- tail(trees2, length(trees2) - col_length[[j]])
		}
	}
	for(k in 1:slices){
		dist_list[[k]] <- rbind(matrix(0, (length(trees) - nrow(dist_list[[k]])), ncol(dist_list[[k]])), dist_list[[k]])
	}
	dist_matrix <- dist_list[[1]]
	for(l in 2:slices){
		dist_matrix <- cbind(dist_matrix, dist_list[[l]])
	}

	mdim <- nrow(dist_matrix)
	dist_matrix <- foreach(m=1:mdim, .combine='rbind')%dopar%{
		x <- as.numeric(dist_matrix[m, ])
		x[c(m:mdim)] <- 0
		x
	}
	colnames(dist_matrix) <- NULL
	dist_matrix <- dist_matrix + t(dist_matrix)
	colnames(dist_matrix) <- NULL
	return(dist_matrix)
}


#' summary.nexdat
#'
#' Function to sumarise an object of class 'nexdat'. 
#' @param x an object of class 'nexdat'.
#' @return An object of class 'summary.nexdat'.
#' @details This function returns the proportion of missing characters and proportion of gaps per taxon.  
#' @examples
#' data(Lavoue2016) 
#' summary(Lavoue2016$dna)
#' 
#' @export 

summary.nexdat <- function(x){
    names(x) <- tolower(names(x))
    res <- list()
    if(any(class(x) == "nexdat") == F){
        warning("Object is not class nexdat")
    }

    # summarise combined data

    x <- lapply(x, tolower)
    comb_data <- x
    missing <- lapply(comb_data, function(z) which(z == "?"| z == "n" | z == "x"))
    gaps <- lapply(comb_data, function(z) which(z == "-"))
    NCHAR <- NA
    NTAX <- NA
    res[[1]] <- c("Ntax" = length(comb_data), "Nchar" = length(comb_data[[1]])) 
    res[[2]] <- as.matrix(data.frame("Ntax" = NTAX, "Nchar" = NCHAR))
    names(res) <- c("data", "partitions")
   
    prop_missing <- unlist(lapply(missing, function(z) round(length(z)/res[[1]][[2]], 2)))
    prop_gaps <- unlist(lapply(gaps, function(z) round(length(z)/res[[1]][[2]], 2)))
    res[[3]] <- as.matrix(data.frame("missing" = prop_missing, "gaps" = prop_gaps))
    names(res)[[3]] <- "total"
    class(res) <- "summary_nexdat"
    return(res)
}




#' summary.multi_nexdat
#'
#' Function to sumarise an object of class 'multi_nexdat'. 
#' @param x an object of class 'multi_nexdat'.
#' @return An object of class 'summary.nexdat'.
#' @details This function returns the proportion of missing characters and proportion of gaps per taxon. If partitions are specidied, the function also returns the proportion of missing data, gaps, ambiguous codings (e.g. polymorphic characters, characters that are coded as uncertain between two or more states) and unamibiguous characters per taxon per partition.  
#' @examples
#' data(Lavoue2016)
#' summary(Lavoue2016)
#' 
#' @export 

summary.multi_nexdat <- function(x){
    # check if object is class 'multi_nexdat'
    if(any(class(x) == "multi_nexdat") == F){
        warning("Object is not class multi_nexdat")
    }

    #list to store results
    res <- list()

    # check if partitions have names.
    if(is.null(names(x))) {
        stop("partition names not included. Partitions should be named 'standard', 'dna' or 'protein'!")
    } else {
        names(x) <- tolower(names(x))
    }

    # check if partitions are names anything other than 'standard', 'dna' or 'protein'.
    Pnames <- length(setdiff(names(x), c("standard", "dna", "protein")))
    if(Pnames > 0){
        warning(paste(Pnames, "partition name(s) not recognised. Partitions should be named 'standard', 'dna' or 'protein'."), sep = "")
    }

    # summarise combined data
    comb_data <- remove_part_info(x)
    missing <- lapply(comb_data, function(z) which(z == "?"| z == "n" | z == "x"))
    gaps <- lapply(comb_data, function(z) which(z == "-"))
    NCHAR <- unlist(lapply(x, function(z) length(z[[1]])))
    NTAX <- unlist(lapply(x, function(z) length(z)))

    res[[1]] <- c("Ntax" = length(comb_data), "Nchar" = length(comb_data[[1]])) 
    res[[2]] <- as.matrix(data.frame("Ntax" = NTAX, "Nchar" = NCHAR))
    names(res) <- c("data", "partitions")
   
    prop_missing <- unlist(lapply(missing, function(z) round(length(z)/res[[1]][[2]], 2)))
    prop_gaps <- unlist(lapply(gaps, function(z) round(length(z)/res[[1]][[2]], 2)))
    res[[3]] <- as.matrix(data.frame("missing" = prop_missing, "gaps" = prop_gaps))
    names(res)[[3]] <- "total"

    # summarise data partitions
    for(i in 1:length(x)){
        x[[i]] <- lapply(x[[i]], tolower)
        missing <- lapply(x[[i]], function(z) which(z == "?"| z == "n" | z == "x"))
        gaps <- lapply(x[[i]], function(z) which(z == "-"))
        if(names(x)[[i]] == "dna" | names(x)[[i]] == "rna"){
            unam_chars <- lapply(x[[i]], function(z) which(z == "a" | z == "c" | z == "g" | z == "t" | z == "u"))
            poly <- lapply(1:length(x[[i]]), function(z) (1:length(x[[i]][[z]]))[-c(missing[[z]], gaps[[z]], unam_chars[[z]])])
            names(poly) <- names(unam_chars)
        }
        if(names(x)[[i]] == "standard"){
            unam_chars <- lapply(x[[i]], function(z) which(is.na(suppressWarnings(as.numeric(z))) == F))
            poly <- lapply(1:length(x[[i]]), function(z) (1:length(x[[i]][[z]]))[-c(missing[[z]], gaps[[z]], unam_chars[[z]])])
            names(poly) <- names(unam_chars)
        }
        if(names(x)[[i]] == "protein"){
            poly <- lapply(1:length(x[[i]]), function(z) which(z == "b" | z == "z" | z == "j" | z == "u" ))
            unam_chars <- lapply(1:length(x[[i]]), function(z) (1:length(x[[i]][[z]]))[-c(missing[[z]], gaps[[z]], poly[[z]])])
        }
        part_Nchar <- length(x[[i]][[1]])
        prop_missing <- unlist(lapply(missing, function(z) round(length(z)/part_Nchar, 2)))
        prop_gaps <- unlist(lapply(gaps, function(z) round(length(z)/part_Nchar, 2)))
        prop_poly <- unlist(lapply(poly, function(z) round(length(z)/part_Nchar, 2)))
        prop_unam <- unlist(lapply(unam_chars, function(z) round(length(z)/part_Nchar, 2)))

        res[[i + 3]] <- as.matrix(data.frame(missing = prop_missing, gaps = prop_gaps, ambiguous = prop_poly, unambiguous = prop_unam))
        names(res)[[i + 3]] <- names(x)[[i]]
    }
    class(res) <- "summary_nexdat"
    return(res)
}


#' Leave-One-Out Cross-Validation
#'
#' Function to conduct leave-one-out cross-validation on an ancestral state estimation model. 
#' @param tree an object of class 'phylo'.
#' @param x character specifying the distance metric. Can be "RF", "quartet", "CID", "MSID" or "SPR".
#' @param model a character string containing the model or a transition model specified in the form of a matrix. See ace for more details.
#' @param fixedQ fixed value of transition matrix Q, if one is desired.
#' @param ... optional arguments, including pi, the prior distribution at the root node (defaults to pi="equal"). Other options for pi include pi="fitzjohn" (which implements the prior distribution of FitzJohn et al. 2009), pi="estimated" (which finds the stationary distribution of state frequencies and sets that as the prior), or an arbitrary prior distribution specified by the user. 
#' @return A value between 0 and 1 representing the mean raw error of the model.
#' @details to do!
#' @examples
#' ## Load data
#' head(vert_data$morph)
## Convert secondary osteon character to tip priors, interpreting inapplicable ('-') as an additional state (i.e. absent).
#' tp <- get_tip_priors(vert_data$morph[,2], extra_state = T)
#' colnames(tp[[1]]) <- c("bone absent", "secondary osteons absent", "secondary osteons present")
## Load parallel packages
#' library(future.apply)
#' plan("multisession")
## Compare mean error of ER and ARD models using loo_cv.
#' loo_cv(tree = vert_data, x = tp[[1]], model = "ER")
#' loo_cv(tree = vert_data, x = tp[[1]], model = "ARD")

#' @export 

loo_cv <- function(tree, x, model, fixedQ=NULL, ...){
    args.x <- list(...)
    if (is.matrix(x)) {
        x <- x[tree$tip.label, ]
        m <- ncol(x)
        states <- colnames(x)
    } else {
        x <- to.matrix(x, sort(unique(x)))
        x <- x[tree$tip.label, ]
        m <- ncol(x)
        states <- colnames(x)
    }
    res <- future.apply::future_lapply(1:nrow(x), function(i){
        true_state <- which(x[i,] > 0)
        x2 <- x
        x2[i,] <- 1/ncol(x2)
        object <- fitMk(tree, x2, model, fixedQ=NULL, args.x)
        cv <- ancr(object, tips=TRUE)
        sum(cv$ace[i,true_state])
    }, future.seed=TRUE)
    return(1-(mean(unlist(res))))
}



