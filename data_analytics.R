library(data.table)
library(dplyr)
library(edgeR)
library(GO.db)
library(org.Mm.eg.db)
library(doParallel)
library(ggplot2)
library(vsn)
library(DESeq2)
library(tidyr)
library(igraph)
library(doParallel)
library(statmod)
library(factoextra)
library(gridExtra)


#' plot plot_boxplot
#' @param gene annotation data
plot_boxplot <- function(sdata, metric = 'logFC', xf = 6, lf = 1, pathname){
  print(pathname)
  
  tiff(paste(pathname,sdata$class, '_', 'genes.tiff', sep = ''),
       width = 1280, height = 780,
       res = 180)
  p = ggplot(sdata , aes(x = time, y = sdata[[metric]],
                      colour = group, fill=group)) +
    geom_boxplot(outlier.colour="black", outlier.shape = 8,
                 outlier.size=1, notch=FALSE) +
    geom_smooth(method = "loess", se=FALSE, aes(group=2), col = 'purple',
                size = 0.9) + 
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = xf),
          axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 6),
          legend.text=element_text(size=4),
          legend.key.size = unit(lf, 'line'),
          legend.position = "none") +
    labs(y = metric) +
    geom_hline(yintercept = c(1, -1), linetype="dashed",
               color = "black", size = 0.2) +
    facet_grid(key ~ group)
  print(p)
  dev.off()
  
}


#' plot_heatmap
#' @param gene annotation data
plot_heatmap <- function(sdata, metric = 'logFC', gname, pathname){
  
  for(g in unique(sdata[[gname]])){
    
    tiff(paste(pathname, sdata$class, '_', gsub(' ', '', g), '.tiff', sep =''), 
         width = 680, height = 780,
         res = 180)
    ss <- subset(sdata, group == g)
    rownames(ss) = 1:nrow(ss)
    
    mine.heatmap <- ggplot(
      data = ss, mapping = aes(x = time,  y = gene_id,
                               fill = logFC)) + geom_tile(height=0.9) +
      xlab(label = "time") + ylab(label = paste(g, "genes")) +
      theme(panel.background = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, size = 6), 
            axis.text.y = element_text(size = 4),
            strip.text = element_text(size = 6)) +
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank()) + 
      scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred") + 
      facet_grid(annotation ~ key,  drop = TRUE, scales = "free", 
                 space = 'free')
    print(mine.heatmap)
    
    dev.off()
    #if(readline(g) == 'q') { break }
    
  }
  
}



p2 <- function(subn_hor, n, metric = 'logFC'){
  
  meta_name <- c("gene_id", "class", "group", "subgroup", 'annotation',
                 "Human-Readable-Description")

  i <- match(c(meta_name), colnames(subn_hor))
  DT.m1 <- melt(subn_hor, id.vars = meta_name,
               measure.vars = names(subn_hor)[-i], variable.name = 'treatment',
               value.name = metric)
  
  dd <- separate(DT.m1, treatment , c("key", "time"), "_", 
                 extra = "merge", remove = FALSE)
  dd$time <- as.factor(dd$time)  
  dd$time <- ordered(dd$time,levels=c("1h", "5h", "24h", "48h", 
                                             "72h",  "120h", "168h"))
  dd$key <- ordered(dd$key,levels=c('S', 'T', 'B'))
  
  #i = which(colnames(dd) == 'subgroup')
  #colnames(dd)[i] = 'annotation'
  return(dd)
 
}


plot_headmap <- function(){
  
  heatmap.2(as.matrix(hormone_genes[1:71,5:15]), 
            na.color = "white", 
            scale="none", 
            Rowv = FALSE, Colv=FALSE,
            col=redgreen)

}

#sdata = subn; gclass = 'hormone'; cl = 3; metric ='FDR';

side_plot <- function(sdata, gclass, cl, metric){
  
  sub_sdata <- subset(sdata, geneclass == gclass)
  # Model Based Clustering
  km.res <- kmeans(sub_sdata[, -match(c("y","Human-Readable-Description",
                                      'geneclass'), colnames(sub_sdata))], 
                   3, nstart = 25)
  o = fviz_cluster(km.res, data = sub_sdata[, -match(c("y","Human-Readable-Description",
                        'geneclass'), colnames(sub_sdata))], 
                   frame.type = "convex", main = gclass,
                   pointsize = 0.5, labelsize = 6, 
                   outlier.color = 'black', outlier.shape = 19)
  
  # Add group cluster
  sub_sdata <- mutate(sub_sdata, cluster = as.numeric(km.res$cluster))
  
  #s <- aggregate(sub_hor, by = list(sub_hor$group), FUN=mean)
  n <- c("y", "Human-Readable-Description", "geneclass",'cluster')
  i = match(n, colnames(sub_sdata))
  DT.m1 = melt(sub_sdata, id.vars = n,
               measure.vars = names(sub_sdata)[-i], variable.name = 'treatment',
               value.name = metric)
  DT.m1[['cluster']] <- as.factor(DT.m1[['cluster']])
  
  p = ggplot(DT.m1 , aes(x = treatment, y = DT.m1[[metric]], 
                         colour = cluster, fill=cluster)) + 
    geom_boxplot(outlier.colour="black", outlier.shape=8,
                 outlier.size=1, notch=TRUE)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = metric) +
    geom_hline(yintercept = c(1, -1), linetype="dashed",
               color = "black", size = 0.2)
  p
  
  grid.arrange(o, p, nrow=2)
  
  return(sub_sdata)
}




#' get_score_data FDR, logFC
#' @param sdata list with gene expression results
#' @param metric score to extract from DB

get_score_data <- function(sdata, metric){
  
  
  m <- data_res[[names(data_res)[1]]][, c('y', 'Human-Readable-Description',
                                          metric)]
  colnames(m)[which(colnames(m) == metric)] <- names(data_res)[1]
  
  for(exp_name in names(data_res)[2: length(names(data_res))]){
    if(data_res[[exp_name]] != 0){
      m <- merge(m, data_res[[exp_name]][, c('y', 
                                             metric)], by = 'y' )
      colnames(m)[which(colnames(m) == metric)] <- exp_name
    }
  }
  
  return(m)
  
}

#' create graph
#' @param exp_mat expression matrix

plot_graph <- function(exp_mat){
  
  #Create a graph adjacency based on correlation distances between 
  #genes in  pairwise fashion.
  g <- graph.adjacency(
    as.matrix(as.dist(cor(t(exp_mat)))),
    mode="undirected",
    weighted=TRUE,
    diag=FALSE
  )
  
  #Simplfy the adjacency object
  #g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
  
  #Colour negative correlation edges as blue
  E(g)[which(E(g)$weight<0)]$color <- "darkblue"
  
  #Colour positive correlation edges as red
  E(g)[which(E(g)$weight>0)]$color <- "darkred"
  
  #Convert edge weights to absolute values
  E(g)$weight <- abs(E(g)$weight)
  
  #Change arrow size
  #For directed graphs only
  #E(g)$arrow.size <- 1.0
  
  #Remove edges below absolute Pearson correlation 0.8
  g <- delete_edges(g, E(g)[which(E(g)$weight < 0.8)])
  
  #Assign names to the graph vertices (optional)
  V(g)$name <- V(g)$name
  
  #Change shape of graph vertices
  V(g)$shape <- "circle"
  
  #Change colour of graph vertices
  V(g)$color <- "lightgreen"
  
  #Change colour of vertex frames
  V(g)$vertex.frame.color <- "white"
  
  #Scale the size of the vertices to be proportional to the level of expression
  # of each gene represented by each vertex
  #Multiply scaled vales by a factor of 10
  scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
  vSizes <- (scale01(apply(exp_mat, 1, mean)) + 1.0) * 10
  
  #Amplify or decrease the width of the edges
  edgeweights <- E(g)$weight * 2.0
  
  #Convert the graph adjacency object into a minimum spanning tree based 
  #on Prim's algorithm
  mst <- mst(g, algorithm="prim")
  
  #Plot the tree object
  plot(
    mst,
    layout=layout.fruchterman.reingold,
    edge.curved=TRUE,
    vertex.size=vSizes,
    vertex.label.dist=-0.5,
    vertex.label.color="black",
    asp=FALSE,
    vertex.label.cex=0.8,
    edge.width=edgeweights,
    edge.arrow.mode=0,
    main="Claviceps "
  )
  return(g)
}  


#' Load claviceps data
load_clavi_data_frame <- function(){
  
    pathname <- 'claviceps/'
    files <- list.files(path = pathname, pattern = '.csv')
    
    clavi <- read.csv(paste(pathname, files[1], sep ='\\/'), header = 
                        TRUE, sep = ',')[,1:3]
    W <- read.csv(paste(pathname, files[7], sep ='\\/'), header = 
                    TRUE, sep = ',')[, 1:4]
    
    g <- merge(clavi, W, by = 'Name')
    
    for(fname in files[7:length(files)]){
      
      g1 <- read.csv(paste(pathname, fname, sep ='\\/'), header = 
                      TRUE, sep = ',')[, c(1, 5:7)]
      g <- merge(g, g1, by = 'Name')
      
    }
    
    rownames(g) <- g[['Name']]
    g <- g[, -1]
    
    return(g)

}

#' Load wheat data
load_claviwheat_data_frame <- function(){
  
  pathname <- 'wheat'
  group_order <- list()
  files <- list.files(path = pathname, pattern = '.csv')
 
  g <- read.csv(paste(pathname, files[1], sep ='\\/'), header = 
                   TRUE, sep = ',')
  group_order[[strsplit(files[1], '_')[[1]][1]]] <- colnames(g)[2:ncol(g)]
 
  for(fname in files[2:length(files)]){
    
    g1 <- read.csv(paste(pathname, fname, sep ='\\/'), header = 
                     TRUE, sep = ',')
    group_order[[strsplit(fname, '_')[[1]][1]]] <- colnames(g1)[2:ncol(g1)]
    g <- merge(g, g1, by = 'Name')
    
  }
  
  rm(g1)
  rownames(g) <- g[['Name']]
  g <- g[, -1]
  g <- setDT(g, keep.rownames = TRUE)
  return(list(data = g, groups = group_order))
  
}

#' @param sdata dataset
plot_fre <- function(sdata){
  
  for(cname in colnames(sdata)[2:ncol(sdata)]){
    par(mfrow = c(1,2))
    hist(sdata[[cname]][ii], main = paste('Claviceps', cname))
    mtext(side = 3, mean(sdata[[cname]][ii]))
    hist(sdata[[cname]][i], main = paste('Wheat', cname))
    mtext(side = 3, mean(sdata[[cname]][i]))
    if(readline(cname) == 'q') { break}
  }
}


look_all <- function(){

  design <- model.matrix(~ treat + time + tissue + 
                         treat:time:tissue, data=sub_exp_settings)
  fit <- glmQLFit(y, design)
  
}

#' @param y reads data
#' @param fot fitted model
#' @param con_vec constrast vector
#' @return results

perform_analytics_fast <- function(g, y, fit, cont_vec, annot_file, 
                                   sub_exp_settings, logcpm){

  qlf <- glmQLFTest(fit, contrast = cont_vec[['vector']])
  print(topTags(qlf, p.value = 0.01))
  #summary(decideTests(qlf))
  plotMD(qlf)
  abline(h=c(-1, 1), col="blue")
  
  
  # Search for expression data for cont_vec experiment in exp settings
  i = match(sub_exp_settings$sample[union(which(sub_exp_settings$p ==
                                                   cont_vec$s1), 
                                              which(sub_exp_settings$p == 
                                                    cont_vec$s2) 
                                                   )], 
        colnames(logcpm))
  
  all_tags <- merge(qlf$table, logcpm[,i], by = 'row.names')
  all_tags <- merge(all_tags, g$rn, by.x = 'Row.names', by.y = 'row.names')
  all_tags <- mutate(all_tags, FDR = p.adjust(PValue, method = 'BH', 
                                              n = nrow(all_tags)))
  
  all_tags <- merge(all_tags, annot_file[, 1:2], by.x = 'y', by.y = 'Gene-ID')
  rownames(all_tags) <- all_tags[['y']]
  m <- c('Row.names')
  all_tags <- all_tags[, -match(m, colnames(all_tags))]
  
  return(all_tags)
  
}


#' @param cnames samples from fit
#' @param math_word_a param1
#' @param math_word_b contrast
search_match <- function(cnames, match_word_a, match_word_b){
  
  
  mat_cont <- rep(0, length(cnames))
  mat_cont[grep(match_word_a, cnames)] <- 1
  mat_cont[grep(match_word_b, cnames)] <- -1
  
  return(list('vector' = mat_cont, 's1' = match_word_a, 's2' = match_word_b))
  
  
}

#' @param fit fitted model

create_con_mat <- function(fit, tp, tr){
  
  con_man_list <- list()
  # # by time point
  # for(t in tp){
  #   
  #   con_man_list[[paste('all_T_allt_', t, 'h', sep='')]] <- 
  #     search_match(colnames(fit), paste('cp_', t, sep=''), 
  #                  paste('mock_', t, sep=''))
  # }
  
  # by timepoint by treatment
  
  for(t in tp){
    for(treat in tr){
      con_man_list[[paste(treat, '_', t, 'h', sep='')]] <-
        search_match(colnames(fit), paste('cp_', t, '_', treat, sep=''), 
                     paste('mock_', t, '_', treat, sep=''))
      
    }
  }
  
  # # by treatment
  # for(treat in tr){
  #   
  #   con_man_list[[paste('all_T_', 't', treat, '_', 'all_h', sep='')]] <-
  #     search_match(colnames(fit), 
  #                  paste(paste('cp_', '\\d{1}', '_', treat, sep=''), '|', 
  #                        paste('cp_', '\\d{2}', '_', treat, sep=''), sep=''), 
  #                  paste(paste('mock_', '\\d{1}', '_', treat, sep=''), '|', 
  #                        paste('mock_', '\\d{2}', '_', treat, sep=''), sep='')
  #                  )
  #   
  # }
  
  return(con_man_list)
  
}



#' Create network
#' @param expression data
#' @param thr cut-off threshold
#' @return updated_g
populate_network <- function(sdata, thr = 0.01){
  
    
      ig <- plot_graph(sdata[, -match(c('y', 'logFC', 'logCPM', 'F', 'PValue', 
                                            'FDR', 'Human-Readable-Description'), 
                                          colnames(sdata))])
      df <- igraph::as_data_frame(ig, 'both')
      df$vertices <- df$vertices %>% 
        left_join(sdata[, match(c('y','logFC', 'logCPM', 'F', 'PValue', 'FDR', 
                                      'Human-Readable-Description'),
                                    colnames(sdata))], c('name'='y'))
      
      updated_g <- graph_from_data_frame(df$edges,
                                         directed = F,
                                         vertices = df$vertices)
    
    
    return(updated_g)
    
}


#' Perform expression analysis
#' @param sdata raw dataset
#' @param annot_file annotation file

perform_analytics <- function(sdata, annot_file, design, filop){
  
  
  y <- DGEList(counts = sub_g[, -1, with=FALSE])
  
  # Remove counts according to threshold
  if(filop == 1){
    keep <- rowSums(cpm(y) > 1) >= 2
  } else{
    keep <- filterByExpr(y)
  }
    # 2 keep <- rowSums(cpm(y) > 1) >= 3
  #
  #keep <- rowSums(cpm(y) > 6)
  y <- y[keep, keep.lib.sizes = FALSE]
  
  # Normalise data
  y <- calcNormFactors(y)
  
  # Estimate dispersion
  y <- estimateDisp(y, design = design, robust=TRUE)
  save(y, file = 'y.dat')
  # Perform counts
  #load('y.dat')
  # Perform quasi-likelihood F-tests
  fit <- glmQLFit(y, design)
  save(fit, file ='fit.dat')
  # qlf <- glmQLFTest(fit, contrast = c(rep(-1, 14), rep(1, 14)))
  # #summary(decideTests(qlf))
  # #plotMD(qlf)
  # #abline(h=c(-1, 1), col="blue")
  # 
  # logcpm <- cpm(y, prior.count = 2, log=TRUE)
  # all_tags <- merge(qlf$table, logcpm, by = 'row.names')
  # all_tags <- merge(all_tags, g[, 1:2], by.x = 'Row.names', by.y = 'row.names')
  # all_tags <- mutate(all_tags, FDR = p.adjust(PValue, method = 'BH', 
  #                                             n = nrow(f)))
  # all_tags <- merge(all_tags, annot_file[, 1:2], by.x = 'rn', by.y = 'Gene-ID')
  # rownames(all_tags) <- all_tags$rn
  # m <- c('Row.names', 'EB.y')
  # all_tags <- all_tags[, -match(m, colnames(all_tags))]
  
  return(list(fit = fit, y = y))
  
}

#' create graphs
#' @param data_res list containing expression data
#' @param gene_group target genes
create_graphs <- function(data_res, gene_group, name_group){
  
  net_res <- net_res_cg <- list()
  foreach(exp_name = names(data_res)) %do% {
    print(exp_name)
    subset_genes <- subset(data_res[[exp_name]], FDR <= 0.03)
    if(dim(subset_genes)[1] > 2){
      net_res[[exp_name]] <- populate_network(subset_genes)
      write_graph(net_res[[exp_name]], 
                  file = paste(exp_name, '_igraph.xml', sep =''),  "graphml")
      gp = subset_genes[subset_genes[['y']]
                                %in% unique(gene_group), ]
      if(!is.null(gp)){
        if(dim(gp)[1] > 2){
        net_res_cg[[exp_name]] <- populate_network(gp)
        write_graph(net_res_cg[[exp_name]], 
                    file = paste(exp_name, '_', name_group, '_ecan_igraph.xml', sep =''),  
                    "graphml")
        }
      }
    } else {
      net_res[[exp_name]] == 0
    }
  }
}



  can_genes_eleni <- read.csv('gene_list.csv', header = TRUE)
  hormone_genes <- read.csv('hormone_heatmap_annotation.csv', header = TRUE)
  defence_genes <- read.csv('defence_heatmap_annotation.csv', header = TRUE)
  annot_file <- fread('iwgsc_refseq_wheat_ALL_HC_LC.txt', header=TRUE,
                      sep = '\t')
  exp_settings <- read.csv('exp_settings.csv', header = TRUE)
  
  test <- 'treat'
  load('wheat.dat')
  #load('fit.dat')
  #load('y.dat')
  
  
  i <- grep('Traes', wheat$data$rn)
  ii <- setdiff(1:nrow(wheat$data), i)
  g <- wheat$data[i, ]
  sub_exp_settings <- filter(exp_settings, time > 0)
  sub_exp_settings <- mutate(sub_exp_settings, p = paste(treat, time, 
                                                         tissue, sep='_'))
  sub_exp_settings$p <- as.factor(sub_exp_settings$p)
  
  sub_g <- g[, c(1, match(sub_exp_settings$sample, colnames(g))), with = FALSE]
  design <- model.matrix(~ 0 + sub_exp_settings$p)
  colnames(design) <- levels(sub_exp_settings$p)
  
  
  pathname <- c('norm_filter/', 'str_filter/')
  
  #Calculate fit
  for(ftype in 1:1){
    p <- perform_analytics(sub_g, annot_file, design, ftype)
    fit <- p$fit
    y <- p$y
    
    
    # Normalise data
    # Compute counts per million (CPM) or reads per kilobase per million (RPKM)
    logcpm <- cpm(y, prior.count = 2, log=TRUE)
    
    
    #Create design matrix
    con_mat <- create_con_mat(fit, unique(sub_exp_settings$time), 
                               unique(sub_exp_settings$tissue))
    
     data_res <- net_res <- net_res_cg <- list()
  
     cl <- makeCluster(detectCores() -1)
     registerDoParallel(cl)
    
     pdf(paste(pathname[ftype], '_FCplot.pdf', sep =''))
     
     foreach(exp_name = names(con_mat)) %do% {
        g = g[,1:2]; cont_vec = con_mat[[exp_name]]
        if(sum(con_mat[[exp_name]][['vector']] == 0) ==
           length(con_mat[[exp_name]][['vector']])){
    
          data_res[[exp_name]] = 0
        } else{
    
          data_res[[exp_name]] <- perform_analytics_fast(g[,1:2], y, fit,
                                     con_mat[[exp_name]],
                                   annot_file, sub_exp_settings, logcpm)
        }
     }
    
    dev.off()
    
    #save(data_res, file = 'data_res.dat')
    
   #Check that columns match
   sapply(con_mat, function(x) colnames(fit)[which(x[['vector']] != 0)])
    
   
   #load('data_res.dat')
  
   n <- get_score_data(data_res, 'logFC')
   hormone_genes <- hormone_genes[, c("class", "group", "subgroup", "annotation", 
                                       "gene_id")]
   defence_genes <- defence_genes[, c("class", "group", "subgroup", "annotation", 
                                      "gene_id")]
   
   gene_table <- rbind(hormone_genes, defence_genes)
   gene_table <- merge(gene_table, n, by.x = 'gene_id', by.y = 'y')
   
   
   
   # hormone
   
  
   sdata_el <- p2(subset(gene_table[, -match('S_120h', colnames(gene_table))], 
                         class == 'hormone'), n, 'logFC')
   plot_boxplot(sdata_el, 'logFC', 6,1, pathname[ftype])
   plot_heatmap(sdata_el, 'logFC', 'group', pathname[ftype])
  
   # defence 
   sdata_el <- p2(subset(gene_table[, -match('S_120h', colnames(gene_table))], 
                         class == 'defence'), n, 'logFC')
   plot_boxplot(sdata_el, 'logFC', 2, 0.5, pathname[ftype] )
   plot_heatmap(sdata_el, 'logFC', 'group', pathname[ftype] )
      

  }
  
 #    break
 # 
 # b <- c('chitinase', 'auxin', 'Rg4', 'glucanase')
 # 
 # 
 # #pattern <- paste0('^\\w*', b[2], '\\w*$', collapse = '|')
 # 
 # pdf('p.pdf')
 # for(tname in b){
 #   i <- grep(tname, n$`Human-Readable-Description`)
 #      if(length(i) > 0){
 #        m <- separate(melt(n[i,]), variable, into=c('tissue', 'time'), '_')
 #        #m[['time']] <- as.numeric(m[['time']])
 #        m[['time']] <- ordered(m[['time']],levels=c("1h", "5h", "24h", "48h", 
 #                                            "72h",  "120h", "168h"))
 #        p <- ggplot(m, aes(x = time, y = value)) + 
 #          geom_boxplot(outlier.colour="black", outlier.shape=8,
 #                       outlier.size=1, notch=FALSE) +
 #          facet_grid(tissue ~.) +
 #          geom_hline(yintercept = c(2, -2), linetype="dashed",
 #                     color = "red", size = 0.2) +
 #          theme(panel.background = element_blank()) +
 #          ggtitle(tname)
 #       
 #        print(p)
 #      }
 #      
 #  }
 # 
 # dev.off()
 # 
 # 
 # 
 # o = unlist(sapply(n$`Human-Readable-Description`[grep('monooxygenase', 
 #                                            n$`Human-Readable-Description`)], 
 #        function(x) grep('P450', x)))
 # 
 # match(names(o), n$`Human-Readable-Description`)
