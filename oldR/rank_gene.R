rank_gene <- function(gene, net, alpha){

    # column normalized
    W <- net$matrix
    W <- W/matrix(colSums(W),net$size,net$size,byrow=TRUE)
    
    # run the RWR algorithm
    p0 <- matrix(0,net$size,1,dimnames=list(net$node,1))
    tmpg <- intersect(gene,net$node)
    p0[match(tmpg,net$node)] <-  1
    
    eng <- sum(p0)
    p0 <- p0/eng
    d <- p0
    p1 <- RWR(W,p0,alpha,d)

    p1
    
}

testtmp <- function(){
    netfile <- "STRINGnetmap.txt"
    net <- STRINGnetwork(netfile)
    net$matrix[net$matrix>0] <- 1
    
    
    tcgat <- read.csv("../genelist/TCGA/SMGg.csv",skip=4)
    tcgag <- unique(tcgat[,1])
    a <- rank_gene(tcgag,net,0.5) 
    output_rank(net$node,a,wf="TCGA_rank.txt")

    
    gene1 <- unlist(read.table("../genelist/Genelist1.txt"))
    a1 <- rank_gene(gene1,net,0.5)
    output_rank(net$node,a1,wf="Panel1_rank.txt")
    
    
}

output_rank <- function(gene,p1,wf){

    tmp <- cbind(gene,p1)
    tmp <- tmp[order(-p1),]
    write.table(tmp,file=wf,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
}

RWR <- function(W,p0,alpha,d){
    # RWR: random walk with restart 
    
    # parameters:
    # W: the edge weights matrix
    # p0 : the start point for iteration
    # alpha: the restart probability
    # d: starting vector
    
    p1 <- (1-alpha)*W%*%p0 + alpha*d;
    iter <- 1
    while (sum(abs(p0-p1))>10^(-10)){
        p0=p1;
        p1 <- (1-alpha)*W%*%p0 + alpha*d;
        iter <- iter + 1
        #print(iter)
    }
    
    #print(iter)
    p1
    
}

STRINGnetwork <- function(file){
    
    
    #file <- "networks/STRING.interaction600.txt"
    links.text <- as.matrix(read.table(file))
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    #links.text[,1] <- gsub("9606.","",links.text[,1])
    #links.text[,2] <- gsub("9606.","",links.text[,2])
    links.text[,1] <- mapT[match(links.text[,1],mapT[,1]),2]
    links.text[,2] <- mapT[match(links.text[,2],mapT[,1]),2]
    subs <- !is.na(links.text[,1]) & !is.na(links.text[,2])
    links.text <- links.text[subs,]
    
    tmp <- links.text[,c(2,1,3)]
    colnames(tmp) <- c("V1","V2","V3")
    net.text <- rbind(links.text,tmp)
    
    net <- read_net(net.text)
    
    net
}

read_net <- function(net.text){
    
    net.node <- unique(union(net.text[,1],net.text[,2]))
    net.node <- net.node[net.node != ""]
    net.size <- length(net.node)
    net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,2]))
    net.edge <- net.edge[net.edge[,2] != "", ]
    net.edge <- net.edge[net.edge[,1] != "", ]
    net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
    net.matrix[net.edge] <- as.numeric(net.text[,3])
    list(size=net.size, node=net.node, matrix=net.matrix)
    
}

top_gener <- function(){
    source("~/.Rprofile")
    file <- "data/TCGA_ToppGeneData_Interaction.csv"
    #file <- "PanelToppGeneData.csv"
    tcga <- read.csv(file)
    tcgar <- tcga[,c("GeneSymbol","Average.Score")]
    tcgar <- tcgar[order(-tcgar[,"Average.Score"]),]
    qwt(tcgar,file="data/TCGA57g.txt")
    #qwt(tcgar,file="data/Panel29g.txt")
   

}