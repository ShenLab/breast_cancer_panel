pedigree <- function(){
# path <- "/home/local/ARCS/qh2159/breast_cancer/pedigree/Family_Pedigree"
# localpath:
path="/Volumes/TwoT/Desktop/breastcancer/PedigreeFiles/150504_Family_Pedigree/150504 Family Pedigree"
files <- list.files(path,pattern=".svg$",full.name=TRUE)

allped <- c()
for(i in 1:length(files)){
    onef <- files[i]
    onet <- onepedigree(onef)
    fam <- gsub(".svg$","",basename(onef))
    onefam <- cbind(fam,onet)
    allped <- rbind(allped,onefam)
}

colnames(allped) <- c("Family ID","Paternal ID","Maternal ID","Individual ID")
write.csv(allped[,c(1,4,2,3)],file="ALL_pedigree.csv",row.names=FALSE)

}

onepedigree <- function(onef){

con <- file(onef,"r")
onetet <- readLines(con,warn=FALSE)
close(con)

famlist <- c()
idsex <- c()
allmem <- c()

for(i in 1:length(onetet)){
    
   if(grepl("<rect class=\"solid\" id=",onetet[i])){
       tmp <- unlist(strsplit(onetet[i]," "))
       maleid <- gsub("\\D","",tmp[which(grepl("^id",tmp))])
       a1 <- c(maleid,"Male")
       idsex <- rbind(idsex,a1)
       
       y1 <- gsub("y=","",tmp[which(grepl("^y=",tmp))])
       y1 <- gsub("\"","",y1)
       allmem <- rbind(allmem,c(a1,y1))
   }
   
   if(grepl("<circle class=\"solid\" id=",onetet[i])){
       tmp <- unlist(strsplit(onetet[i]," "))
       femaleid <- gsub("\\D","",tmp[which(grepl("^id",tmp))])
       a2 <- c(femaleid,"Female")
       idsex <- rbind(idsex,a2)

       y2 <- gsub("cy=","",tmp[which(grepl("^cy=",tmp))])
       y2 <- gsub("\"","",y2)
       allmem <- rbind(allmem,c(a2,y2))
   }
   
   if(grepl("class=\"mating\"",onetet[i])){
       tmp <- unlist(strsplit(onetet[i]," "))
       parid <- gsub("id=","",tmp[which(grepl("^id",tmp))])
       parid <- gsub("\"","",parid)
       ids <- unlist(strsplit(parid,":"))
       
       pars <- 1:2
       if(all(ids %in% idsex[,1])){
           pars[1] <- intersect(ids,idsex[idsex[,2]=="Male",1])
           pars[2] <- intersect(ids,idsex[idsex[,2]=="Female",1])
       }else if(sum(ids %in% idsex[,1])==1){
           if(idsex[idsex[,1] %in% ids,2]=="Male"){
               pars[1] <- intersect(ids,idsex[idsex[,2]=="Male",1])
               pars[2] <- setdiff(ids,idsex[,1])
           }else{
               pars[1] <- setdiff(ids,idsex[,1])
               pars[2] <- intersect(ids,idsex[idsex[,2]=="Female",1])
           }
       }else{ print("Error!!!");}
       
       y <- gsub("y1=","",tmp[which(grepl("^y1=",tmp))])
       y <- gsub("\"","",y)
       
       onetri <- c(pars,y)
       famlist <- rbind(famlist,onetri)
       allmem <- rbind(allmem,onetri)
   }
   
}

trios <- c()
allV <- paste(allmem[,1],allmem[,2],allmem[,3],sep="_")
fams <- paste(famlist[,1],famlist[,2],famlist[,3],sep="_")
npar <- dim(famlist)[1]

for(i in 1:npar){
    sub <- which(allV==fams[i])
    
    if(i < npar){ 
        j <- which(allV==fams[i+1])-1;
        for(k in (sub+1):j){
            if(length(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3]))) > 0){
                psub <- max(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3])))
                onetrio <- c(famlist[psub,1:2],allmem[k,1])
                if(!(allmem[k,1] %in% famlist[i+1,1:2])){
                    trios <- rbind(trios,onetrio)
                }else if(!(setdiff(famlist[i+1,1:2],allmem[k,1]) %in% trios[trios[,1]==famlist[psub,1] & trios[,2]==famlist[psub,2],3])){
                    trios <- rbind(trios,onetrio)
                }
            }
        }
    }
    
    if(i == npar){ 
        j <- length(allV);
        for(k in (sub+1):j){
            psub <- max(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3])))
            onetrio <- c(famlist[psub,1:2],allmem[k,1])
            trios <- rbind(trios,onetrio)
        }
    }

}


trios
}

phenomap <- function(){
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    idmap <- pheno[,2:3]
    
    allped <- read.csv("ALL_pedigree.csv")
    allpedmap <- cbind(allped,"","","")
    allpedmap[,5] <- idmap[match(allped[,2],idmap[,1]),2]
    allpedmap[,6] <- idmap[match(allped[,3],idmap[,1]),2]
    allpedmap[,7] <- idmap[match(allped[,4],idmap[,1]),2]
    
    allpedmap[is.na(allpedmap)] <- ""
    colnames(allpedmap) <- c("Family ID","Individual ID","Paternal ID","Maternal ID","Individual Subject ID","Paternal Subject ID","Maternal Subject ID")
    write.csv(allpedmap,file="ALL_pedigree_map.csv",row.names=FALSE)
}

