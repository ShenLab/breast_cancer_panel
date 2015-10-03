#pedigree <- function(){}

path <- "/home/local/ARCS/qh2159/breast_cancer/pedigree/Family_Pedigree"
onef <- "/home/local/ARCS/qh2159/breast_cancer/pedigree/Family_Pedigree/100018.svg"
onef <- "/Volumes/TwoT/Desktop/breastcancer/PedigreeFiles/150504_Family_Pedigree/150504 Family Pedigree/100018.svg"
con <- file(onef,"r")
onetet <- readLines(con,warn=FALSE)
close(con)

famlist <- c()
idsex <- c()
u <- 1
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
       u <- u+1
   }
}

trios <- c()
allV <- paste(allmem[,1],allmem[,2],allmem[,3],sep="_")
fams <- paste(famlist[,1],famlist[,2],famlist[,3],sep="_")
npar <- dim(famlist)[1]

for(i in 1:npar){
    sub <- which(allV==fams[i])
    if(i < npar){
        j <- which(allV==fams[i+1])-1
        if(j == sub+1){
            onetrio <- c(allmem[sub,1:2],allmem[j,1])
            trios <- rbind(trios,onetrio)
        }else if(j >= sub+2){
            for(k in (sub+1):(j-1)){
                if(as.numeric(allmem[k,3])-100 > as.numeric(allmem[sub,3])){
                    onetrio <- c(allmem[sub,1:2],allmem[k,1])
                    trios <- rbind(trios,onetrio)
                }
                
                if(as.numeric(allmem[k,3])-100 < as.numeric(allmem[sub,3])){
                    if(length(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3]))) > 0){
                        psub <- max(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3])))
                        onetrio <- c(famlist[psub,1:2],allmem[k,1])
                        trios <- rbind(trios,onetrio)
                    }
                }
            }
            
            k <- j
            if(!all(c(famlist[i+1,1],famlist[i+1,2]) %in% allmem[1:j,1])){
                if(as.numeric(allmem[k,3])-100 > as.numeric(allmem[sub,3])){
                    onetrio <- c(allmem[sub,1:2],allmem[k,1])
                    trios <- rbind(trios,onetrio)
                }
                if(as.numeric(allmem[k,3])-100 < as.numeric(allmem[sub,3])){
                    psub <- max(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3])))
                    onetrio <- c(famlist[psub,1:2],allmem[k,1])
                    trios <- rbind(trios,onetrio)
                }
            }
            
        }
    }
    
    if(i == npar){
        j <- length(allV)
        for(k in (sub+1):j){
            if(as.numeric(allmem[k,3])-100 > as.numeric(allmem[sub,3])){
                onetrio <- c(allmem[sub,1:2],allmem[k,1])
                trios <-  rbind(trios,onetrio)
            }else{
                psub <- max(which(as.numeric(allmem[k,3])-100 > as.numeric(famlist[1:i,3])))
                onetrio <- c(famlist[psub,1:2],allmem[k,1])
                trios <- rbind(trios,onetrio)
            }
        }
    }
    
}

