args <- commandArgs(T)
kk <- as.numeric(args[1])

source("Faminfo.R")
if(kk <= 20){
    run_famSKAT(kk,FALSE,1)
}

if(kk > 20){
    kk <- kk-20
    run_famSKAT(kk,TRUE,1)
}

