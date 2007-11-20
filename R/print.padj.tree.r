print.padj.tree <- function(x,klimt=FALSE,Data=NULL,file="FromR.tree",...) {
    n <- i <- 1
    no.split <- FALSE
    if(!is.list(x[[1]]) && x[[1]] == 0){
        Obs <- x[[3]]
        no.split <- TRUE
    }
    else{ 
        Obs <- x[[1]]$split[[1]]$Obs
        }
    pr.baum <- function(baum,i) {
        anf <- n
        node <- paste( n,")",sep = "" )
        abst <- paste( rep(" ",i))
        split.var <- baum$split[[1]]$Splitvar
        split.p <- baum$split[[1]]$Splitp
        leaf <- "*"
        if (i==1) {
            cat("n=", Obs, "\n\n")
            cat( "node), split, n, p-value for nodes and deviance for leafs, yval, (yprob)\n" )
            cat(" * denotes terminal node\n\n")
            if(no.split){
                dev <- 0
                yprob <- paste(format(as.vector(round(x[[5]],4)),
                    nsmall=1),sep="")
                yval <- paste(x[[6]],"(")
                split.v <- "root"
                cat(node,split.v,Obs," ",dev,yval,yprob,")",leaf,"\n")
            }
            else{
                dev <- format(1-baum$split[[1]]$Dev,4)
                yprob <- paste(format(as.vector(round(baum$split[[1]]$Prob,4)),
                    nsmall=1),sep="")
                yval <- paste(baum$split[[1]]$Pred.class,"(")
                split.v <- "root"
                cat(node,split.v,Obs," ",dev,yval,yprob,")","\n")
            }
        }
        if (length(baum$left$split) > 0 && !no.split) {
            if(length(split.p) > 1) {
                split.pl<-names(split.p)[split.p == 1]
                n <<- 2*anf
                node <- paste( n,")",sep = "" )
                yval <- paste(baum$left$split[[1]]$Pred.class,"(")
                Obs <- baum$left$split[[1]]$Obs
                dev <- format(1-baum$left$split[[1]]$Dev,4)
                yprob <- paste(format(as.vector(round(baum$left$split[[1]]$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"= ")
                if(length(split.pl)>1)
                    cat(split.pl,sep=",")
                else
                    cat(split.pl)
                cat(" ",Obs,dev,yval,yprob,")","\n")
                pr.baum(baum$left,i<<-i+1)
            }
            else {
                n <<- 2*anf
                node <- paste( n,")",sep = "" )
                yval <- paste(baum$left$split[[1]]$Pred.class,"(")
                Obs <- baum$left$split[[1]]$Obs
                dev <- format(1-baum$left$split[[1]]$Dev,4)
                yprob <- paste(format(as.vector(round(baum$left$split[[1]]$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"<",split.p,Obs,dev,yval,yprob,")","\n")
                pr.baum(baum$left,i<<-i+1)
            }
        }
        else if(!no.split){
            if(length(split.p) > 1){
                split.pl<-names(split.p)[split.p == 1]
                n <<- 2*anf
                node <- paste( n,")",sep = "" )
                Obs <- baum$left$Obs
                yval<-paste(names(baum$left$Prob[sort(baum$left$Prob,
                    index.return=TRUE,decreasing=TRUE)$ix[1]]))
                dev <- baum$left$Dev
                yprob <- paste(format(as.vector(round(baum$left$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"= ")
                if(length(split.pl)>1)
                    cat(split.pl,sep=",")
                else
                    cat(split.pl)
                cat(" ",Obs," ",dev,yval,"(",yprob,")",leaf,"\n")
            }
            else{
                n <<- 2*anf
                node <- paste( n,")",sep = "" )
                Obs <- baum$left$Obs
                yval<-paste(names(baum$left$Prob[sort(baum$left$Prob,
                    index.return=TRUE,decreasing=TRUE)$ix[1]]))
                dev <- baum$left$Dev
                yprob <- paste(format(as.vector(round(baum$left$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"<",split.p,Obs," ",
                    dev,yval,"(",yprob,")",leaf,"\n")
                }
            }
        if (length(baum$right$split) > 0 && !no.split){
            if(length(split.p) > 1){
                split.pr<-names(split.p)[split.p == 0]
                n <<- 2*anf+1
                node <- paste( n,")",sep = "" )
                yval <- paste(baum$right$split[[1]]$Pred.class,"(")
                Obs <- baum$right$split[[1]]$Obs
                dev <- format(1-baum$right$split[[1]]$Dev,4)
                yprob<-paste(format(as.vector(round(baum$right$split[[1]]$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"= ")
                if(length(split.pr)>1)
                    cat(split.pr,sep=",")
                else
                    cat(split.pr)
                cat(" ",Obs,dev,yval,yprob,")","\n")
                pr.baum(baum$right,i<<-i+1)
            }
            else{
                n <<- 2*anf+1
                node <- paste( n,")",sep = "" )
                yval <- paste(baum$right$split[[1]]$Pred.class,"(")
                Obs <- baum$right$split[[1]]$Obs
                dev <- format(1-baum$right$split[[1]]$Dev,4)
                yprob<-paste(format(as.vector(round(baum$right$split[[1]]$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,">=",split.p,Obs,
                    dev,yval,yprob,")","\n")
                pr.baum(baum$right,i<<-i+1)
            }
        }
        else if(!no.split){
            if(length(split.p) >1) {
                split.pr<-names(split.p)[split.p == 0]
                n <<- 2*anf+1
                node <- paste( n,")",sep = "" )
                Obs <- baum$right$Obs
                yval<-paste(names(baum$right$Prob[sort(baum$right$Prob,
                    index.return=TRUE,decreasing=TRUE)$ix[1]]))
                dev <- baum$right$Dev
                yprob <- paste(format(as.vector(round(baum$right$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,"=")
                if(length(split.pr)>1)
                    cat(split.pr,sep=",")
                else
                    cat(split.pr)
                cat(" ",Obs," ",dev,yval,"(",yprob,")",leaf,"\n")
            }
            else {
                n <<- 2*anf+1
                node <- paste( n,")",sep = "" )
                Obs <- baum$right$Obs
                yval<-paste(names(baum$right$Prob[sort(baum$right$Prob,
                    index.return=TRUE,decreasing=TRUE)$ix[1]]))
                dev <- baum$right$Dev
                yprob <- paste(format(as.vector(round(baum$right$Prob,4)),
                    nsmall=1),sep="")
                cat(abst,node,split.var,">=",split.p,Obs," ",
                    dev,yval,"(",yprob,")",leaf,"\n")
            }
        }
    }
    if (klimt) {
        if (!is.data.frame(Data)) warning("\n Data is not a data.frame!\n")
        Data <- na.omit(Data)
        write.table(Data,file,row.names=FALSE,sep="\t",quote=FALSE)
        sink(file,TRUE)
        pr.baum(x[[1]],1)
        print(x[[2]])
        sink()
        system(paste("java -jar Klimt.jar",file))
    }
    else {
        pr.baum(x[[1]],1)
        cat("\n")
        print(x[[2]])
    }
    invisible(x)
}

 
