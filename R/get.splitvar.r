get.splitvar <- function(x,sq=1:length(x$trees),parm="Splitvar") {
    get.ID <- function (m.tr,n=1,id=NULL,which="Splitvar") {
        i <- k <- m <- 1
        ID<-vector()
        if(class(m.tr)[1] == "TWIX" )
            m.tree <- m.tr$multitree
        else if (class(m.tr)[1] == "bootTWIX")
            m.tree <- m.tr$multitree[[n]]
        if( is.null(id))
            tree.id <- m.tr$trees[[n]]$id
        else
            tree.id <- id
        ausgabe <- function(m.tree){
            k<-m
            root <- i;
            if("Splitvar" == which)
                ID[m] <<-m.tree$split[tree.id[root]][[1]]$Splitvar
            if("Dev" == which)
                ID[m] <<-m.tree$split[tree.id[root]][[1]]$Dev
			if("Split" == which)
                ID[m] <<-m.tree$split[tree.id[root]][[1]]$Splitp
            if (tree.id[i<<-i+1]!=0){
                m<<-2*k
                ausgabe(m.tree$left[[tree.id[root]]])
            }
            if (tree.id[i<<-i+1]!=0){
                m<<-2*k+1
                ausgabe(m.tree$right[[tree.id[root]]])
            }
        }
        ausgabe(m.tree)
        ID
    }
    m<-DT<-list()
    p<-1
    for(i in sq){
        m[[p]]<-get.ID(x,n=i,which=parm)
        p<-p+1
        }
    MM<-matrix(,length(m),max(sapply(m,length)))
    for(i in 1:length(m)){
        for(j in 1:length(m[[i]])) {
            MM[i,j]<-m[[i]][[j]]
        }
    }
    #data.frame(MM)
    MM
}
