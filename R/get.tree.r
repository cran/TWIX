get.tree <- function (m.tr,n=1,id=NULL) {
    i <- 1
    if(class(m.tr) == "TWIX" )
        m.tree <- m.tr$multitree
    else if (class(m.tr) == "bootTWIX")
        m.tree <- m.tr$multitree[[n]]
    if( is.null(id))
        tree.id <- m.tr$trees[[n]]$id
    else
        tree.id <- id
    tree <- list()
    ausgabe <- function(m.tree){
    root <- i
    tree <- list(split=m.tree$split[tree.id[root]],
                left=if (tree.id[i<<-i+1]!=0)
                        ausgabe(m.tree$left[[tree.id[root]]])
                     else
                        if(is.null(m.tree$left[[tree.id[root]]]$split))
                            Prob=m.tree$left[[tree.id[root]]]
                        else
                            Prob=m.tree$left[[tree.id[root]]]$split[[1]]
                ,right=if (tree.id[i<<-i+1]!=0)
                        ausgabe(m.tree$right[[tree.id[root]]])
                     else
                        if(is.null(m.tree$right[[tree.id[root]]]$split))
                            Prob=m.tree$right[[tree.id[root]]]
                        else
                            Prob=m.tree$right[[tree.id[root]]]$split[[1]]
                )
        tree
     }
    if(tree.id[1] != 0  ){
        tree <- list(ausgabe(m.tree),m.tr$formula)
    }
    else{
        Obs <-m.tree$split[[1]]$Obs + length(m.tr$Bad.id)
        prob <- m.tree$split[[1]]$Prob
        yval <- m.tree$split[[1]]$Pred.class
        tree <- list(0,m.tr$formula,Obs,0,prob,yval)
    }
    class(tree) <- "single.tree"
    tree
}
