export <- function(x,sq=1,directory="ForKlimt") {
    if (inherits(x, "TWIX") | inherits(x, "bootTWIX")) {
        if(Sys.info()[[1]] == "Windows"){
            shell(paste("md",directory))
            #system(paste("mkdir",directory))
        } else {
            system(paste("mkdir",directory))
            }
        dir <- getwd()
        setwd(paste(dir,paste("/",directory,sep=""),sep=""))
        for(i in sq) {
            write.table(NULL,paste(as.character(i),".tree",sep=""),row.names=FALSE)
            sink(paste(as.character(i),".tree",sep=""),TRUE)
            print(get.tree(x,i))
            sink()
        }
        setwd(dir)
    }
}
 
