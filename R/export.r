export<-function (x, sq = 1, directory = "ForKlimt"){
	dir <- getwd()
    if(Sys.info()[[1]] == "Windows"){
		shell(paste("md", directory))
	}
	else{
		system(paste("mkdir", directory))
	}
    setwd(paste(dir, paste("/", directory, sep = ""), sep = ""))
    if(inherits(x, "TWIX") | inherits(x, "bootTWIX")){
        for(i in sq){
			write.table(NULL, paste(as.character(i), ".tree", 
                sep = ""), row.names = FALSE, quote = FALSE)
            sink(paste(as.character(i), ".tree", sep = ""), TRUE)
            print(get.tree(x, i))
            sink()
        }
        setwd(dir)
    }
   	if(inherits(x,"single.tree")){
   		write.table(NULL, paste("FromR.tree",
   			sep = ""), row.names = FALSE, quote = FALSE)
        sink(paste("FromR.tree", ".tree", sep = ""), TRUE)
        print(x)
        sink()
        setwd(dir)
	}
}
