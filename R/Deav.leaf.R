    Dev.leaf <- function(x)
    {
        TD <- 0
        CCR <- table(x)
        s <- sum(CCR)
        TD <- -sum(s*sapply(CCR,function(x,y){ if(x/y == 0) 0 else (x/y)*log(x/y)},y=s))
        round(TD,digits=6)
    }
