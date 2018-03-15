myfun.sample <- function(x, ...){
        if(all(is.na(x))){
	            return(NA)
    }
    return(sample(x[!is.na(x)], ...))
}
