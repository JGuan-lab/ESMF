`normalize.col.in.matrix` <-
function (m,col.value=1,method="mean")
{

	
	if (method=="mean") {
		mean.diff.c <- c()
		for (col in 1: ncol(m)){
		
			mean.old <- mean(m[,col],na.rm = TRUE)
			mean.diff <- mean.old / col.value
		
			index <-  m[,col] != 0
            	#shift the total value line by the diff value, doing this, the new mean value
			#will be the new average value
			m[index,col] <- m[index,col]/ mean.diff
			mean.diff.c <- c(mean.diff.c,mean.diff )
		}
		r.list <- list()
		r.list[["Matrix"]] <- m
		r.list[["Values"]] <- mean.diff.c
		return (r.list)
	}
	if (method=="sum") {
		#use the mean method because the sum = mean * nrows
		col.value <- col.value / nrow(m)
		return (normalize.col.in.matrix(m,col.value,method="mean"))

	}
}

