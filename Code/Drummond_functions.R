#from supplement for Drummond DA, Raval A, Wilke CO (2006) A single determinant dominates the rate of yeast protein evolution. Molecular Biology and Evolution 23: 327–337.

partial.cor.test <- function(v1, v2, v3) {
  # eliminate NAs
  z12 = !is.na(v1) & !is.na(v2)
  z13 = !is.na(v1) & !is.na(v3)
  z23 = !is.na(v2) & !is.na(v3)
  #calculate correlation
  c12 <- cor(v1[z12], v2[z12])
  c23 <- cor(v2[z23], v3[z23])
  c13 <- cor(v1[z13], v3[z13])

  r <- (c12-(c13*c23))/(sqrt(1-(c13^2)) * sqrt(1-(c23^2)))

  # be conservative with the length estimate
  n <- length(v1[!is.na(v1) & !is.na(v2) & !is.na(v3)])
  t <- r*sqrt((n-3)/(1-r^2))
  p.value <- 2*(1-pnorm(abs(t)))
  l <- data.frame(estimate=r, p.value=p.value, n=n, c12=c12, c13=c13, c23=c23)
  l
}

mysum <- function( g ){
  # function to normalize the components
  # prints actually the squared values, to get percentages
  f <- function(i, g) {round( g$projection[,i]*g$projection[,i]/(g$projection[,i]%*% g$projection[,i]), 3)}
  # print the summary
  h <- summary( g )
  htmp <- h
  # number of components
  comps = length( attr(h, "dimnames")[[2]] )
  # print summary with percentage differences
  for ( i in 2:comps ){
    h[[1,i]] <- htmp[[1,i]]-htmp[[1,i-1]]
    h[[2,i]] <- htmp[[2,i]]-htmp[[2,i-1]]
  }
  cat( "\nPercentage differences:\n" )
  print( round( h, 2 ), print.gap=2 )

  cat( "\nPercentage contributions to components:\n" )

  # print the normalized projection
  p = g$projection
  for ( i in 1:comps ){
    x <- p[,i]
    p[,i] <- x^2/t(x) %*% x
  }
  print( round( p, 3 ) )

 	# print total variance contributions
	cat( "\nTotal variance explained by each variable\n")
	total_var_explained = as.array(0*(1:comps))
	#print(dimnames(total_var_explained))
	#dimnames(total_var_explained) <- paste('comp', 1:comps)
	for (i in 1:comps) {
		total_var_explained[i] <- sum(h[2,] %*% p[i,])
	}
	print( round(total_var_explained, 3) )
	
	#bars <- as.vector(h[[2,]]) %*% p
	#bars
	bars <- p
	for (i in 1:comps) {
		bars[,i] <- h[2,i] * p[,i]
	}
	for (i in 1:comps) {
		bars[,i] <- sort(bars[,i])
	}
	bars
}