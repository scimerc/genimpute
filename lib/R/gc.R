gc_pval <- function ( pvalstat, filter=NULL ) {
    chisqstat <- qchisq( pvalstat, 1, lower.tail = F )
    if( is.null( filter ) ) {
      gc.lambda <- median( chisqstat, na.rm = T ) / qchisq( 0.5, 1 )
    } else {
      gc.lambda = c()
      Nf = dim(filter)[2]
      if ( is.na(Nf) ) Nf = 1
      for ( i in 1:Nf ) {
	gc.lambda <- c( gc.lambda, median( chisqstat[ filter[,i] & is.finite( chisqstat ) ], na.rm=T ) / qchisq( 0.5, 1 ) )
      }
      gc.lambda = median( gc.lambda, na.rm=T )
    }
    #print( gc.lambda )
    return( pchisq( chisqstat / gc.lambda, 1, lower.tail = F ) )
}
