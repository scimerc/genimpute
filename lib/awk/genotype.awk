# @include ( "nucleocode.awk" )
function genotype(alleleA,alleleB)
{
  genobase = 1 + sepcode();
  genocode = 0;
  genotuple[1] = alleleA;
  genotuple[2] = alleleB;
  asort( genotuple );
  delete( genovec );
  split( genotuple[1], genovec, "" )
  for ( k = 1; k <= length(genovec); k++ )
    genocode += nucleocode(genovec[k]) * genobase^(k-1);
  genocode += sepcode() * genobase^(k-1);
  delete( genovec );
  split( genotuple[2], genovec, "" )
  for ( h = 1; h <= length(genovec); h++ )
    genocode += nucleocode(genovec[h]) * genobase^(k+h-1);
  return( genocode );
}
