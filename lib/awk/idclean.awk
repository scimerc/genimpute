function idclean(x)
{
  gsub( "[][)(}{/\\,.;:|!?@#$%^&*~=_><+-]+", "_", x )
  return(x)
}
