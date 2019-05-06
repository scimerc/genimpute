# sample mean
function smean ( x )
{
  cnt = 0
  sum = 0
  for ( t in x ) {
    cnt++
    sum += x[t]
  }
  return sum / cnt
}

# sample standard deviation
function ssd ( x )
{
  cnt = 0
  sum = 0
  for ( t in x ) {
    cnt++
    sum += ( x[t] - smean(x) )^2
  }
  return sum / cnt
}

# pairwise minimum
function pmin ( x, y )
{
  if ( x < y ) return x
  else return y
}

# pairwise maximum
function pmax ( x, y )
{
  if ( x > y ) return x
  else return y
}

# sample minimum
function min ( x )
{
  cnt = 0
  for ( t in x ) {
    if ( cnt == 0 ) tmp = x[t]
    else tmp = min( tmp, x[t] )
    cnt++
  }
}

# sample maximum
function max ( x )
{
  cnt = 0
  for ( t in x ) {
    if ( cnt == 0 ) tmp = x[t]
    else tmp = max( tmp, x[t] )
    cnt++
  }
}

