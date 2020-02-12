# @include "nucleocode.awk"
# @include "genotype.awk"
function gflip(a1,a2,b1,b2)
{
  # store variant nucleotide components in vectors
  a1len = split( a1, a1vec, "" );
  a2len = split( a2, a2vec, "" );
  b1len = split( b1, b1vec, "" );
  b2len = split( b2, b2vec, "" );
  # no flip if the two variant versions have different cumulative lengths
  if ( a1len + a2len != b1len + b2len ) return(0);
  # initialize variant code strings
  a1code = "";
  a2code = "";
  ac1code = "";
  ac2code = "";
  b1code = "";
  b2code = "";
  bc1code = "";
  bc2code = "";
  # assemble variant code strings and complement strings
  for ( k = 1; k <= a1len; k++ ) {
    if ( nucleocode(a1vec[k]) == 0 ) return(0);
    a1code = a1code nucleocode(a1vec[k]);
    ac1code = ac1code comp_nucleocode(a1vec[k]);
  }
  for ( k = 1; k <= a2len; k++ ) {
    if ( nucleocode(a2vec[k]) == 0 ) return(0);
    a2code = a2code nucleocode(a2vec[k]);
    ac2code = ac2code comp_nucleocode(a2vec[k]);
  }
  for ( k = 1; k <= b1len; k++ ) {
    if ( nucleocode(b1vec[k]) == 0 ) return(0);
    b1code = b1code nucleocode(b1vec[k]);
    bc1code = bc1code comp_nucleocode(b1vec[k]);
  }
  for ( k = 1; k <= b2len; k++ ) {
    if ( nucleocode(b2vec[k]) == 0 ) return(0);
    b2code = b2code nucleocode(b2vec[k]);
    bc2code = bc2code comp_nucleocode(b2vec[k]);
  }
  # get genotype codes (order-naive)
  gen_a = genotype(a1code,a2code);
  gen_ac = genotype(ac1code,ac2code);
  gen_b = genotype(b1code,b2code);
  gen_bc = genotype(bc1code,bc2code);
  # return true if the genotype codes are complementary (but not self-complementary)
  if ( gen_a == gen_bc && gen_ac == gen_b && gen_a != gen_ac && gen_b != gen_bc ) return(1);
  return(0);
}

function gflipx(a1,a2,b1,b2)
{
  # store variant nucleotide components in vectors
  a1len = split( a1, a1vec, "" );
  a2len = split( a2, a2vec, "" );
  b1len = split( b1, b1vec, "" );
  b2len = split( b2, b2vec, "" );
  # no flip if the two variant versions have different cumulative lengths
  if ( a1len + a2len != b1len + b2len ) return(0);
  # initialize variant code strings
  a1code = "";
  a2code = "";
  b1code = "";
  b2code = "";
  # assemble variant code strings
  for ( k = 1; k <= a1len; k++ )
    a1code = a1code nucleocode(a1vec[k])
  for ( k = 1; k <= a2len; k++ )
    a2code = a2code nucleocode(a2vec[k])
  for ( k = 1; k <= b1len; k++ )
    b1code = b1code nucleocode(b1vec[k])
  for ( k = 1; k <= b2len; k++ )
    b2code = b2code nucleocode(b2vec[k])
  # b1 and b2 are reversed compared to a1 and a2 if the lengths of a1 and b1 don't match
  if ( a1len != b1len ) {
    dmcode = b1code;
    b1code = b2code;
    b2code = dmcode;
  }
  # initialize match flags (assume match by default)
  bmatch = 1;
  brmatch = 1;
  # initialize flip flags (assume flip by default)
  bflip = 1;
  binfo = 0;
  # initialize reverse flip flags (b1 and b2 could still be reversed if they are equally long
  brflip = 1;
  brinfo = 0;
  # store the whole variant strings in vectors
  split( a1code a2code, avec, "" );
  split( b1code b2code, bvec, "" );
  split( b2code b1code, brvec, "" );
  for ( k = 1; k <= a1len + a2len; k++ ) {
    if ( avec[k] != 0 && bvec[k] != 0 ) {
      binfo = 1; # we have full information at position k
      # rule out flip if the two corresponding nuclotides do not complement each other
      if ( nucleocode(avec[k]) != comp_nucleocode(bvec[k]) ) bflip = 0;
      if ( nucleocode(avec[k]) != nucleocode(bvec[k]) ) bmatch = 0;
    }
    if ( avec[k] != 0 && brvec[k] != 0 ) {
      brinfo = 1; # we have full reverse information at position k
      # rule out reverse flip if the two corresponding nuclotides do not complement each other
      if ( nucleocode(avec[k]) != comp_nucleocode(brvec[k]) ) brflip = 0;
      if ( nucleocode(avec[k]) != nucleocode(brvec[k]) ) brmatch = 0;
    }
  }
  # return true if flip or reverse flip wasn't ruled out but matches were
  if ( bmatch == 0 && brmatch == 0 && ( bflip == 1 && binfo == 1 || brflip == 1 && brinfo == 1 ) ) return(1);
  return(0);
}

