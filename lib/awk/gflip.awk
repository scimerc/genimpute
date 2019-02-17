# @include "nucleocode.awk"
# @include "genotype.awk"
function gflip(a1,a2,b1,b2)
{
  a1len = split( a1, a1vec, "" );
  a2len = split( a2, a2vec, "" );
  b1len = split( b1, b1vec, "" );
  b2len = split( b2, b2vec, "" );
  if ( a1len + a2len != b1len + b2len ) return(0);
  else {
    a1code = "";
    a2code = "";
    ac1code = "";
    ac2code = "";
    b1code = "";
    b2code = "";
    bc1code = "";
    bc2code = "";
    for ( k = 1; k <= a1len; k++ ) {
      if ( nucleocode(a1vec[k]) == 0 ) return(0);
      else {
        a1code = a1code nucleocode(a1vec[k]);
        ac1code = ac1code comp_nucleocode(a1vec[k]);
      }
    }
    for ( k = 1; k <= a2len; k++ ) {
      if ( nucleocode(a2vec[k]) == 0 ) return(0);
      else {
        a2code = a2code nucleocode(a2vec[k]);
        ac2code = ac2code comp_nucleocode(a2vec[k]);
      }
    }
    for ( k = 1; k <= b1len; k++ ) {
      if ( nucleocode(b1vec[k]) == 0 ) return(0);
      else {
        b1code = b1code nucleocode(b1vec[k]);
        bc1code = bc1code comp_nucleocode(b1vec[k]);
      }
    }
    for ( k = 1; k <= b2len; k++ ) {
      if ( nucleocode(b2vec[k]) == 0 ) return(0);
      else {
        b2code = b2code nucleocode(b2vec[k]);
        bc2code = bc2code comp_nucleocode(b2vec[k]);
      }
    }
    gen_a = genotype(a1code,a2code);
    gen_ac = genotype(ac1code,ac2code);
    gen_b = genotype(b1code,b2code);
    gen_bc = genotype(bc1code,bc2code);
    if ( gen_a == gen_bc && gen_ac == gen_b && gen_a != gen_ac && gen_b != gen_bc ) return(1);
    else return(0);
  }
}

function gflipx(a1,a2,b1,b2)
{
  a1len = split( a1, a1vec, "" );
  a2len = split( a2, a2vec, "" );
  b1len = split( b1, b1vec, "" );
  b2len = split( b2, b2vec, "" );
  if ( a1len + a2len != b1len + b2len ) return(0);
  else {
    a1code = "";
    a2code = "";
    b1code = "";
    b2code = "";
    for ( k = 1; k <= a1len; k++ )
      a1code = a1code nucleocode(a1vec[k])
    for ( k = 1; k <= a2len; k++ )
      a2code = a2code nucleocode(a2vec[k])
    for ( k = 1; k <= b1len; k++ )
      b1code = b1code nucleocode(b1vec[k])
    for ( k = 1; k <= b2len; k++ )
      b2code = b2code nucleocode(b2vec[k])
    if ( a1len != b1len ) {
      dmcode = b1code;
      b1code = b2code;
      b2code = dmcode;
    }
    bflip = 1;
    brflip = 1;
    binfo = 0;
    brinfo = 0;
    split( a1code a2code, avec, "" );
    split( b1code b2code, bvec, "" );
    split( b2code b1code, brvec, "" );
    for ( k = 1; k <= a1len + a2len; k++ ) {
      if ( avec[k] != 0 && bvec[k] != 0 ) {
        binfo = 1;
        if ( nucleocode(avec[k]) != comp_nucleocode(bvec[k]) ) bflip = 0;
      }
      if ( avec[k] != 0 && brvec[k] != 0 ) {
        brinfo = 1;
        if ( nucleocode(avec[k]) != comp_nucleocode(brvec[k]) ) brflip = 0;
      }
    }
    if ( bflip == 1 && binfo == 1 || brflip == 1 && brinfo == 1 ) return(1);
    else return(0);
  }
}
