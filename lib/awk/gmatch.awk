# @include "nucleocode.awk"
# @include "genotype.awk"
function gmatch(a1,a2,b1,b2)
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
    c1code = "";
    c2code = "";
    for ( k = 1; k <= a1len; k++ ) {
      if ( nucleocode(a1vec[k]) == 0 ) return(0);
      else a1code = a1code nucleocode(a1vec[k]);
    }
    for ( k = 1; k <= a2len; k++ ) {
      if ( nucleocode(a2vec[k]) == 0 ) return(0);
      else a2code = a2code nucleocode(a2vec[k]);
    }
    for ( k = 1; k <= b1len; k++ ) {
      if ( nucleocode(b1vec[k]) == 0 ) return(0);
      else {
        b1code = b1code nucleocode(b1vec[k]);
        c1code = c1code comp_nucleocode(b1vec[k]);
      }
    }
    for ( k = 1; k <= b2len; k++ ) {
      if ( nucleocode(b2vec[k]) == 0 ) return(0);
      else {
        b2code = b2code nucleocode(b2vec[k]);
        c2code = c2code comp_nucleocode(b2vec[k]);
      }
    }
    gen_a = genotype(a1code,a2code);
    gen_b = genotype(b1code,b2code);
    gen_c = genotype(c1code,c2code);
    if ( gen_a == gen_b || gen_a == gen_c ) return(1);
    else return(0);
  }
}

function gmatchx(a1,a2,b1,b2)
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
    c1code = "";
    c2code = "";
    for ( k = 1; k <= a1len; k++ )
      a1code = a1code nucleocode(a1vec[k]);
    for ( k = 1; k <= a2len; k++ )
      a2code = a2code nucleocode(a2vec[k]);
    for ( k = 1; k <= b1len; k++ ) {
      b1code = b1code nucleocode(b1vec[k]);
      c1code = c1code comp_nucleocode(b1vec[k]);
    }
    for ( k = 1; k <= b2len; k++ ) {
      b2code = b2code nucleocode(b2vec[k]);
      c2code = c2code comp_nucleocode(b2vec[k]);
    }
    if ( a1len != b1len ) {
      dmcode = b1code;
      b1code = b2code;
      b2code = dmcode;
      dmcode = c1code;
      c1code = c2code;
      c2code = dmcode;
    }
    bmatch = 0;
    cmatch = 0;
    brmatch = 0;
    crmatch = 0;
    split( a1code a2code, avec, "" );
    split( b1code b2code, bvec, "" );
    split( c1code c2code, cvec, "" );
    split( b2code b1code, brvec, "" );
    split( c2code c1code, crvec, "" );
    for ( k = 1; k <= a1len + a2len; k++ ) {
      if ( avec[k] != 0 && bvec[k] != 0 && avec[k] != bvec[k] ) bmatch = -1;
      if ( avec[k] != 0 && cvec[k] != 0 && avec[k] != cvec[k] ) cmatch = -1;
      if ( avec[k] != 0 && brvec[k] != 0 && avec[k] != brvec[k] ) brmatch = -1;
      if ( avec[k] != 0 && crvec[k] != 0 && avec[k] != crvec[k] ) crmatch = -1;
    }
    if ( bmatch == 0 || cmatch == 0 || brmatch == 0 || crmatch == 0 ) return(1);
    else return(0);
  }
}
