function nucleocode(nucleotide)
{
  if ( nucleotide ~ /^[aA1]$/ ) return 1;
  else if ( nucleotide ~ /^[cC2]$/ ) return 2;
  else if ( nucleotide ~ /^[gG3]$/ ) return 3;
  else if ( nucleotide ~ /^[tT4]$/ ) return 4;
  else return 0;
}

function comp_nucleocode(nucleotide)
{
  if ( nucleotide ~ /^[aA1]$/ ) return 4;
  else if ( nucleotide ~ /^[cC2]$/ ) return 3;
  else if ( nucleotide ~ /^[gG3]$/ ) return 2;
  else if ( nucleotide ~ /^[tT4]$/ ) return 1;
  else return 0;
}

function A_to_i(nucleotide)
{
    if ( nucleotide == "A" ) return 1;
    else if ( nucleotide == "C" ) return 2;
    else if ( nucleotide == "G" ) return 3;
    else if ( nucleotide == "T" ) return 4;
    else return 0;
}

function a_to_i(nucleotide)
{
    if ( nucleotide == "a" ) return 1;
    else if ( nucleotide == "c" ) return 2;
    else if ( nucleotide == "g" ) return 3;
    else if ( nucleotide == "t" ) return 4;
    else return 0;
}

function i_to_A(nucleotide)
{
    if ( nucleotide == "1" ) return "A";
    else if ( nucleotide == "2" ) return "C";
    else if ( nucleotide == "3" ) return "G";
    else if ( nucleotide == "4" ) return "T";
    else return "0";
}

function i_to_a(nucleotide)
{
    if ( nucleotide == "1" ) return "a";
    else if ( nucleotide == "2" ) return "c";
    else if ( nucleotide == "3" ) return "g";
    else if ( nucleotide == "4" ) return "t";
    else return "0";
}

function sepcode()
{
  return 5;
}
