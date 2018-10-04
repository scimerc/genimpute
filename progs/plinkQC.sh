#!/bin/bash

# exit on error
trap 'exit' ERR

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

AWKPATH="${AWKPATH}:${BASEDIR}/lib/awk"
AWKLOCINCLUDE=$( printf -- '-f %s\n' $( echo ${AWKINCLUDE} \
    abs.awk \
    nucleocode.awk \
    genotype.awk \
    gflip.awk \
    gmatch.awk \
    round.awk \
) | sort -u )
alias join='join --check-order'
fdrscript="R --slave --file=${BASEDIR}/progs/fdr.Rscript"
hvmscript="R --slave --file=${BASEDIR}/progs/het_vs_miss.Rscript"
plinkexe=$( ( which plink || true ) 2> /dev/null )
if [[ "${plinkexe}" == "" ]] ; then
    echo "plink is required by $( basename $0 ). plink source codes and builds can be found at"
    echo "www.cog-genomics.org. note that some of the functionalities needed by $( basename $0 )"
    echo "were not implemented in plink2 at the time of writing."
    exit
fi

opt_parser="${BASEDIR}/lib/sh/opt_parser.bash"
opt_list=("b=" "dups" "g=" "i=" "maf=" "nohvm" "o=" "p=" "r=" "s=" "v=" "x" "h")
get_opt ()
{
    case $1 in
        "b" )
            blacklist=`echo "$2" | sed "s/^=\+//"` ;;
        "dups" )
            keepdups='on' ;;
        "g" )
            genome=`echo "$2" | sed "s/^=\+//"` ;;
        "i" )
            samplefile=`echo "$2" | sed "s/^=\+//"` ;;
        "maf" )
            myfreq_std=`echo "$2" | sed "s/^=\+//"` ;;
        "nohvm" )
            hvm='off' ;;
        "o" )
            myprefix=`echo "$2" | sed "s/^=\+//"` ;;
        "p" )
            phenotypes=`echo "$2" | sed "s/^=\+//"` ;;
        "r" )
            refalleles=`echo "$2" | sed "s/^=\+//"` ;;
        "s" )
            samplemiss=`echo "$2" | sed "s/^=\+//"` ;;
        "v" )
            varmiss=`echo "$2" | sed "s/^=\+//"` ;;
        "x" )
            mini="yes" ;;
        "h" )
            helpme='yes' ;;
    esac
}
helpme=''
blacklist=''
genome='hg19'
hvm='on'
hweneglogp=12
hweneglogp_ctrl=4
hwtag='midp include-nonctrl'
keepdups='off'
mincnt=100
mini='no'
myfreq_hq=0.2
myfreq_std=0.05
myprefix='plink'
phenotypes=''
pihat=0.9
refalleles=''
samplefile=''
samplemiss=0.05
varmiss=0.05
tmpargs=''
while [[ "${tmpargs}" == "" ]] ; do
    tmpargs=$( mktemp .tmpXXXXXXXX )
done
source ${opt_parser} > ${tmpargs} || true
mybatches=$( cat ${tmpargs} )
rm -f ${tmpargs}

if [[ "${mybatches}" == "" || "${helpme}" != "" ]] ; then

    if [[ "${helpme}" == "" ]] ; then echo "missing input."; fi

    echo -e "\n USAGE:"
    echo -e "     $( basename $0 ) [OPTIONS] <bed|bcf|vcf file(s)>\n"
    echo -e "   where\n"
    echo -e "     <bed|bcf|vcf file(s)> are the genotype files to be merged."
    echo -e "\n OPTIONS:"
    echo -e "     -b <black list>             recommended list of bad (high LD) regions (bed format)"
    echo -e "     -dups                       keep duplicate individuals [default: off]"
    echo -e "     -g <genome version>         genome version [default: hg19]"
    echo -e "     -i <sample file>            optional sample selection file"
    echo -e "     -maf <freq>                 minor allele frequency filter"
    echo -e "     -nohvm                      turns off heterozygosity VS missingness tests"
    echo -e "     -o <output prefix>          optional output prefix [default: 'plink']"
    echo -e "     -p <phenotype file>         optional file with 'affected'/'control' status tags"
    echo -e "     -r <reference alleles>      tab separated reference alleles [format: <CHR:POS> <A1> <A2>]"
    echo -e "     -s <sample missingness>     maximum fraction of variants missing in an individual genotype"
    echo -e "     -v <variant missingness>    maximum fraction of genotypes missing for an individual variant"
    echo -e "     -h                          print help\n"

else

    module load R

    echo -e "==== preMaCH.sh -- $(date)\n"
    echo -e "==== options in effect:\n"
    echo "     blacklist=${blacklist}"
    echo "     keepdups=${keepdups}"
    echo "     genome=${genome}"
    echo "     samplefile=${samplefile}"
    echo "     hvm=${hvm}"
    echo "     myprefix=${myprefix}"
    echo "     phenotypes=${phenotypes}"
    echo "     refalleles=${refalleles}"
    echo "     samplemiss=${samplemiss}"
    echo "     varmiss=${varmiss}"
    echo -e "\n==================================\n"
    echo -e "batch files:\n$( ls ${mybatches} )"
    echo -e "\n==================================\n"

    join --version
    echo
    R --version

    mydir=$( dirname ${myprefix} )
    mkdir -p ${mydir}

    bedhex='1b6c ff01'
    genhex='1234 1234'
    bcfhex='4342 0246'
    vcfhex='2323 6966'
    fflag='--bfile'
    myarchvec=( $mybatches )
    mybatchvec=( $mybatches )
    outprefix=${myprefix}
    if (( ${#mybatchvec[*]} > 1 )) ; then
        outprefix=${myprefix}_merged
    fi
    excludefile=${outprefix}.exclude
    extractfile=${outprefix}.extract
    echo -n > ${excludefile}
    echo -n > ${extractfile}
    init=1
    for (( i = 0; i < ${#mybatchvec[@]}; i++ )); do
        myarchvec[$i]=$( zcat -f ${mybatchvec[$i]} | hexdump | head -n 1 | cut -d ' ' -f 2,3 )
        if [[ "${mini}" == "yes" ]] ; then
            case ${myarchvec[$i]} in
                "${bedhex}" )
                    fflag='--bim' ;;
                "${bcfhex}" )
                    fflag='--bcf' ;;
                "${vcfhex}" )
                    fflag='--vcf' ;;
            esac
            extmp=$( mktemp -u .tmpXXXXXXX )
            plink ${fflag} ${mybatchvec[$i]/%.bed/.bim} --make-just-bim --out ${extmp}
            if [[ -s "${extmp}.bim" ]] ; then
                if (( init == 1 )) ; then
                    awk '{ OFS="\t"; print( $1, $4 - 1, $4 ); }' ${extmp}.bim > ${extmp}.bed
                else
                    bedtools intersect \
                        -a ${extractfile} \
                        -b <( awk '{ OFS="\t"; print( $1, $4-1, $4, $2 ); }' ${extmp}.bim ) \
                    > ${extmp}.bed
                    mv ${extmp}.bed ${extractfile}
                    if [[ -f "${extmp}.bim" ]] ; then
                        rm -f ${extmp}.bim
                    fi
                fi
                init=0
            fi
        fi
    done
    mybatchfile=${myprefix}_batches
    ls ${mybatches} | sed -r 's/[.]bed$//g;' > ${mybatchfile}
    mybatchvec=( $( cat ${mybatchfile} ) )

    if [[ "${phenotypes}" == "" || ! -f "${phenotypes}" ]] ; then
        hweneglogp=${hweneglogp_ctrl}
    fi

    plinkflagsdef=""
    if [[ "${samplefile}" != "" && -f "${samplefile}" ]] ; then
        plinkflagsdef="--keep ${samplefile}"
    fi

    if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
        procdir=${myprefix}_qc.proc
        mkdir -p ${procdir}
        goflag=""
        while [[ "${goflag}" == "" ]] ; do
            goflag="go"
            tmpbatchfile=$( mktemp .tmpbatch.XXXXXX )
            echo "created temporary batch file ${tmpbatchfile}."
            echo -n > ${tmpbatchfile}
            for batch in $( cat ${mybatchfile} ) ; do
                echo "processing batch ${batch}.."
                tmpfile=$( mktemp .tmpctrl.XXXXXX )
                echo "created temporary variant control file ${tmpfile}."
                batchctrl=${procdir}/$( basename ${batch} ).fatto
                batchexcludefile=${procdir}/$( basename ${batch} ).exclude
                batchflipfile=${procdir}/$( basename ${batch} ).flip
                if [[ -s "${excludefile}" ]] ; then
                    sort -u ${excludefile} >> ${batchexcludefile}
                fi
                j=-1
                tmpbatch=${batch}
                echo -n "current batch is number "
                for (( i = 0; i < ${#mybatchvec[@]}; i++ )); do
                    if [[ "${tmpbatch}" == "${mybatchvec[$i]}" ]] ; then
                        j=${i}
                    fi
                done
                echo $j
                fflag='--bfile'
                if (( j > 0 )) ; then
                    case ${myarchvec[$j]} in
                        "${bedhex}" )
                            fflag='--bfile' ;;
                        "${bcfhex}" )
                            fflag='--bcf' ;;
                        "${vcfhex}" )
                            fflag='--vcf' ;;
                    esac
                    unset myarchvec[$j]
                    unset mybatchvec[$j]
                fi
                if [[ ! -f "${batchctrl}" ]] ; then
                    if [[ -s "${tmpbatch}.bim" ]] ; then
                        for chr in $( cut -f 1 ${tmpbatch}.bim | sort -u ) ; do
                            exflag="--chr ${chr}"
                            if [[ -s "${extractfile}" ]] ; then
                                exflag="${exflag} --extract range ${extractfile}"
                            fi
                            tmpvardups=$( mktemp .tmpchr${chr}.XXXXXX )
                            chrbatchprefix=${procdir}/$( basename ${tmpbatch} )_chr${chr}
                            while [[ ! -f "${chrbatchprefix}.tped" ]] ; do
                                plink $fflag ${tmpbatch} ${exflag} --recode transpose --out ${chrbatchprefix}
                            done
                            # extract marker information corresponding to unique position-specific genotype series:
                            # identical position-specific genotype series will thereby be ignored as harmless here.
                            sort -t ' ' -u -k 4 ${chrbatchprefix}.tped | awk '{
                                delete catalog
                                for ( k = 5; k <= NF; k++ ) {
                                    if ( $k > 0 ) catalog[$k] = 1
                                }
                                asorti( catalog )
                                printf( "%s %s %s %s", $1, $2, $3, $4 )
                                for ( allele in catalog ) printf( "_%s", catalog[allele] )
                                print( "" )
                            }' | sort -t ' ' -k 4,4 > ${chrbatchprefix}.gp
                            echo "extracting marker information for unique position-specific genotype series.."
                            cut -d ' ' -f 3,4 ${chrbatchprefix}.gp | sort -t ' ' | uniq -d | sort -t ' ' -k 2,2 \
                            | join -t ' ' -1 4 -2 2 ${chrbatchprefix}.gp - | cut -d ' ' -f 3 | sort -u \
                            > ${tmpvardups}
                            echo -n "$( wc -l ${tmpvardups} | cut -d ' ' -f 1 ) "
                            echo "non-coherent duplicate variants marked for deletion in chromosome ${chr}."
                            cat ${tmpvardups} >> ${batchexcludefile}
                            sort -t ' ' -k 2,2 ${chrbatchprefix}.tped | awk '{
                                delete catalog
                                for ( k = 5; k <= NF; k++ ) {
                                    if ( $k > 0 ) catalog[$k] = 1
                                }
                                asorti( catalog )
                                printf( "%s %s %s", $2, $3, $4 )
                                for ( allele in catalog ) printf( "_%s", catalog[allele] )
                                print( "" )
                            }' | join -t ' ' -v1 - ${tmpvardups} | sort -t ' ' -k 3,3 > ${chrbatchprefix}.gpz
                            echo "extracting coherent variants from the set of duplicate variants.."
                            cut -d ' ' -f 3 ${chrbatchprefix}.gpz | uniq -d \
                            | join -t ' ' -2 3 - ${chrbatchprefix}.gpz | sort -t ' ' -k 3,3r \
                            | sort -t ' ' -u -k 1,1 | cut -d ' ' -f 2 | sort -u > ${tmpvardups}
                            echo -n "$( wc -l ${tmpvardups} | cut -d ' ' -f 1 ) "
                            echo "unique coherent duplicate variants retained in chromosome ${chr}."
                            cut -d ' ' -f 3 ${chrbatchprefix}.gpz | uniq -d \
                            | join -t ' ' -2 3 - ${chrbatchprefix}.gpz | cut -d ' ' -f 2 | sort \
                            | join -t ' ' -v1 - ${tmpvardups} >> ${batchexcludefile}
                            rm ${tmpvardups}
                        done
                    fi
                    touch ${batchflipfile} ${batchexcludefile}
                    if [[ "${refalleles}" != "" && -s "${tmpbatch}.bim" ]] ; then
                        echo "matching variants to reference.."
                        awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${tmpbatch}.bim | sort -t $'\t' -k 7,7 \
                        | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' ${refalleles} - | awk -F $'\t' \
                        $AWKLOCINCLUDE -v batchexcludefile=${batchexcludefile} -v batchflipfile=${batchflipfile} \
                        --lint=true --source 'BEGIN{
                            total_miss = 0
                            total_mism = 0
                            total_flip = 0
                        } {
                            if ( $5 == "-" || $6 == "-" ) {
                                print( $1 ) >>batchexcludefile
                                total_miss++
                            }
                            else {
                                if ( !gmatchx( $2, $3, $5, $6 ) ) {
                                    print( $5 ) >>batchexcludefile
                                    total_mism++
                                }
                                else if ( gflip( $2, $3, $5, $6 ) ) {
                                    print( $5 ) >>batchflipfile
                                    total_flip++
                                }
                            }
                        } END{
                            print( "total missing: ", total_miss )
                            print( "total mismatch: ", total_mism )
                            print( "total flipped: ", total_flip )
                        }'
                    fi
                    touch ${batchctrl}
                fi
                plinkflags=${plinkflagsdef}
                sort -u ${batchexcludefile} > $tmpfile
                mv $tmpfile ${batchexcludefile}
                sort -u ${batchflipfile} > $tmpfile
                mv $tmpfile ${batchflipfile}
                echo "$( wc -l ${batchexcludefile} ) variants to be excluded."
                echo "$( wc -l ${batchflipfile} ) variants to be flipped."
                if [[ -s "${batchexcludefile}" ]] ; then
                    plinkflags="${plinkflags} --exclude ${batchexcludefile}"
                fi
                if [[ -s "${batchflipfile}" ]] ; then
                    plinkflags="${plinkflags} --flip ${batchflipfile}"
                fi
                if [[ "${plinkflags}" != "" ]] ; then
                    mytag=""
                    echo "polishing batch ${batch}.."
                    while [[ -f "${tmpbatch}.bed" ]] ; do
                        tmpbatch=${procdir}/$( basename ${batch} )${mytag}
                        mytag="${mytag}Z"
                    done
                    plink --bfile ${batch} ${plinkflags} --make-bed --out ${tmpbatch}
                    if [[ -f "${tmpbatch}.bim" ]] ; then
                        perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.bim
                    fi
                    if [[ -f "${tmpbatch}.fam" ]] ; then
                        perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.fam
                    fi
                    tmpbatchctrl=${procdir}/$( basename ${tmpbatch} ).fatto
                    if [[ -f "${batchctrl}" && "${batchctrl}" != "${tmpbatchctrl}" ]] ; then
                        cp ${batchctrl} ${tmpbatchctrl}
                    fi
                fi
                if [[ -f "${tmpbatch}.bim" ]] ; then
                    perl -p -e 's/[ \t]+/\t/g' ${tmpbatch}.bim | awk -F $'\t' '{
                        OFS="\t";
                        a[1] = $5; a[2] = $6; asort(a);
                        $2 = $1":"$4"_"a[1]"_"a[2];
                        print;
                    }' ${tmpbatch}.bim > ${tmpfile}
                    mv ${tmpfile} ${tmpbatch}.bim
                fi
                if [[ -f "${tmpbatch}.bed" ]] ; then
                    myarchvec[$j]=${bedhex}
                    mybatchvec[$j]=${tmpbatch}
                    echo "${tmpbatch}" >> ${tmpbatchfile}
                fi
            done
            if [[ -s "${tmpbatchfile}" ]] ; then
                mv ${tmpbatchfile} ${mybatchfile}
            elif [[ -f "${tmpbatchfile}" ]] ; then
                rm ${tmpbatchfile}
            fi
            if (( ${#mybatchvec[*]} > 1 )) ; then
                misfile=${outprefix}.missnp
                mischrfile=${outprefix}.mischr
                plink --merge-list ${mybatchfile} --out ${outprefix}
                if [[ -f "${outprefix}.bim" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
                fi
                if [[ -f "${outprefix}.fam" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
                fi
                egrep '^Warning: Multiple' ${outprefix}.log | cut -d ' ' -f 7 | tr -d "'." \
                | sort -u >> ${mischrfile}
                if [[ -s "${mischrfile}" ]] ; then
                    sort -u ${mischrfile} >> ${misfile}
                fi
                if [[ -f "${misfile}" ]] ; then
                    mv ${misfile} ${excludefile}
                    goflag=""
                fi
            else
                plink --bfile $( cat ${mybatchfile} ) --make-bed --out ${outprefix}
                if [[ -f "${outprefix}.bim" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
                fi
                if [[ -f "${outprefix}.fam" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
                fi
            fi
        done
    fi
    myprefix=${outprefix}

    excludeopt=""
    if [[ "${blacklist}" != "" ]] ; then
        excludeopt="--exclude range ${blacklist}"
    fi
    hqprefix=${myprefix}_hq
    hqprefixunique=${hqprefix}_unique
    cleanfile=${hqprefix}.clean.id
    uniquefile=${hqprefixunique}.rel.id
    prunehqprefix=${hqprefix}_LDpruned
    if [[ ! -f "${prunehqprefix}.bed" || ! -f "${prunehqprefix}.bim" || ! -f "${prunehqprefix}.fam" ]] ; then
        plink --bfile ${myprefix} ${excludeopt} --geno ${varmiss} --maf ${myfreq_hq} --make-bed --out ${hqprefix}
        if [[ -s "${hqprefix}.bed" ]] ; then
            plink --bfile ${hqprefix} --indep-pairphase 500 5 0.2 --out ${hqprefix}_LD
            plink --bfile ${hqprefix} --extract ${hqprefix}_LD.prune.in --make-bed --out ${prunehqprefix}
        fi
    fi
    usex=''
    checkprefix=${prunehqprefix}
    parcount=$( awk '$1 == 25' ${checkprefix}.bim | wc -l )
    if (( parcount == 0 )) ; then
        plink --bfile ${prunehqprefix} --split-x ${genome} no-fail --make-bed --out ${prunehqprefix}_splitx
        checkprefix=${prunehqprefix}_splitx
    fi
    nosex=''
    xcount=$( awk '$1 == 23' ${checkprefix}.bim | wc -l )
    if (( xcount > mincnt )) ; then
        plink --bfile ${checkprefix} --impute-sex --make-bed --out ${prunehqprefix}_isex
        if [[ -s "${prunehqprefix}_isex.nosex" ]] ; then
            nosex="--remove ${prunehqprefix}_isex.nosex"
        fi
        if [[ -s "${prunehqprefix}_isex.fam" ]] ; then
            plink --bfile ${prunehqprefix}_isex ${nosex} --hwe 1.E-${hweneglogp} ${hwtag} \
                --make-just-bim --out ${prunehqprefix}_hwesex
            if [[ -s "${prunehqprefix}_hwsex.bim" ]] ; then
                plink --bfile ${prunehqprefix}_isex --extract <( cut -f2 ${prunehqprefix}_hwsex.bim ) \
                    --impute-sex --make-bed --out ${prunehqprefix}_isex
            fi
            nosex=''
            usex="--update-sex ${prunehqprefix}_isex.fam 3"
            awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' \
                ${prunehqprefix}_isex.fam > ${prunehqprefix}_isex.nosex
            if [[ -s "${prunehqprefix}_isex.nosex" ]] ; then
                nosex="--remove ${prunehqprefix}_isex.nosex"
            fi
        fi
    fi
    keeptag=''
    echo "checking sample quality.."
    if [[ ! -s "${uniquefile}" ]] ; then
        cp ${prunehqprefix}.fam ${cleanfile}
        if [[ "${hvm}" == "on" ]] ; then
            echo "computing individual heterozygosity and missing rates.."
            plink --bfile ${prunehqprefix} --het --out ${hqprefix}_sq
            plink --bfile ${prunehqprefix} --missing --out ${hqprefix}_sq
            ${hvmscript} --args -m ${hqprefix}_sq.imiss -h ${hqprefix}_sq.het -o ${hqprefix}
        fi
        if [[ -s "${cleanfile}" ]] ; then
            echo "identifying duplicate individuals.."
            if [[ "${keepdups}" == "on" ]] ; then pihat=1 ; fi
            plink --bfile ${prunehqprefix} --keep ${cleanfile} --genome gz --out ${myprefix}
            plink --bfile ${prunehqprefix} --keep ${cleanfile} --cluster --read-genome ${myprefix}.genome.gz \
                --rel-cutoff ${pihat} --out ${hqprefixunique}
        fi
    fi
    if [[ -s "${uniquefile}" ]] ; then
        keeptag="--keep ${uniquefile}"
    fi
    outprefix=${myprefix}_clean
    if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
        plink --bfile ${myprefix} ${usex} ${keeptag} --make-bed --out ${outprefix}
        if [[ -f "${outprefix}.bim" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
        fi
        if [[ -f "${outprefix}.fam" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
        fi
    fi
    myprefix=${outprefix}

    outprefix=${myprefix}_varQC
    if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
        tmpvarmiss=${varmiss}
        N=$( wc -l ${myprefix}.fam | cut -d ' ' -f 1 )
        if (( N < mincnt )) ; then tmpvarmiss=0.1 ; fi
        hweneglogp_sex=$(( hweneglogp*2 ))
        if (( hweneglogp_sex > 12 )) ; then hweneglogp_sex=12 ; fi
        plink --bfile ${myprefix} ${nosex} --not-chr 23,24 --geno ${tmpvarmiss} --maf ${myfreq_std} \
            --hwe 1.E-${hweneglogp} ${hwtag} --make-just-bim --out ${outprefix}_nonsex
        plink --bfile ${myprefix} ${nosex} --chr 23,24 --geno ${tmpvarmiss} --maf ${myfreq_std} \
            --hwe 1.E-${hweneglogp_sex} ${hwtag} --make-just-bim --out ${outprefix}_sex
        if [[ -s "${outprefix}_nonsex.bim" || -s "${outprefix}_sex.bim" ]] ; then
            cut -f 2 ${outprefix}_*sex.bim | sort -u > ${outprefix}.mrk
        fi
        if [[ -s "${outprefix}.mrk" ]] ; then
            plink --bfile ${myprefix} --extract ${outprefix}.mrk --make-bed --out ${outprefix}
            if [[ -f "${outprefix}.bim" ]] ; then
                perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
            fi
            if [[ -f "${outprefix}.fam" ]] ; then
                perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
            fi
        fi
    fi
    myprefix=${outprefix}

    outprefix=${myprefix}_sampleQC
    if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
        tmpsamplemiss=${samplemiss}
        N=$( wc -l ${myprefix}.bim | cut -d ' ' -f 1 )
        if (( N < mincnt )) ; then tmpsamplemiss=0.1 ; fi
        plink --bfile ${myprefix} --mind ${tmpsamplemiss} --make-bed --out ${outprefix}
        if [[ -f "${outprefix}.bim" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
        fi
        if [[ -f "${outprefix}.fam" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
        fi
    fi
    myprefix=${outprefix}

    if [[ "${phenotypes}" != "" && -f "${phenotypes}" ]] ; then
        outprefix=${myprefix}_pheno
        if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
            plink --bfile ${myprefix} --make-pheno ${phenotypes} affected --make-bed \
                --out ${outprefix}
            if [[ -f "${outprefix}.bim" ]] ; then
                perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
            fi
            if [[ -f "${outprefix}.fam" ]] ; then
                perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
            fi
        fi
        myprefix=${outprefix}
        outprefix=${myprefix}_hwe_ctrl
        if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
            hweneglogp_ctrl_sex=$(( hweneglogp_ctrl/2 ))
            if (( hweneglogp_ctrl_sex < 12 )) ; then hweneglogp_ctrl_sex=12 ; fi
            plink --bfile ${myprefix} ${nosex} --not-chr 23,24 --hwe 1.E-${hweneglogp_ctrl} \
                --make-just-bim --out ${outprefix}_nonsex
            plink --bfile ${myprefix} ${nosex} --chr 23,24 --hwe 1.E-${hweneglogp_ctrl_sex} \
                --make-just-bim --out ${outprefix}_sex
            if [[ -s "${outprefix}_nonsex.bim" || -s "${outprefix}_sex.bim" ]] ; then
                cut -f 2 ${outprefix}_*sex.bim | sort -u > ${outprefix}.mrk
            fi
            if [[ -s "${outprefix}.mrk" ]] ; then
                plink --bfile ${myprefix} --extract ${outprefix}.mrk --make-bed --out ${outprefix}
                if [[ -f "${outprefix}.bim" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
                fi
                if [[ -f "${outprefix}.fam" ]] ; then
                    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
                fi
            fi
        fi
        myprefix=${outprefix}
        mycontrols=${mydir}/$( basename ${phenotypes} .txt )_ctrl.txt
        mycontrols_with_batch=${mydir}/$( basename ${phenotypes} .txt )_ctrl.batch
        awk -F $'\t' '{ if ( tolower($3) == "control" ) print; }' ${phenotypes} > ${mycontrols}
    else
        mycontrols=${myprefix}_ctrl.txt
        mycontrols_with_batch=${myprefix}_ctrl.batch
        cut -f 1,2,6 ${myprefix}.fam | perl -p -e 's/[ \t]+/\t/g' > ${mycontrols}
    fi

    if (( ${#mybatchvec[*]} > 1 )) ; then
        tmpfile=$( mktemp .tmpbcXXXXX )
        echo "created temporary batch control file ${tmpfile}."
        echo -n > ${mycontrols_with_batch}
        for batch in $( cat ${mybatchfile} ) ; do
            tmpbatch=${mydir}/$( basename ${batch} )_ctrl
            plink --bfile ${batch} --keep ${mycontrols} --make-bed --out ${tmpbatch}
            if [[ -f "${tmpbatch}.bim" && -f "${tmpbatch}.fam" ]] ; then
                perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.bim ${tmpbatch}.fam
                awk -F $'\t' -v batch=$( basename ${batch} ) '{ OFS="\t"; print( $1, $2, batch ) }' \
                    ${tmpbatch}.fam >> ${mycontrols_with_batch}
            fi
        done
        sort -u -k 1,2 ${mycontrols_with_batch} > ${tmpfile}
        mv ${tmpfile} ${mycontrols_with_batch}
        tmpfile=$( mktemp .tmpmodXXXXX )
        echo "created temporary model file ${tmpfile}."
        bevarfile=${myprefix}.exclude
        echo -n > ${bevarfile}
        for bbatch in $( cat ${mybatchfile} ) ; do
            batch=$( basename ${bbatch} )
            echo "assessing batch effects for '${batch}'.."
            plinkfile=${mydir}/plink_${batch}
            plink --bfile ${myprefix} --make-pheno ${mycontrols_with_batch} ${batch} --model \
                --out ${plinkfile}
            if [[ -f "${plinkfile}.model" ]] ; then
                perl -p -e 's/^[ \t]+//g;' ${plinkfile}.model | perl -p -e 's/[ \t]+/\t/g' > ${tmpfile}
                mv ${tmpfile} ${plinkfile}.model
                for atest in $( cut -f 5 ${plinkfile}.model | tail -n +2 | sort -u ) ; do
                    echo "summarizing ${atest} tests.."
                    awk -F $'\t' -v atest=${atest} 'NR > 1 && $5 == atest' ${plinkfile}.model \
                    | cut -f 2,10 | ${fdrscript} > ${tmpfile}
                    awk -F $'\t' '{ if ( $2 < 0.9 ) print( $1 ); }' ${tmpfile} >> ${bevarfile}
                done
            fi
        done
        sort -u ${bevarfile} > ${tmpfile}
        mv ${tmpfile} ${bevarfile}
        outprefix=${myprefix}_nbe
        plink --bfile ${myprefix} --exclude ${bevarfile} --make-bed --out ${outprefix}
        if [[ -f "${outprefix}.bim" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
        fi
        if [[ -f "${outprefix}.fam" ]] ; then
            perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
        fi
    fi
    if [[ -f "${outprefix}.fam" ]] ; then
        cp ${outprefix}.fam ${outprefix}.fam.org
        awk '{ OFS="\t"; for ( k = 1; k < 5; k++ ) gsub( "[/+-]", "_", $k ); print; }' ${outprefix}.fam.org \
        > ${outprefix}.fam
    fi

    hqprefix=${outprefix}_hq
    plink --bfile ${outprefix} ${excludeopt} --maf ${myfreq_hq} --make-bed --out ${hqprefix}
    plink --bfile ${hqprefix} --indep-pairphase 500 5 0.2 --out ${hqprefix}_LD
    plink --bfile ${hqprefix} --extract ${hqprefix}_LD.prune.in --make-bed --out ${hqprefix}_LDpruned
    plink --bfile ${hqprefix}_LDpruned --genome gz --out ${hqprefix}
    plink --bfile ${hqprefix}_LDpruned --cluster --read-genome ${hqprefix}.genome.gz --pca header tabs var-wts \
        --out ${hqprefix}

    echo -e "\nall done. check your output files out."
    echo -e "\n================================================================================\n"

fi

