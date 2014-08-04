#!/bin/bash
s1=t1.00
s2=r0.20
#s1=r0.10
#s2=r0.90
unset PRUNING_SINGLE_OUTGROUP ;
function runsim {
    maxr=$1
    tag="$2"
    cfg="$3"
    rm "grand_summary_${tag}.txt"
    for ((i=0 ; i < $maxr ; ++i))
    do
        unset REVERSE_MAPPING ;
        python ../summary_stats.py "${tag}${i}" > "summary_${tag}${i}" || exit
        export REVERSE_MAPPING=1 ;
        python ../summary_stats.py "${tag}${i}" > "rev_summary_${tag}${i}" || exit
        sv1=$(grep ^${s1} "summary_${tag}${i}" | awk '{print $3}')
        sv2=$(grep ^${s2} "summary_${tag}${i}" | awk '{print $3}')
        rsv1=$(grep ^${s1} "rev_summary_${tag}${i}" | awk '{print $3}')
        rsv2=$(grep ^${s2} "rev_summary_${tag}${i}" | awk '{print $3}')
        echo $i ${sv1} ${sv2} ${rsv1} ${rsv2}>> "grand_summary_${tag}.txt"
    done
}

nreps=10
runsim $nreps c "co-speciation.cfg"
runsim $nreps n "non-phylo-jump.cfg"
runsim $nreps g "geo-restricted.cfg"
runsim $nreps r "related-hosts.cfg"

