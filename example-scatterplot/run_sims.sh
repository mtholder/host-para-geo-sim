#!/bin/bash
s1=t1.00
s2=r0.20
#s1=r0.10
#s2=r0.90

function runsim {
    maxr=$1
    tag="$2"
    cfg="$3"
    rm "grand_summary_${tag}.txt"
    for ((i=0 ; i < $maxr ; ++i))
    do
        python ../sim-host-parasite.py $RANDOM $RANDOM 50 40 50 "${cfg}"  >"${tag}${i}" || exit
        python ../summary_stats.py "${tag}${i}" > "summary_${tag}${i}" || exit
        sv1=$(grep ^${s1} "summary_${tag}${i}" | awk '{print $3}')
        sv2=$(grep ^${s2} "summary_${tag}${i}" | awk '{print $3}')
        echo $i ${sv1} ${sv2} >> "grand_summary_${tag}.txt"
    done
}

nreps=10
runsim $nreps c "co-speciation.cfg"
runsim $nreps n "non-phylo-jump.cfg"
runsim $nreps g "geo-restricted.cfg"
runsim $nreps r "related-hosts.cfg"

