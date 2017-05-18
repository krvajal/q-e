#!/bin/bash
elem="h he li be b c n o f ne na mg"
for i in $elem
do 
    echo "$i"
    ../src/ld1.x < ../examples/all-electron/kli/$i.kli.in > $i.kli.out
    grep  "Etot" $i.kli.out
    cp "kli_pot.up" $i"_kli_pot.up"
    cp "kli_pot.dw" $i"_kli_pot.dw"
    paste "kli_pot.up" "kli_pot.dw" > $i."_kli_pot.ud"
done