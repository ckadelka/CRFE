- Map between GO_ID and term_name: from go.obo (04Dec22 archive)
grep -e '^id:' -e '^name:' ~/RotationData/Friedberg/gaf_files/OBO_Archives/22go.obo | sed 's/^id: \|^name: //g' | paste - - -d'\t' > GO_human_name_ID.txt

- Get GO annotations (BPO only): from .gaf (04Dec22 archive)
awk '$1 !~ /\!/' ../gaf_files/All_Annotations_Archives/04Dec22goa_human.gaf | cut -f 3,5,9 | awk '$3 ~ /P/ {print $1,$2}' | sort | uniq | less > ../../Kadelka/CRFE/gith/GOhuman_04Dec22.txt

- Group genes by term to get annotation file (ID) for CRFE: 
awk '
        NR == 1 {
                print $2 "\t" $1
        }
        NR > 1 {
                A[$2] = A[$2] ? A[$2] OFS $1: $1
        }
        END {
                for ( k in A )
                        print k "\t" A[k]
        }
' GOhuman_04Dec22.txt > GOhuman_ID.txt

- Full join by GO_ID to get annotation file (named) for CRFE:
join -j 1 -o 1.1,1.2,2.2 <(sort -k1 GO_human_name_ID.txt) <(sort -k1 data/GOhuman_ID.txt) -t $'\t' | cut -f 2,3 > data/GOhuman_named.txt

