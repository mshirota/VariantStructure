grep -v '>' $1 | perl -pe "s/[\n\r]//g" | perl -pe "s/(.{10})/\$1 /g" | perl -pe "s/$/\n/" | perl -pe "s/(.{55})/\$1\n/g" | awk '{printf "%5d\t%s\n", (NR-1)*50+1, $0}'
