fn=${1}
chain=${2}
resid=${3}
chain_resid=$(printf "%s%4d" ${2} ${3})
awk '$0 ~ /^ATOM  .{7}CA.{6}'"$chain_resid"'/' $fn | head -n 1 | cut -b 61-66

