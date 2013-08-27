fn=$1
pt=$2
cat ConvertAdapterTemplate.h.in | sed -e "s/%function%/$fn/" -e "s/%param%/$pt/" | tee $fn.h
cat ConvertAdapterTemplate.cxx.in | sed -e "s/%fn%/$fn/" -e "s/%param%/$pt/" | tee $fn.cxx
