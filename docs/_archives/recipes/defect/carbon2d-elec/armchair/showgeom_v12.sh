for dir in vacancy[12]_density; do
  echo "Converting geo in $dir"
  cd $dir
  gen2xyz geo.gen
  cd -
done

for dir in vacancy[12]_density; do
  jmol -L $dir/geo.xyz &
done

