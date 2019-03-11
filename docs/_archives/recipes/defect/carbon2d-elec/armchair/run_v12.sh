for dir in vacancy[12]_density; do
  echo "calculating in $dir"
  cd $dir
  ./run.sh
  cd -
done
