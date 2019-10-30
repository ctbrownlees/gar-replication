for file in $(ls *.R)
  do
   R --slave --no-restore --file=$file
  done

