for i in $(find . -name "*.nim");
do
   if [[ -f "${i%.nim}" ]];
   then 
     echo removing ${i%.nim}
     rm ${i%.nim}
   fi
done

