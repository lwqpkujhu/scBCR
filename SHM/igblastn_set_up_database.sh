## set up igblast
mkdir /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast/database
/home/weiqiang/software/kleinstein-immcantation-73434b4fb26c/scripts/fetch_igblastdb.sh -o /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast
cp -r /home/weiqiang/software/ncbi-igblast-1.14.0/internal_data /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast
cp -r /home/weiqiang/software/ncbi-igblast-1.14.0/optional_file /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast

/home/weiqiang/software/kleinstein-immcantation-73434b4fb26c/scripts/fetch_imgtdb.sh -o /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt
/home/weiqiang/software/kleinstein-immcantation-73434b4fb26c/scripts/imgt2igblast.sh -i /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt -o /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast

