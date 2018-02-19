cut -d',' -f8 curgen_db_001.txt | sort | uniq > UnValues
grep -c "." UnValues
