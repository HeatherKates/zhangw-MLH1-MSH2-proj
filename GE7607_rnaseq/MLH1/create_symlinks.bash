# Recreate symlinks in ../data/ to point to the original files
find /orange/zhangw/GE7607/GE-7607-WZhang-10B-22MNM5LT3-Lane2 -type f -name "*MLH1*.fastq.gz" | while read -r filepath; do
    filename=$(basename "$filepath")
    ln -sf "$filepath" "data/$filename"
done

