for file in $(find /orange/zhangw/NS3842/NS-3842-WZhang-15B-223353LT1-Lane1-2 -type f -name "*MLH1*"); do
    ln -s $file data/
done

