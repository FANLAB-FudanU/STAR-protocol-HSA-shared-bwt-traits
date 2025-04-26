#!/usr/bin/bash
if [[ $# -eq 0 ]]; then
    echo "Usage: cal_used_time.sh starttime endtime"
    echo 'Example: cal_used_time.sh "$starttime" "$endtime"'
    echo "Warning: Must use double quote to avoid space splitting !!!!!"
    exit 
fi

starttime=$1
endtime=$2

# echo "starttime is $starttime"
# echo "endtime is $endtime"

start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数
end_seconds=$(date --date="$endtime" +%s);

echo "$start_seconds" "$end_seconds" | awk '{total_sec=$2-$1; min=int(total_sec/60); lef_sec=total_sec%60; printf("Total time: %d minutes, %d seconds.\n", min, lef_sec); }'

