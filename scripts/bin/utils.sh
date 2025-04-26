#!/usr/bin/bash 

function loginfo() {
    time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[ INFO $time ] [ $1 ]"
}

function logerror() {
    time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[ ERROR $time ] [ $1 ]"
    return 1
}


function get_chrom_idx()
{
    if [ $# -eq 0 ]; then
        echo "Usage: "
        echo 'get_chrom_idx $loci_file'
        return 0
    fi

    loci_file=$1
    chrom_idx=$(head -1 "$loci_file" |sed 's/\t/\n/g' |awk '$1=="chr" || $1=="chrom" || $1=="CHR" || $1=="CHROM" || $1=="lead_chrom"{print NR}')
    if [[ -z $chrom_idx ]]; then
        return 1
    fi
    echo "$chrom_idx"
    return 0
}

function get_start_idx()
{
    if [ $# -eq 0 ]; then
        echo "Usage: "
        echo 'get_start_idx $loci_file'
        return 0
    fi

    loci_file=$1
    start_idx=$(head -1 "$loci_file" |sed 's/\t/\n/g' |awk '$1=="locus_start" || $1=="start"{print NR}')
    if [[ -z $start_idx ]]; then
        return 1
    fi
    echo "$start_idx"
    return 0
}

function get_stop_idx()
{
    if [ $# -eq 0 ]; then
        echo "Usage: "
        echo 'get_stop_idx $loci_file'
        return 0
    fi

    loci_file=$1
    stop_idx=$(head -1 "$loci_file" |sed 's/\t/\n/g' |awk '$1=="locus_stop" || $1=="stop"{print NR}')
    if [[ -z $stop_idx ]]; then
        return 1
    fi
    echo "$stop_idx"
    return 0
}

function get_ld_prefix_idx()
{
    if [ $# -eq 0 ]; then
        echo "Usage: "
        echo 'get_ld_prefix_idx $loci_file'
        return 0
    fi

    loci_file=$1
    ld_prefix_idx=$(head -1 "$loci_file" |sed 's/\t/\n/g' |awk '$1=="ld_prefix" || $1=="LD_prefix"{print NR}')
    if [[ -z $ld_prefix_idx ]]; then
        return 1
    fi
    echo "$ld_prefix_idx"
    return 0
}

function get_locus_info_idxs()
{
    if [ $# -eq 0 ]; then
        echo "Usage: "
        echo 'get_locus_info_idxs $loci_file $return_ld_prefix'
        echo 'read -r chrom_idx start_idx stop_idx ld_prefix_idx <<< "$(get_locus_info_idxs "$loci_file" True)"'
        echo 'read -r chrom_idx start_idx stop_idx <<< "$(get_locus_info_idxs "$loci_file" False)"'
        return 0
    fi

    loci_file=$1
    return_ld_prefix=$2
    # get chrom_idx
    if ! chrom_idx=$(get_chrom_idx "$loci_file"); then
        echo "Error! $loci_file does not contain chrom column"
        return 1
    fi

    # get start_idx
    if ! start_idx=$(get_start_idx "$loci_file"); then
        echo "Error! $loci_file does not contain start column"
        return 1
    fi

    # get stop_idx
    if ! stop_idx=$(get_stop_idx "$loci_file"); then
        echo "Error! $loci_file does not contain stop column"
        return 1
    fi

    if [[ $return_ld_prefix == "True" ]]; then
        # get ld_prefix_idx
        if ! ld_prefix_idx=$(get_ld_prefix_idx "$loci_file"); then
            echo "Error! $loci_file does not contain ld_prefix column"
            return 1
        fi

        echo -e "$chrom_idx\t$start_idx\t$stop_idx\t$ld_prefix_idx"

    else
        echo -e "$chrom_idx\t$start_idx\t$stop_idx"
        
    fi
    return 0
}
