
acc=SRS6067521
# Get headers
esearch \
    -db sra \
    -query ERS3350306 \
    | efetch \
    -format runinfo \
    |head -n1 \
> test.txt

#    acc=SRS6151291
    echo $acc

    esearch \
        -db sra \
        -query $acc \
        | efetch \
        -format runinfo \
        |tail -n +2 \
        >> test.txt
