python /maintenance/liftOver_catalog.py \
    --species RatNor \
    --map Littrelletal2018_rn6 \
    --chainFile rn6ToRn7.over.chain.gz \
    --validationChain rn7ToRn6.over.chain.gz \
    --winLen 1000 \
    --useAdjacentAvg \
    --retainIntermediates \
    --gapThresh 1000000
