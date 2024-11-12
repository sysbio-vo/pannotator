#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

process hello {

    output:
        stdout

    script:
    """
    echo 'I do desire that we may be better strangers.'
    """
}

workflow {

    hello()
}