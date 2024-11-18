process MODULE_NAME {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container ''

    input:
    tuple val(meta)

    output:
    tuple val(meta)

    script:
    """
    """
}