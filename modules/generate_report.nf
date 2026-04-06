process GENERATE_REPORT {
    tag "report_generation"
    label "report_generation"

    publishDir params.outdir, mode: 'copy'

    input:
    path(gff3_annotations)

    output:
    path("pannotator_report.html")

    script:
    """
    generate_report.py \
        --gff3_files ${gff3_annotations.join(' ')} \
        --out pannotator_report.html
    """
}
