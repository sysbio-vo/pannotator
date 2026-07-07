process DOWNLOAD_BAKTA_DB {
    cache false
    tag "setup"
    label "setup"
    label 'bakta'
    publishDir params.bakta_db_dir, mode: 'copy'

    input:
    val(bakta_db_type)

    output:
    path(relative_db_path)

    script:
    if ( bakta_db_type == 'light' ) {
        relative_db_path = "db-light"
    } else if ( bakta_db_type == 'full' ) {
        relative_db_path = "db"
    }
    """
    bakta_db download --output . --type ${bakta_db_type}
    """
}
