//
// Check input samplesheet and get read channels
//
def sampleIdFromName = {name -> name.replaceFirst(~/(\.[^\.]+)+$/, '')}


workflow DIR_PARSE {
    take: 
    indir // indir paramter 

    main:
    Channel.fromPath("${params.indir}/*${params.infile_extension}")
        .map { asm -> tuple([id: sampleIdFromName(asm.name)], asm) }
        .set{ infiles }

    emit:
    infiles

}

workflow MANIFEST_PARSE {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty { error "Cannot find path file ${samplesheet}" }
        .splitCsv ( header:true, sep:',' )
        .flatMap { row -> create_assembly_channels(row) }
        .set{ infiles }

    emit:
    infiles // tuple([id: sample_id], assembly_file), same as DIR_PARSE
}

// Function to get list of tuple([ id: ... ], assembly_file)
// 'assembly' can be a file, or a directory of assemblies (ids then come from the file names, as in DIR_PARSE)
def create_assembly_channels(Map row) {
    if ( !row.ID || !row.assembly ) {
        error "ERROR: Please check input samplesheet -> 'ID' and 'assembly' fields are required: ${row}"
    }

    def path = file(row.assembly)
    if ( !path.exists() ) {
        error "ERROR: Please check input samplesheet -> Assembly file does not exist!\n${row.assembly}"
    }

    if ( path.isDirectory() ) {
            error "ERROR: Please check input samples sheet -> Assembly entry is a directory ${row.assembly}"
        }

    return [ tuple([id: row.ID], path) ]
    }
