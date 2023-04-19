include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'

workflow PREPARE_GENOMES {
    take:
        genomes  // [meta["genome_name": genome_name], val(igenome), file(row.fasta), file(row.index)]

    main:
        ch_versions = Channel.empty()

        genomes
            .branch {
                igenome: it[1].toString().toUpperCase() != 'NA'
                index_avail: it[3].toString().toUpperCase() != 'NA' && it[1].toString().toUpperCase() == 'NA'
                no_index_avail: it[3].toString().toUpperCase() == 'NA' && it[1].toString().toUpperCase() == 'NA'
            }
            .set { genomes_fork }

        BOWTIE2_BUILD (
            genomes_fork.no_index_avail.map {
                meta, igenome, fasta, index -> [meta, path(fasta)]
            }
        )

        genomes_fork.no_index_avail.map {
                meta, igenome, fasta, index -> [meta, fasta]
        }.join(
            BOWTIE2_BUILD.out.index
        ).mix(
            genomes_fork.index_avail.map {
                meta, igenome, fasta, index -> [meta, path(fasta), path(index)]
            },
            genomes_fork.igenome.map {
                meta, igenome, fasta, index -> [meta, path(params.genomes[ igenome ][ fasta ]), path(params.genomes[ igenome ][ bowtie2 ])]
            }
        ).set (ch_genomes)

        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    emit:
        genomes = ch_genomes // [meta["genome_name": genome_name], path(fasta), path(index)]
        versions = ch_versions
}
