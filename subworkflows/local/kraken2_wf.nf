include { KRAKEN2      } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN_PARSE } from '../../modules/local/kraken_parse'
include { KRAKEN_MERGE } from '../../modules/local/kraken_merge'

workflow KRAKEN2_WF {
    take:
        path reads
        path kraken2_db
    main:
        ch_versions = Channel.empty()
        KRAKEN2(reads, kraken2_db)
        KRAKEN_PARSE(KRAKEN2.out.report)
        KRAKEN_MERGE(KRAKEN_PARSE.out.kraken_read_count.collect())
        ch_versions = ch_versions.mix(KRAKEN_MERGE.out.versions.first())

    emit:
        kraken_merged_report: KRAKEN_MERGE.out.kraken_merged_report
        versions: ch_versions
}
