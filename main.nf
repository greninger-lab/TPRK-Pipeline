nextflow.enable.dsl=1

def helpMessage() {
  log.info"""
  tprK pipeline, with nextflow!

  An example command for running the pipeline is as follows:
  nextflow run michellejlin/tprk -r nextflow --INPUT ./ --OUTDIR output/ -resume -with-docker ubuntu:18.04 -with-trace \\

  Mandatory Arguments
  --INPUT         Input folder where all fastqs are located.
  ./ can be used for current directory.
  --OUTDIR        Output directory.
  --METADATA      Metadata file formatted in a .csv with columns: SampleName, Illumina, PacBio.
  If running with --PACBIO or --ILLUMINA simply leave those columns blank (but make sure to
    have commas as appropriate). Names should be in PB_<sample_name>.fastq or Ill_<sample_name>.fastq format.
    See the example metadatas in the example/ folder for a sample.

    Input Specifications
    --PACBIO        Write this flag to specify that there are only PacBio files here.
    Comparison figures to Illumina will not be generated.
    --ILLUMINA      Write this flag to specify that there are only Illumina files here.
    Comparison figures to PacBio will not be generated.
    --REFERENCE     Specify Illumina sample name (not file), to compare others to for dot-line plots. Can be used in tandem with --LARGE.

    Filtering Options
    --RF_FILTER         Optional flag for specifying what relative frequency
    an additional filtered final merged table and visualizations should be sorted at.
    By default this is set to 0.00001.
    --COUNT_FILTER      Optional flag for specifying what count
    an additional filtered final merged table and visualizations should be sorted at.
    By default this is set to 0.

    """.stripIndent()
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  /*                                                    */
  /*          SET UP CONFIGURATION VARIABLES            */
  /*                                                    */
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  // Show help message
  params.help = false
  if (params.help){
    helpMessage()
    exit 0
  }

  params.INPUT = false
  params.OUTDIR= false
  params.RF_FILTER = 0.00001
  params.COUNT_FILTER = 0
  params.ILLUMINA_FILTER = false
  params.PACBIO = false
  params.ILLUMINA = false
  params.METADATA = false
  params.LARGE = false
  params.REFERENCE = false

  INPUT_TYPE = "both"
  PACBIO_FLAG = ""
  ILLUMINA_FLAG = ""
  INPUT_TYPE2 = ""

  SYPH_R = file("${baseDir}/bin/syph_r.py")
  COMPARE_DF = file("${baseDir}/bin/compare_df.R")
  FILTER_ALL_READS = file("${baseDir}/bin/filterAllReads.py")
  RECALCULATE_FREQUENCY = file("${baseDir}/bin/recalculate_frequency.R")
  RAD_FREQUENCY = file("${baseDir}/bin/RAD_Frequency.R")
  VARIABLE_REGION_COMPARE = file("${baseDir}/bin/Variable_region_compare.R")

  /////////////////////////////
  /*    VALIDATE INPUTS      */
  /////////////////////////////

  // if METADATA not set
  if (params.METADATA == false) {
    println("Must provide metadata file input as .csv format with three columns: \
    SampleName, Illumina, PacBio. Use --METADATA flag.")
    exit(1)
    } else{
      METADATA_FILE = file(params.METADATA)
    }
    // if INPUT not set
    if (params.INPUT == false) {
      println( "Must provide an input directory with --INPUT")
      exit(1)
    }
    // Make sure INPUT ends with trailing slash
    if (!params.INPUT.endsWith("/")){
      params.INPUT = "${params.INPUT}/"
    }
    // if OUTDIR not set
    if (params.OUTDIR == false) {
      println( "Must provide an output directory with --OUTDIR")
      exit(1)
    }
    // Make sure OUTDIR ends with trailing slash
    if (!params.OUTDIR.endsWith("/")){
      params.OUTDIR = "${params.OUTDIR}/"
    }
    // Figure out input of PACBIO or ILLUMINA
    if(((params.PACBIO) && (params.ILLUMINA))){
      println("--PACBIO and --ILLUMINA cannot be used together. Please specify only one, \
      or do not use these flags if you want to run the default way of comparing both PacBio and Illumina files.")
      exit(1)
      } else if (params.PACBIO) {
        INPUT_TYPE = "pacbio"
        PACBIO_FLAG = "--pacbio"
        println("--PACBIO indicated. Will not compare to Illumina files or generate figures.")
        } else if (params.ILLUMINA) {
          INPUT_TYPE = "illumina"
          ILLUMINA_FLAG = "--illumina"
          println("--ILLUMINA indicated. Will not compare to PacBio files or generate figures.")
        }

        // Helpful messages
        println("Will filter final products for >${params.RF_FILTER} relative frequency \
        and >${params.COUNT_FILTER} count.")
        if (params.ILLUMINA_FILTER){
          println("Will only include PacBio reads supported by Illumina reads that pass the filter, \
          as specified by --ILLUMINA_FILTER.")
        }

        /////////////////////////////
        /*    METADATA PARSING     */
        /////////////////////////////

        // Reads in Illumina/PacBio pairs from metadata
        if (INPUT_TYPE == "both") {
          input_pacbio_ch = Channel
          .fromPath(METADATA_FILE)
          .splitCsv(header:true)
          .map{ row-> tuple(row.SampleName, file(row.PacBio)) }
          illumina_ch = Channel
          .fromPath(METADATA_FILE)
          .splitCsv(header:true)
          .map{ row-> tuple(row.SampleName, file(row.Illumina)) }
          metadata_ch = Channel
          .fromPath(METADATA_FILE)
          .splitCsv(header:true)
          .map{ row-> tuple(row.SampleName, file(row.Illumina), file(row.PacBio), INPUT_TYPE) }
          } else if (INPUT_TYPE == "illumina") {
            illumina_ch = Channel
            .fromPath(METADATA_FILE)
            .splitCsv(header:true)
            .map{ row-> tuple(row.SampleName, file(row.Illumina)) }
            metadata_ch = Channel
            .fromPath(METADATA_FILE)
            .splitCsv(header:true)
            .map{ row-> tuple(row.SampleName, file(row.Illumina), file(METADATA_FILE), INPUT_TYPE) }
            } else if (INPUT_TYPE == "pacbio") {
              input_pacbio_ch = Channel
              .fromPath(METADATA_FILE)
              .splitCsv(header:true)
              .map{ row-> tuple(row.SampleName, file(row.PacBio)) }

              metadata_ch = Channel
              .fromPath(METADATA_FILE)
              .splitCsv(header:true)
              .map{ row-> tuple(row.SampleName, file(row.PacBio), INPUT_TYPE) }
            }

            sample_name_ch = Channel
            .fromPath(METADATA_FILE)
            .splitCsv(header:true)
            .map{ row-> row.SampleName }

            ////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////
            /*                                                    */
            /*        CREATE FREQUENCY TABLES PER SAMPLE          */
            /*                                                    */
            ////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////

            //
            // PacBio Section
            //
            if(INPUT_TYPE != "illumina") {
              // Denoises PacBio files with RAD.
              process denoisePacBioFiles {

                container "cave42/denoise_pac"

                // Retry on fail at most three times
                //errorStrategy 'retry'
                //maxRetries 0

                input:
                tuple val(base), file(PACBIO_FILE) from input_pacbio_ch
                file(RAD_FREQUENCY)

                output:
                tuple val(base), file("PB_${base}.noprimers.filtered.RAD.nolines.fix.fasta") into pacbio_ch
                tuple val(base), file("PB_${base}.noprimers.filtered.RAD.nolines.fix.fasta") into pacbio_ch2
                val(base) into pacbio_sample_name_ch
                publishDir "${params.OUTDIR}/denoised_fastas", mode: 'copy', pattern: '*.fix.fasta'

                script:
                """
                gunzip -d --force *.gz

                Rscript ${RAD_FREQUENCY} -s ${baseDir} -d ${params.INPUT} -m ${METADATA_FILE} -a PB_${base}.fastq -c ${task.cpus}
                """
              }

              // Create frequency tables for each PacBio sample.
              process createFrequencyTables_PacBio {
                container "quay.io/greninger-lab/tprk"

                // Retry on fail at most three times
                errorStrategy 'retry'
                maxRetries 0

                publishDir "${params.OUTDIR}/Tables/Frequency_Tables/", mode: 'copy', pattern: '*.csv'

                input:
                file(PACBIO_FILE)from pacbio_ch.collect()
                file(METADATA_FILE)
                file(COMPARE_DF)
                file(SYPH_R)

                output:
                file("*final_data.csv") into pacbio_final_data_ch
                file("*final_data.csv") into final_data_ch_pb
                file("*final_AA_data.csv") into final_dna_data_ch_pacbio
                file "all_assignments.csv" into all_assignments_ch1
                file("allreads.csv") optional true into allreads_ch
                file "compare_pacbio_df.csv" into compare_pacbio_ch
                file("*summary_statistics.csv") into summary_stats_pacbio_ch

                script:
                """
                Rscript ${COMPARE_DF} -s ${SYPH_R} -m ${METADATA_FILE} -d ./ --pacbio -c ${task.cpus}
                """
              }
            }

            //
            // Illumina Section
            //

            // If --ILLUMINA, no all_assignments.csv will have been created in
            // createFrequencyPlots_PacBio. Creates it here.
            if (INPUT_TYPE == "illumina") {
              process createAllAssignments{
                input:
                output:
                file("all_assignments.csv") into all_assignments_ch1
                file("compare_pacbio_df.csv") into compare_pacbio_ch

                script:
                """
                touch all_assignments.csv
                touch compare_pacbio_df.csv
                """
              }
            }

            if (INPUT_TYPE != "pacbio") {
              // Create frequency tables for each Illumina sample.
              // Also grabs frequency tables from PacBio samples and merges the two,
              // creating allreads.csv file.
              process createFrequencyTables_Illumina {

                container "quay.io/greninger-lab/tprk"

                publishDir "${params.OUTDIR}/Tables/", mode: 'copy', pattern: '*.csv'

                // Retry on fail at most three times
                // errorStrategy 'retry'
                // maxRetries 3

                input:
                file(ILLUMINA_FILE) from illumina_ch.collect()
                file "all_assignments.csv" from all_assignments_ch1
                file("compare_pacbio_df.csv") from compare_pacbio_ch
                file(METADATA_FILE)
                file(COMPARE_DF)
                file(SYPH_R)

                output:
                file("*final_data.csv") into illumina_final_data_ch
                file("*final_data.csv") into final_data_ch_ill
                file("*final_AA_data.csv") into final_dna_data_ch
                file("all_assignments.csv") into all_assignments_ch2
                file("allreads.csv") into allreads_ch
                file("*summary_statistics.csv") into summary_stats_ill_ch
                file("*") into illumina_frequency_tables_ch

                script:
                """
                gunzip -d --force *.gz

                Rscript ${COMPARE_DF} -s ${SYPH_R} -m ${METADATA_FILE} -d ./ --illumina -c ${task.cpus}
                """
              }

              process summaryStats_Illumina {

                input:
                file("*summary_statistics.csv") from summary_stats_ill_ch.collect()

                output:
                file("*.csv") into summary_stats_ill_ch2
                publishDir "${params.OUTDIR}/Tables/", mode: 'copy', pattern: 'all_summary_stats.csv'

                shell:
                '''
                echo Sample,Total Input Reads,V1_Reads,V2_Reads,V3_Reads,V4_Reads,V5_Reads,V6_Reads,V7_Reads > all_summary_stats.csv
                for file in *summary_statistics.csv;do tail -n +2 $file>>all_summary_stats.csv;echo "" >> all_summary_stats.csv; done
                inputreads=$(awk -F, '{ total += $2 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V1=$(awk -F, '{ total += $3 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V2=$(awk -F, '{ total += $4 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V3=$(awk -F, '{ total += $5 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V4=$(awk -F, '{ total += $6 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V5=$(awk -F, '{ total += $7 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V6=$(awk -F, '{ total += $8 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)
                V7=$(awk -F, '{ total += $9 } END { print total / (( NR - 1 ))}' all_summary_stats.csv)

                echo "" >> all_summary_stats.csv
                echo "mean:,"$inputreads","$V1","$V2","$V3","$V4","$V5","$V6","$V7"" >> all_summary_stats.csv

                '''
              }
          }

            if (INPUT_TYPE != "illumina") {
              process summaryStats_PacBio {

                input:
                file("*summary_statistics.csv") from summary_stats_pacbio_ch.collect()

                output:
                file("*.csv") into summary_stats_pb_ch2
                publishDir "${params.OUTDIR}/Tables/", mode: 'copy', pattern: 'all_summary_stats.csv'

                shell:
                '''
                echo Sample,Total Input Reads,V1_Reads,V2_Reads,V3_Reads,V4_Reads,V5_Reads,V6_Reads,V7_Reads > all_summary_stats_pacbio.csv
                for file in *summary_statistics.csv;do tail -n +2 $file>>all_summary_stats_pacbio.csv;echo "" >> all_summary_stats_pacbio.csv; done
                inputreads=$(awk -F, '{ total += $2 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V1=$(awk -F, '{ total += $3 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V2=$(awk -F, '{ total += $4 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V3=$(awk -F, '{ total += $5 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V4=$(awk -F, '{ total += $6 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V5=$(awk -F, '{ total += $7 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V6=$(awk -F, '{ total += $8 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)
                V7=$(awk -F, '{ total += $9 } END { print total / (( NR - 1 ))}' all_summary_stats_pacbio.csv)

                echo "" >> all_summary_stats_pacbio.csv
                echo "mean:,"$inputreads","$V1","$V2","$V3","$V4","$V5","$V6","$V7"" >> all_summary_stats_pacbio.csv
                '''
              }
            }

            //
            // All reads section
            //

            // Filters allreads.csv based on set parameters and
            // recalculates relative frequencies after filter.
            process filterReads {
              container "quay.io/greninger-lab/tprk"

              // Retry on fail at most three times
              errorStrategy 'retry'
              maxRetries 3

              publishDir params.OUTDIR, mode: 'copy'

              input:
              file "allreads.csv" from allreads_ch
              file METADATA_FILE
              file FILTER_ALL_READS
              file RECALCULATE_FREQUENCY

              output:
              file "allreads_filtered.csv" into allreads_filt_ch
              file "allreads_filtered_heatmap.csv" into allreads_filt_heatmap_ch
              file "allreads.csv" into allreads_ch2
              file "allreads_filtered.csv" into filt_subset_ch

              script:
              """
              # Creates allreads_filtered.csv and recalculates the relative frequencies.
              python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a "allreads.csv"
              Rscript ${RECALCULATE_FREQUENCY} -f "allreads_filtered.csv" -m ${METADATA_FILE}

              # Creates allreads_filtered_heatmap.csv and recalculates the relative frequencies. This csv includes samples under the count/relative freq filters
              # if another sample shares the same read.
              python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a "allreads.csv" -is_heatmap
              Rscript ${RECALCULATE_FREQUENCY} -f "allreads_filtered_heatmap.csv" -m ${METADATA_FILE}

              """
            }