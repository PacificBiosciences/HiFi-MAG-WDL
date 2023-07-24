version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "bam_to_fastq/bam_to_fastq.wdl" as BamToFastq
import "assembly/assembly.wdl" as Assembly
import "completeness_aware_binning/completeness_aware_binning.wdl" as CompletenessAwareBinning
import "coverage/coverage.wdl" as Coverage
import "binning/binning.wdl" as Binning


workflow metagenomics {
	input {
		String sample_id
		File? bam
		File? fastq
		File db

		Int min_length = 500000
		Int min_completeness = 93

		# Backend configuration
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn
		String? container_registry

		Boolean preemptible
	}

	call BackendConfiguration.backend_configuration {
		input:
			backend = backend,
			zones = zones,
			aws_spot_queue_arn = aws_spot_queue_arn,
			aws_on_demand_queue_arn = aws_on_demand_queue_arn,
			container_registry = container_registry
	}

	RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

	Boolean run_bam_to_fastq = if (defined(bam)) then true else false

	if (run_bam_to_fastq) {
		call BamToFastq.bam_to_fastq {
			input:
				sample_id = sample_id,
				bam = select_first([bam]),
				default_runtime_attributes = default_runtime_attributes
		}
	}

	call Assembly.assembly {
		input:
			sample_id = sample_id,
			fastq = select_first([bam_to_fastq.converted_fastq, fastq]),
			default_runtime_attributes = default_runtime_attributes
	}

	call CompletenessAwareBinning.completeness_aware_binning {
		input:
			sample_id = sample_id,
			contigs_fasta = assembly.primary_contig_fasta,
			min_length = min_length,
			min_completeness = min_completeness,
			db = db,
			default_runtime_attributes = default_runtime_attributes
	}

	call Coverage.coverage {
		input:
			sample_id = sample_id,
			contigs_fasta = assembly.primary_contig_fasta,
			reads_fasta = assembly.reads_fasta,
			key = completeness_aware_binning.key,
			default_runtime_attributes = default_runtime_attributes
	}

	call Binning.binning {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = completeness_aware_binning.incomplete_contigs,
			filtered_depth = coverage.filtered_depth,
			sorted_bam = coverage.sorted_bam,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# BAM to FASTQ output
		File? converted_fastq = bam_to_fastq.converted_fastq

		# Assembly output
		File primary_contig_graph = assembly.primary_contig_graph
		File primary_contig_fasta = assembly.primary_contig_fasta
		File reads_fasta = assembly.reads_fasta

		# Completeness-aware binning output
		File key = completeness_aware_binning.key
		File incomplete_contigs = completeness_aware_binning.incomplete_contigs
		File? report = completeness_aware_binning.report
		File? passed_bins = completeness_aware_binning.passed_bins
		File? scatterplot = completeness_aware_binning.scatterplot
		File? histogram = completeness_aware_binning.histogram

		# Coverage output
		File sorted_bam = coverage.sorted_bam
		File sorted_bam_index = coverage.sorted_bam_index
		File filtered_depth = coverage.filtered_depth

		# Binning output
		Array[File] metabat2_reconstructed_bins_fastas = binning.metabat2_reconstructed_bins_fastas
		File metabat2_bin_sets_tsv = binning.metabat2_bin_sets_tsv
		File semibin2_bins_tsv = binning.semibin2_bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = binning.semibin2_reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = binning.semibin2_bin_sets_tsv
		Array[File] dastool_bins = binning.dastool_bins

	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		# BAM to FASTQ
		bam: {help: "Sample BAM to convert to fastq format"}
		# Assembly
		fastq: {help: "Sample in fastq format"}
		# Completeness-aware binning
		db: {help: "CheckM2 DIAMOND reference database Uniref100/KO"}
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		# Backend configuration
		backend: {help: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS']"}
		zones: {help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"}
		aws_spot_queue_arn: {help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"}
		aws_on_demand_queue_arn: {help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"}
		container_registry: {help: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used"}
		preemptible: {help: "Where possible, run tasks preemptibly"}
	}
}
