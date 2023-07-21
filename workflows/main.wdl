version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "assembly/assembly.wdl" as Assembly
import "completeness_aware_binning/completeness_aware_binning.wdl" as CompletenessAwareBinning
import "depth/depth.wdl" as Depth


workflow metagenomics {
	input {
		String sample
		#File contig_gfa
		#File db

		# TODO - remove
		File contig
		File reads_fasta
		File passed

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

	#call Assembly.assembly {
	#	input:
	#		sample = sample,
	#		fastq = fastq,
	#		default_runtime_attributes = default_runtime_attributes
	#}

	#call CompletenessAwareBinning.completeness_aware_binning {
	#	input:
	#		sample = sample,
	#		contig = assembly.primary_contig_fasta,
	#		min_length = min_length,
	#		min_completeness = min_completeness,
	#		db = db,
	#		default_runtime_attributes = default_runtime_attributes
	#}

	call Depth.depth {
		input:
			sample = sample,
			contig = contig,
			reads_fasta = reads_fasta,
			passed = passed,
			#contig = assembly.primary_contig_fasta,
			#reads_fasta = assembly.reads_fasta,
			#passed = completeness_aware_binning.key,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# Assembly output
		#File primary_contig_graph = assembly.primary_contig_graph
		#File primary_contig_fasta = assembly.primary_contig_fasta
		#File reads_fasta = assembly.reads_fasta

		# Completeness output
		#File key = completeness_aware_binning.key

		# Depth output
		File bam = depth.bam
		File bam_index = depth.bam_index
		File filtered_depth = depth.filtered_depth
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contig_gfa: {help: "Contigs"} #TODO
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		key: {help: ""} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}
