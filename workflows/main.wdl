version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "completeness/completeness.wdl" as Completeness

workflow metagenomics {
	input {
		String sample
		File contig

		Int min_length = 500000
		#Int min_completeness = 93

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

	call Completeness.completeness {
		input:
			sample = sample,
			contig = contig,
			min_length = min_length,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# Completeness output
		File key = completeness.key
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contig: {help: "Contigs"} #TODO
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		key: {help: ""} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}
