version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "assembly/assembly.wdl" as Assembly

workflow metagenomics {
	input {
		String sample
		File fastq

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

	call Assembly.assembly {
		input:
			sample = sample,
			fastq = fastq,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# Assembly output
	}

	parameter_meta {
		sample: {help: "Sample name"}
		fastq: {help: "Sample"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}
