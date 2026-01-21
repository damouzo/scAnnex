/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility Functions for scAnnex Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Functions to ensure that resource requirements don't go beyond maximum limits
----------------------------------------------------------------------------------------
*/

class Utils {
    
    /**
     * Check that a value doesn't exceed the maximum allowed resource
     * @param obj The resource value to check (memory, time, or cpus)
     * @param type The type of resource ('memory', 'time', or 'cpus')
     * @param params The params object containing max_memory, max_time, max_cpus
     * @return The value capped at the maximum if exceeded, or the original value
     */
    static def checkMax(obj, type, params) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min(obj, params.max_cpus as int)
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }
}
