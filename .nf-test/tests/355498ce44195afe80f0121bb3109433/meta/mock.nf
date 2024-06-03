import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { SEQKIT_STATS } from '/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/local/seqkit/stats/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .excludeNulls()  // Do not include fields with value null..
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                input[0] = [
                    [id: 'test'],
                    file("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/assets/test-datasets/fastas/MEGAHIT-BUFF_01_B.contigs.fa.gz", checkIfExists: true)
                ]
                
    //----

    //run process
    SEQKIT_STATS(*input)

    if (SEQKIT_STATS.output){

        // consumes all named output channels and stores items in a json file
        for (def name in SEQKIT_STATS.out.getNames()) {
            serializeChannel(name, SEQKIT_STATS.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = SEQKIT_STATS.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonOutput.toJson(result)
    
}
