__author__ = "Alice Iliasova"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

configfile: "configCNV.yaml"

# include functions.py modules
include: "SCRIPTS/functions_CNV.py"

# extend the ALL rule using python extend list function
ALL = []

#TARGETS

#PREP_INTERVALS=["path/to/preprocessed.intervals_list.interval_list"]
#COUNTS=expand("NEXTERA_AROS/{sample}.counts.hdf5", sample=SAMPLES)
#PON=["NEXTERA_AROS/CNV_PON_NEXTERA_AROS.pon.hdf5"]
#SCR=["TSVs/LYM2.standardizedCR.tsv"],
#DCR=["TSVs/LYM2.denoisedCR.tsv"]
#PLOTS=["PLOTS"]
#NAC=["TSVs/LYM2.NormalAllelicCounts.tsv"]
#TAC=["TSVs/LYM2.TumourAllelicCounts.tsv"]
#SEG=["SEG"]
#SEG_TO=["SEG_TO"]
#CALLED_CR=["SEG_TO/LYM2.called.seg"]
BAF_PLOTS=["BAF_PLOTS"]
#BAF_PLOTS_TO=["BAF_PLOTS_TO"]


#ALL.extend(PREP_INTERVALS)
#ALL.extend(COUNTS)
#ALL.extend(PON)
#ALL.extend(SCR)
#ALL.extend(DCR)
#ALL.extend(PLOTS)
#ALL.extend(NAC)
#ALL.extend(TAC)
#ALL.extend(SEG)
#ALL.extend(SEG_TO)
#ALL.extend(CALLED_CR)
ALL.extend(BAF_PLOTS)

rule ALL:
        input:
                ALL,

rule PreprocessIntervals:
        input:
                INTERVALS=config["INTERVALS"]
        output: PreprocessedINTERVALS="path/to/preprocessed.intervals_list.interval_list"
        params:
                REF=config["REF"]


        threads: 4

        log:"LOGS/CNV/PreprocessIntervals.log"

        benchmark:"LOGS/CNV/PreprocessIntervals.txt"

        message:"Running GATK PreprocessIntervals for {input} using {threads} threads and saving into {output}"

        shell:"""
        gatk PreprocessIntervals \
                -L {input.INTERVALS} \
                -R {params.REF} \
                --bin-length 0 \
                --interval-merging-rule OVERLAPPING_ONLY \
                -O {output.PreprocessedINTERVALS}
        """

rule CollectCounts:
        input:
                unpack(ReadCounts_input)
        output:
                Counts="NEXTERA_AROS/{sample}.counts.hdf5"
        params:
                PreprocessedINTERVALS="path/to/preprocessed.intervals_list.interval_list"


        threads: 4

        log:"LOGS/CNV/CollectCounts.{sample}.log"

        benchmark:"LOGS/CNV/CollectCounts.{sample}.txt"

        message:"Running GATK CollectFragmentCounts for {input} and saving in {output}"

        shell:"""
        gatk CollectReadCounts \
                -I {input.File} \
                -L {params.PreprocessedINTERVALS} \
                --interval-merging-rule OVERLAPPING_ONLY \
                -O {output.Counts}
        """


rule CreateReadCountPanelOfNormals:
	input:
		PON_Input="PON_Input.tsv",
        
	output:
		CNVponC="NEXTERA_AROS/CNV_PON_NEXTERA_AROS.pon.hdf5"
	

	threads: 4

	log:"LOGS/CNV/CNVponC.log"

	benchmark:"LOGS/CNV/CNVponC.txt"

	message: "Running GATK CollectFragmentCounts for {input} using {threads} threads and saving as {output}"

	shell:"""		
	gatk --java-options "-Xmx6500m" CreateReadCountPanelOfNormals \
		--arguments_file {input.PON_Input} \
		--minimum-interval-median-percentile 5.0 \
		-O {output.CNVponC}
	"""

rule DenoiseReadCounts:
	input:
		I="NEXTERA_AROS/LYM2.counts.hdf5",
		#PREFIX=config["PREFIX"]
	output:
		StandardisedCopyRatios="TSVs/LYM2.standardizedCR.tsv",
		DenoisedCopyRatios="TSVs/LYM2.denoisedCR.tsv"
	params:
		CNVponC="NEXTERA_CfL/CNV_PON_NEXTERA_CfL.pon.hdf5",
		
		
	threads: 4

	log:"LOGS/CNV/DenoisedCopyRatios.log"

	benchmark:"LOGS/CNV/DenoisedCopyRatios.txt"

	message: "Running GATK DenoiseReadCounts for {input} using {threads} threads and saving as {output}"

	shell:"""		
	gatk --java-options "-Xmx12g" DenoiseReadCounts \
		-I {input.I} \
		--count-panel-of-normals {params.CNVponC} \
		--standardized-copy-ratios {output.StandardisedCopyRatios} \
		--denoised-copy-ratios {output.DenoisedCopyRatios}
    """
rule PlotDenoisedCopyRatios:
	input:
		StandardisedCopyRatios="TSVs/LYM2.standardizedCR.tsv",
		DenoisedCopyRatios="TSVs/LYM2.denoisedCR.tsv"
	output:
		PLOTS=directory("PLOTS")
	params:
		REF_DICT=config["REF_DICT"],
		
		
	threads: 4

	log:"LOGS/CNV/PlotDenoisedCopyRatios.log"

	benchmark:"LOGS/CNV/PlotDenoisedCopyRatios.txt"

	message: "Running GATK PlotDenoisedCopyRatios for {input} using {threads} threads and saving as {output}"

	shell:"""			
	gatk PlotDenoisedCopyRatios \
		--standardized-copy-ratios {input.StandardisedCopyRatios} \
		--denoised-copy-ratios {input.DenoisedCopyRatios} \
		--sequence-dictionary {params.REF_DICT} \
		--minimum-contig-length 46709983 \
		--output {output.PLOTS} \
		--output-prefix LYM2
	"""

rule CollectAllelicCounts_N:
	input:
        	NORMAL="/path/to/PoN/bam//Batch2_70_LYM6.recal.bam"
	output: 
		NormalAllelicCounts="TSVs/LYM6.NormalAllelicCounts.tsv"
	params:
		INTERVALS=config["COMMON_SNPS"],
		REF=config["REF"],
		
	threads: 8

	log:"LOGS/CNV/AllelicCounts_N.log"

	benchmark:"LOGS/CNV/AllelicCounts_N.txt"

	message: "Running GATK CollectAllelicCounts  for {input} using {threads} threads and saving as {output}"

	shell:"""
	gatk --java-options "-Xmx36g -Djava.io.tmpdir=/scratch" CollectAllelicCounts \
		-L {params.INTERVALS} \
		-I {input.NORMAL} \
		-R {params.REF} \
		-O {output.NormalAllelicCounts} 2>>{log}
    """

rule CollectAllelicCounts_T:
	input:
		TUMOUR="/nobackup/proj/sbchtrbl/BSU/LYM/bam/62_LYM2_Vikki_Rand.recal.bam"
	output:
		TumourAllelicCounts="TSVs/LYM2.TumourAllelicCounts.tsv"
	params:
		INTERVALS=config["COMMON_SNPS"],
		REF=config["REF"],

	threads: 8

	log:"LOGS/CNV/AllelicCountsT.log"

	benchmark:"LOGS/CNV/AllelicCounts.txt"

	message: "Running GATK  CollectAllelicCounts for {input} using {threads} threads and saving as {output}"
	
	shell:"""
	gatk --java-options "-Xmx36g -Djava.io.tmpdir=/scratch" CollectAllelicCounts \
		-L {params.INTERVALS} \
		-I {input.TUMOUR} \
		-R {params.REF} \
		-O {output.TumourAllelicCounts} 2>>{log}
	"""

rule ModelSegments:
	input:
		NormalAllelicCounts="TSVs/LYM2.NormalAllelicCounts.tsv",
		TumourAllelicCounts="TSVs/LYM2.TumourAllelicCounts.tsv",
		DenoisedCopyRatios="TSVs/LYM2.denoisedCR.tsv"
	output: 
		Segments=directory("SEG")
	#params:
		#Out_Prefix="LYM4"
		
		
	threads: 4

	log:"LOGS/CNV/ModelSegments.log"

	benchmark:"LOGS/CNV/ModelSegments.txt"

	message: "Running GATK ModelSegments for {input} using {threads} threads and saving as {output}"

	shell:"""  
	gatk --java-options "-Xmx4g" ModelSegments \
		--denoised-copy-ratios {input.DenoisedCopyRatios} \
		--allelic-counts {input.TumourAllelicCounts} \
		--normal-allelic-counts {input.NormalAllelicCounts} \
		--output {output.Segments}\
		--output-prefix LYM2
    """

rule ModelSegments_TumourOnly:
		input:
			TumourAllelicCounts="TSVs/LYM2.TumourAllelicCounts.tsv",
			DenoisedCopyRatios="TSVs/LYM2.denoisedCR.tsv"
		output:
			Segments=directory("SEG_TO")
        #params:
                #Out_Prefix="LYM4"


		threads: 4

		log:"LOGS/CNV/ModelSegments.log"

		benchmark:"LOGS/CNV/ModelSegments.txt"
		
		message: "Running GATK ModelSegments for {input} using {threads} threads and saving as {output}"

		shell:"""
		gatk --java-options "-Xmx4g" ModelSegments \
			--denoised-copy-ratios {input.DenoisedCopyRatios} \
			--allelic-counts {input.TumourAllelicCounts} \
			--output-prefix LYM2 \
			--output {output.Segments}
		"""

rule CallCopyRatioSegments:
	input:
		CR_SEG="SEG_TO/LYM2.cr.seg"
       
	output: 
		Called_CR="SEG_TO/LYM2.called.seg"
		
	threads: 4

	log:"LOGS/CNV/CallCopyRatioSegments.log"

	benchmark:"LOGS/CNV/CallCopyRatioSegments.txt"

	message: "Running GATK CallCopyRatioSegments for {input} using {threads} threads and saving as {output}"

	shell:""" 
	gatk CallCopyRatioSegments \
		--input {input.CR_SEG} \
		--output {output.Called_CR}
    """

rule PlotModeledSegments:
	input:
		DenoisedCopyRatios="TSVs/LYM2.denoisedCR.tsv",
		AllelicCounts="SEG_TO/LYM2.hets.tsv",
		Segments="SEG_TO/LYM2.modelFinal.seg",
		Dict=config["REF_DICT"]
	output: 
		BAF_PLOTS=directory("BAF_PLOTS")
	
	threads: 4

	log:"LOGS/CNV/PlotModeledSegments.log"

	benchmark:"LOGS/CNV/PlotModeledSegments.txt"

	message: "Running GATK PlotModeledSegments for {input} using {threads} threads and saving as {output}"

	shell:"""   
	gatk PlotModeledSegments \
		--denoised-copy-ratios {input.DenoisedCopyRatios} \
		--allelic-counts {input.AllelicCounts} \
		--segments {input.Segments} \
		--sequence-dictionary {input.Dict} \
		--minimum-contig-length 46709983 \
		--output {output.BAF_PLOTS} \
		--output-prefix LYM2
    """
rule PlotModeledSegments_TO:
		input:
			AllelicCounts="SEG_TO/LYM2.hets.tsv",
			Segments="SEG_TO/LYM2.modelFinal.seg",
			Dict=config["REF_DICT"]
		output:
			BAF_PLOTS_TO=directory("BAF_PLOTS_TO")

		threads: 4
		log:"LOGS/CNV/PlotModeledSegments.log"

		benchmark:"LOGS/CNV/PlotModeledSegments.txt"

		message: "Running GATK PlotModeledSegments for {input} using {threads} threads and saving as {output}"

		shell:"""
		gatk PlotModeledSegments \
			--allelic-counts {input.AllelicCounts} \
			--segments {input.Segments} \
			--sequence-dictionary {input.Dict} \
			--minimum-contig-length 46709983 \
			--output {output.BAF_PLOTS_TO} \
			--output-prefix LYM2
	"""
