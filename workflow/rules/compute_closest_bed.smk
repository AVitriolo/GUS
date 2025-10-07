rule compute_closest_bed:
    output:
        ""
    input:
        input_1 = "",
        input_2 = ""
    shell:
        """
	     closestBed -a {input.input_1} -b {input.input_2} -d > {output}
        """



