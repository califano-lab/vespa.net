# number of ARACNe bootstrap iterations
seed = list(range(1,201))

rule all:
    input:
        "results/ddpi_substrate_regulon/regulon.rds", "results/hsm_substrate_regulon/regulon.rds", "results/dpi_activity_regulon/regulon.rds"

# prepare substrate regulon data
rule prepare_substrate_regulon:
    input:
        phospho = "phospho.rds",
        proteo = "proteo.rds",
        ref = "phospho.rds"
    output:
        kinases = "results/prepare_substrate_regulon/kinases.txt",
        kinases_phosphatases = "results/prepare_substrate_regulon/kinases_phosphatases.txt",
        targets = "results/prepare_substrate_regulon/targets.txt",
        phosphointeractions = "results/prepare_substrate_regulon/phosphointeractions.txt",
        peptides = "results/prepare_substrate_regulon/peptides.txt",
        matrix = "results/prepare_substrate_regulon/matrix.txt",
    script:
        "scripts/prepare_substrate_regulon.R"

# generate dDPI regulon
rule ddpi_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/ddpi_substrate_regulon/fwer_computed.txt"
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s 1 -t && touch {output}"

rule ddpi_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.ddpi_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/ddpi_substrate_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule ddpi_substrate_regulon_consolidate:
    input:
        iteration = expand("results/ddpi_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/ddpi_substrate_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -o $(dirname {output}) -c"

rule ddpi_substrate_regulon_generate:
    input:
        network = rules.ddpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/ddpi_substrate_regulon/regulon.rds"
    script:
        "scripts/generate_regulon.R"

# generate HSM/D regulon
rule hsm_substrate_regulon_mit:
    input:
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/hsm_substrate_regulon/fwer_computed.txt"
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -o $(dirname {output}) -s 1 -t && touch {output}"

rule hsm_substrate_regulon_bs:
    input:
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.hsm_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/hsm_substrate_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_substrate_regulon_consolidate:
    input:
        iteration = expand("results/hsm_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/hsm_substrate_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -o $(dirname {output}) -c"

rule hsm_substrate_regulon_generate:
    input:
        network = rules.hsm_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/hsm_substrate_regulon/regulon.rds"
    script:
        "scripts/generate_regulon.R"

# prepare activity regulon data
rule prepare_activity_regulon:
    input:
        phospho = "phospho.rds",
        proteo = "proteo.rds",
        ref = "phospho.rds",
        ddpi_substrate_regulon = rules.ddpi_substrate_regulon_generate.output.regulon,
        hsm_substrate_regulon = rules.hsm_substrate_regulon_generate.output.regulon

    output:
        kinases = "results/prepare_activity_regulon/kinases.txt",
        kinases_phosphatases = "results/prepare_activity_regulon/kinases_phosphatases.txt",
        targets = "results/prepare_activity_regulon/targets.txt",
        phosphointeractions = "results/prepare_activity_regulon/phosphointeractions.txt",
        peptides = "results/prepare_activity_regulon/peptides.txt",
        matrix = "results/prepare_activity_regulon/matrix.txt"
    script:
        "scripts/prepare_activity_regulon.R"

# generate DPI regulon
rule dpi_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix
    output:
        mit = "results/dpi_activity_regulon/fwer_computed.txt"
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s 1 -t && touch {output}"

rule dpi_activity_regulon_bs:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix,
        mit = rules.dpi_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/dpi_activity_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule dpi_activity_consolidate:
    input:
        iteration = expand("results/dpi_activity_regulon/{seed}", seed=seed)
    output:
        network = "results/dpi_activity_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx12G -jar java/aracne.jar -o $(dirname {output}) -c"

rule dpi_activity_regulon_generate:
    input:
        network = rules.dpi_activity_consolidate.output.network,
        matrix = rules.prepare_activity_regulon.output.matrix,
        peptides = rules.prepare_activity_regulon.output.peptides
    output:
        regulon = "results/dpi_activity_regulon/regulon.rds"
    script:
        "scripts/generate_regulon.R"
