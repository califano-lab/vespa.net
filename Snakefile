# number of ARACNe bootstrap iterations
seed = list(range(1,201))

# obtain dataset ids from phospho data
dsids, = glob_wildcards("{dsid}_phospho.rds")

rule all:
    input:
        "results/ddpi_substrate_site_regulon.rds", "results/hsm_substrate_site_regulon.rds", "results/pc_substrate_site_regulon.rds"

# prepare substrate regulon data
rule prepare_substrate_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        ref = "reference.rds"
    output:
        kinases = "results/{dsid}/prepare_substrate_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_substrate_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_substrate_regulon/targets.txt",
        phosphointeractions = "results/{dsid}/prepare_substrate_regulon/phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_substrate_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_substrate_regulon/matrix.txt",
    script:
        "scripts/prepare_substrate_regulon.R"

# prepare substrate regulon data
rule prepare_pcsubstrate_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        ref = "reference.rds"
    output:
        kinases = "results/{dsid}/prepare_pcsubstrate_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_pcsubstrate_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_pcsubstrate_regulon/targets.txt",
        phosphointeractions = "results/{dsid}/prepare_pcsubstrate_regulon/phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_pcsubstrate_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_pcsubstrate_regulon/matrix.txt",
    script:
        "scripts/prepare_pcsubstrate_regulon.R"

# generate dDPI regulon
rule ddpi_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/ddpi_substrate_regulon/fwer_computed.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule ddpi_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.ddpi_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/ddpi_substrate_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule ddpi_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/ddpi_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/ddpi_substrate_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -o $(dirname {output}) -c -j {threads}"

rule ddpi_substrate_regulon_generate:
    input:
        network = rules.ddpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/ddpi_substrate_regulon.rds"
    script:
        "scripts/generate_regulon.R"

# generate HSM/D regulon
rule hsm_substrate_regulon_mit:
    input:
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/hsm_substrate_regulon/fwer_computed.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule hsm_substrate_regulon_bs:
    input:
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.hsm_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/hsm_substrate_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/hsm_substrate_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -o $(dirname {output}) -c -j {threads}"

rule hsm_substrate_regulon_generate:
    input:
        network = rules.hsm_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_substrate_regulon.rds"
    script:
        "scripts/generate_regulon.R"


# generate PC regulon
rule pc_substrate_regulon_mit:
    input:
        phosphointeractions = rules.prepare_pcsubstrate_regulon.output.phosphointeractions,
        targets = rules.prepare_pcsubstrate_regulon.output.targets,
        matrix = rules.prepare_pcsubstrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/pc_substrate_regulon/fwer_computed.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule pc_substrate_regulon_bs:
    input:
        phosphointeractions = rules.prepare_pcsubstrate_regulon.output.phosphointeractions,
        targets = rules.prepare_pcsubstrate_regulon.output.targets,
        matrix = rules.prepare_pcsubstrate_regulon.output.matrix,
        mit = rules.pc_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/pc_substrate_regulon/{seed}")
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -e {input.matrix} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule pc_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/pc_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/pc_substrate_regulon/network.txt"
    threads: 1
    shell:
        "java -Xmx4G -jar java/aracne.jar -o $(dirname {output}) -c -j {threads}"

rule pc_substrate_regulon_generate:
    input:
        network = rules.pc_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_pcsubstrate_regulon.output.matrix,
        peptides = rules.prepare_pcsubstrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/pc_substrate_regulon.rds"
    script:
        "scripts/generate_regulon.R"

# generate ddpi substrate regulons
rule meta_ddpi_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/ddpi_substrate_regulon.rds", dsid=dsids)
    output:
        meta_site_regulons = "results/ddpi_substrate_site_regulon.rds",
        meta_protein_regulons = "results/ddpi_substrate_protein_regulon.rds",
    threads: 4
    script:
        "scripts/generate_meta_regulon.R"

rule meta_hsm_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/hsm_substrate_regulon.rds", dsid=dsids)
    output:
        meta_site_regulons = "results/hsm_substrate_site_regulon.rds",
        meta_protein_regulons = "results/hsm_substrate_protein_regulon.rds",
    threads: 4
    script:
        "scripts/generate_meta_regulon.R"

rule meta_pc_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/pc_substrate_regulon.rds", dsid=dsids)
    output:
        meta_site_regulons = "results/pc_substrate_site_regulon.rds",
        meta_protein_regulons = "results/pc_substrate_protein_regulon.rds",
    threads: 4
    script:
        "scripts/generate_meta_regulon.R"
