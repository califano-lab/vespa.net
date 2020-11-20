# number of ARACNe bootstrap iterations
seed = list(range(1,201))

# obtain dataset ids from phospho data
dsids, = glob_wildcards("{dsid}_phospho.rds")

rule all:
    input:
        "results/meta_substrate_site_regulon.rds", "results/meta_activity_site_regulon.rds", "results/meta_substrate_protein_regulon.rds", "results/meta_activity_protein_regulon.rds",
        "results/ddpimeta_substrate_site_regulon.rds", "results/dpimeta_activity_site_regulon.rds", "results/ddpimeta_substrate_protein_regulon.rds", "results/dpimeta_activity_protein_regulon.rds",
        "results/hsmmeta_substrate_site_regulon.rds", "results/hsmmeta_activity_site_regulon.rds", "results/hsmmeta_substrate_protein_regulon.rds", "results/hsmmeta_activity_protein_regulon.rds"

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
    params:
        hsm_threshold = 0.05,
        restrict_peptides = False
    singularity:
        "phosphoviper.simg"
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
        mit = "results/{dsid}/ddpi_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

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
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule ddpi_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/ddpi_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/ddpi_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads}"

rule ddpi_substrate_regulon_generate:
    input:
        network = rules.ddpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/ddpi_substrate_regulon.rds"
    params:
        likelihood_threshold = 0.5,
        priors = False
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_regulon.R"

# generate HSM/D regulon
rule hsm_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/hsm_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -i {input.phosphointeractions} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule hsm_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.hsm_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/hsm_substrate_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar --noDPI -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -i {input.phosphointeractions} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/hsm_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads}"

rule hsm_substrate_regulon_generate:
    input:
        network = rules.hsm_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_substrate_regulon.rds"
    params:
        likelihood_threshold = 0.5,
        priors = False
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_regulon.R"

# generate meta substrate regulons
rule meta_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = [expand("results/{dsid}/ddpi_substrate_regulon.rds", dsid=dsids), expand("results/{dsid}/hsm_substrate_regulon.rds", dsid=dsids)]
    output:
        meta_site_regulons = "results/meta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/meta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate dDPI-meta substrate regulons
rule ddpimeta_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/ddpi_substrate_regulon.rds", dsid=dsids)
    output:
        meta_site_regulons = "results/ddpimeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/ddpimeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate HSM-meta substrate regulons
rule hsmmeta_substrate_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/hsm_substrate_regulon.rds", dsid=dsids)
    output:
        meta_site_regulons = "results/hsmmeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/hsmmeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"

# prepare activity regulon data
rule prepare_activity_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        meta_substrate_regulons = rules.meta_substrate_regulon_generate.output.meta_protein_regulons,
        fasta = "library.fasta"
    output:
        kinases = "results/{dsid}/prepare_activity_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_activity_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_activity_regulon/targets.txt",
        phosphointeractions = "results/{dsid}/prepare_activity_regulon/phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_activity_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_activity_regulon/matrix.txt"
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        hsm_threshold = 0.05,
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/prepare_activity_regulon.R"

# generate DPI regulon
rule dpi_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix
    output:
        mit = "results/{dsid}/dpi_activity_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule dpi_activity_regulon_bs:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix,
        mit = rules.dpi_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/dpi_activity_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule dpi_activity_consolidate:
    input:
        iteration = expand("results/{{dsid}}/dpi_activity_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/dpi_activity_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads}"

rule dpi_activity_regulon_generate:
    input:
        network = rules.dpi_activity_consolidate.output.network,
        matrix = rules.prepare_activity_regulon.output.matrix,
        peptides = rules.prepare_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/dpi_activity_regulon.rds"
    params:
        likelihood_threshold = 0.5,
        priors = False
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_regulon.R"

# generate HSM/P regulon
rule hsm_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_activity_regulon.output.phosphointeractions,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix
    output:
        mit = "results/{dsid}/hsm_activity_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule hsm_activity_regulon_bs:
    input:
        kinases_phosphatases = rules.prepare_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_activity_regulon.output.phosphointeractions,
        targets = rules.prepare_activity_regulon.output.targets,
        matrix = rules.prepare_activity_regulon.output.matrix,
        mit = rules.hsm_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/hsm_activity_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar --noDPI -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_activity_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_activity_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/hsm_activity_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads}"

rule hsm_activity_regulon_generate:
    input:
        network = rules.hsm_activity_consolidate.output.network,
        matrix = rules.prepare_activity_regulon.output.matrix,
        peptides = rules.prepare_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_activity_regulon.rds"
    params:
        likelihood_threshold = 0.5,
        priors = False
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_regulon.R"

# generate meta activity regulons
rule meta_activity_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = rules.meta_substrate_regulon_generate.output.meta_protein_regulons,
        regulons = [expand("results/{dsid}/dpi_activity_regulon.rds", dsid=dsids), expand("results/{dsid}/hsm_activity_regulon.rds", dsid=dsids)],
        fasta = "library.fasta"
    output:
        meta_site_regulons = "results/meta_activity_site_regulon.rds",
        meta_protein_regulons = "results/meta_activity_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate DPI-meta activity regulons
rule dpimeta_activity_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = rules.meta_substrate_regulon_generate.output.meta_protein_regulons,
        regulons = expand("results/{dsid}/dpi_activity_regulon.rds", dsid=dsids),
        fasta = "library.fasta"
    output:
        meta_site_regulons = "results/dpimeta_activity_site_regulon.rds",
        meta_protein_regulons = "results/dpimeta_activity_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate HSM-meta activity regulons
rule hsmmeta_activity_regulon_generate:
    input:
        ref = rules.prepare_substrate_regulon.input.ref,
        substrate_regulons = rules.meta_substrate_regulon_generate.output.meta_protein_regulons,
        regulons = expand("results/{dsid}/hsm_activity_regulon.rds", dsid=dsids),
        fasta = "library.fasta"
    output:
        meta_site_regulons = "results/hsmmeta_activity_site_regulon.rds",
        meta_protein_regulons = "results/hsmmeta_activity_protein_regulon.rds",
    params:
        minimum_targets = 10,
        maximum_targets = 500,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 10,
        ct_penalty = 20
    threads: 4
    singularity:
        "phosphoviper.simg"
    script:
        "scripts/generate_meta_regulon.R"
