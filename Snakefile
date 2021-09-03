# vespa.net version 1.0.2

# number of ARACNe bootstrap iterations
seed = list(range(1,201))

# obtain dataset ids from phospho data
dsids, = glob_wildcards("{dsid}_phospho.rds")

rule all:
    input:
        "results/ddpimeta_substrate_site_regulon.rds", "results/ddpimeta_substrate_protein_regulon.rds",
        "results/dpimeta_substrate_site_regulon.rds", "results/ddimeta_substrate_protein_regulon.rds",
        "results/nodpimeta_substrate_site_regulon.rds", "results/nodpimeta_substrate_protein_regulon.rds",
        "results/hsmmeta_substrate_site_regulon.rds", "results/hsmmeta_substrate_protein_regulon.rds",
        "results/lpmeta_substrate_site_regulon.rds", "results/lpmeta_substrate_protein_regulon.rds"

# prepare substrate regulon data
rule prepare_substrate_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        #ref = "reference.rds" # Uncomment if restrict_peptides = True
    output:
        kinases = "results/{dsid}/prepare_substrate_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_substrate_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_substrate_regulon/targets.txt",
        hsm_phosphointeractions = "results/{dsid}/prepare_substrate_regulon/hsm_phosphointeractions.txt",
        pc_phosphointeractions = "results/{dsid}/prepare_substrate_regulon/pc_phosphointeractions.txt",
        lp_phosphointeractions = "results/{dsid}/prepare_substrate_regulon/lp_phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_substrate_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_substrate_regulon/matrix.txt",
    params:
        hsm_threshold = 0,
        restrict_peptides = False
    singularity:
        "vespa.simg"
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
        "java -Xmx28G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

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
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule ddpi_substrate_regulon_generate:
    input:
        network = rules.ddpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/ddpi_substrate_regulon.rds"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

# generate DPI regulon
rule dpi_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/dpi_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx28G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule dpi_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.dpi_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/dpi_substrate_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule dpi_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/dpi_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/dpi_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule dpi_substrate_regulon_generate:
    input:
        network = rules.dpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/dpi_substrate_regulon.rds"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

# generate noDPI regulon
rule nodpi_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/nodpi_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx28G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} --noDPI -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule nodpi_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.nodpi_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/nodpi_substrate_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} --noDPI -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule nodpi_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/nodpi_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/nodpi_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule nodpi_substrate_regulon_generate:
    input:
        network = rules.nodpi_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/nodpi_substrate_regulon.rds"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

# generate HSM/D regulon
rule hsm_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.hsm_phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/hsm_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule hsm_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.hsm_phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.hsm_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/hsm_substrate_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/hsm_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule hsm_substrate_regulon_generate:
    input:
        network = rules.hsm_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_substrate_regulon.rds"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

# generate LP regulon
rule lp_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.lp_phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/lp_substrate_regulon/fwer_computed.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule lp_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_substrate_regulon.output.lp_phosphointeractions,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.lp_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/lp_substrate_regulon/{seed}")
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule lp_substrate_regulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/lp_substrate_regulon/{seed}", seed=seed)
    output:
        network = "results/{dsid}/lp_substrate_regulon/network.txt"
    threads: 1
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule lp_substrate_regulon_generate:
    input:
        network = rules.lp_substrate_regulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/lp_substrate_regulon.rds"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

# generate meta substrate regulons
rule meta_substrate_regulon_generate:
    input:
        ref = "reference.rds",
        substrate_regulons = [],
        regulons = [expand("results/{dsid}/ddpi_substrate_regulon.rds", dsid=dsids), expand("results/{dsid}/dpi_substrate_regulon.rds", dsid=dsids), expand("results/{dsid}/nodpi_substrate_regulon.rds", dsid=dsids), expand("results/{dsid}/hsm_substrate_regulon.rds", dsid=dsids), expand("results/{dsid}/lp_substrate_regulon.rds", dsid=dsids)]
    output:
        meta_redundantsite_regulons = "results/meta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/meta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/meta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate dDPI-meta substrate regulons
rule ddpimeta_substrate_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/ddpi_substrate_regulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/ddpimeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/ddpimeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/ddpimeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate DPI-meta substrate regulons
rule ddpimeta_substrate_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/dpi_substrate_regulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/dpimeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/dpimeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/dpimeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate noDPI-meta substrate regulons
rule ddpimeta_substrate_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/nodpi_substrate_regulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/nodpimeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/nodpimeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/nodpimeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate HSM-meta substrate regulons
rule hsmmeta_substrate_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/hsm_substrate_regulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/hsmmeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/hsmmeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/hsmmeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate LP-meta substrate regulons
rule lpmeta_substrate_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = [],
        regulons = expand("results/{dsid}/lp_substrate_regulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/lpmeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/lpmeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/lpmeta_substrate_protein_regulon.rds",
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20,
        orthogonal_cutoff = 0.5,
        transform = "zscore"
    threads: 4
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"
