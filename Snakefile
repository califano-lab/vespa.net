# vespa.net version 1.0.2

# number of ARACNe bootstrap iterations
seed = list(range(1,201))

# obtain dataset ids from phospho data
dsids, = glob_wildcards("{dsid}_phospho.rds")

rule all:
    input:
        "results/stdpimeta_substrate_site_regulon.rds", "results/stdpimeta_substrate_redundantsite_regulon.rds", "results/stdpimeta_substrate_protein_regulon.rds", "results/dpimeta_activity_site_regulon.rds", "results/dpimeta_activity_redundantsite_regulon.rds", "results/dpimeta_activity_protein_regulon.rds",
        "results/hsmmeta_substrate_site_regulon.rds", "results/hsmmeta_substrate_redundantsite_regulon.rds", "results/hsmmeta_substrate_protein_regulon.rds", "results/hsmmeta_activity_site_regulon.rds", "results/hsmmeta_activity_redundantsite_regulon.rds", "results/hsmmeta_activity_protein_regulon.rds",
        "results/lpmeta_substrate_site_regulon.rds", "results/lpmeta_substrate_redundantsite_regulon.rds", "results/lpmeta_substrate_protein_regulon.rds", "results/lpmeta_activity_site_regulon.rds", "results/lpmeta_activity_redundantsite_regulon.rds", "results/lpmeta_activity_protein_regulon.rds",
        "results/meta_substrate_site_regulon.rds", "results/meta_substrate_redundantsite_regulon.rds", "results/meta_substrate_protein_regulon.rds", "results/meta_activity_site_regulon.rds", "results/meta_activity_redundantsite_regulon.rds", "results/meta_activity_protein_regulon.rds"

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
        lp_phosphointeractions = "results/{dsid}/prepare_substrate_regulon/lp_phosphointeractions.txt",
        pc_phosphointeractions = "results/{dsid}/prepare_substrate_regulon/pc_phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_substrate_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_substrate_regulon/matrix.txt",
    params:
        hsm_threshold = 0,
        restrict_peptides = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/prepare_substrate_regulon.R"

# generate stDPI substrate regulon
rule stdpi_substrate_regulon_mit:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix
    output:
        mit = "results/{dsid}/stdpi_substrate_regulon/fwer_computed.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx28G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule stdpi_substrate_regulon_bs:
    input:
        kinases = rules.prepare_substrate_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_substrate_regulon.output.kinases_phosphatases,
        targets = rules.prepare_substrate_regulon.output.targets,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        mit = rules.stdpi_substrate_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/stdpi_substrate_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule stdpi_substrate_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/stdpi_substrate_regulon/{seed}/{seed}",
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/stdpi_substrate_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/stdpi_substrate_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule stdpi_substrate_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/stdpi_substrate_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/stdpi_substrate_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule stdpi_substrate_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/stdpi_substrate_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/stdpi_substrate_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule stdpi_substrate_siteregulon_generate:
    input:
        network = rules.stdpi_substrate_siteregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/stdpi_substrate_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule stdpi_substrate_proteinregulon_generate:
    input:
        network = rules.stdpi_substrate_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/stdpi_substrate_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule stdpimeta_substrate_regulon_generate:
    input:
        ref = "reference.rds",
        substrate_regulons = [],
        siteregulons = expand("results/{dsid}/stdpi_substrate_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/stdpi_substrate_proteinregulon.rds", dsid=dsids)
    output:
        meta_redundantsite_regulons = "results/stdpimeta_substrate_redundantsite_regulon.rds",
        meta_site_regulons = "results/stdpimeta_substrate_site_regulon.rds",
        meta_protein_regulons = "results/stdpimeta_substrate_protein_regulon.rds",
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate HSM substrate regulon
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
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
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
        iteration = temp("results/{dsid}/hsm_substrate_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_substrate_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/hsm_substrate_regulon/{seed}/{seed}",
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/hsm_substrate_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/hsm_substrate_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule hsm_substrate_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_substrate_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/hsm_substrate_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"

    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule hsm_substrate_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_substrate_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/hsm_substrate_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule hsm_substrate_siteregulon_generate:
    input:
        network = rules.hsm_substrate_siteregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_substrate_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule hsm_substrate_proteinregulon_generate:
    input:
        network = rules.hsm_substrate_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_substrate_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule hsmmeta_substrate_regulon_generate:
    input:
        ref = "reference.rds",
        substrate_regulons = [],
        siteregulons = expand("results/{dsid}/hsm_substrate_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/hsm_substrate_proteinregulon.rds", dsid=dsids)
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate LP substrate regulon
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
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
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
        iteration = temp("results/{dsid}/lp_substrate_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule lp_substrate_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/lp_substrate_regulon/{seed}/{seed}",
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/lp_substrate_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/lp_substrate_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule lp_substrate_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/lp_substrate_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/lp_substrate_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule lp_substrate_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/lp_substrate_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/lp_substrate_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule lp_substrate_siteregulon_generate:
    input:
        network = rules.lp_substrate_siteregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/lp_substrate_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule lp_substrate_proteinregulon_generate:
    input:
        network = rules.lp_substrate_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_substrate_regulon.output.matrix,
        peptides = rules.prepare_substrate_regulon.output.peptides
    output:
        regulon = "results/{dsid}/lp_substrate_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule lpmeta_substrate_regulon_generate:
    input:
        ref = "reference.rds",
        substrate_regulons = [],
        siteregulons = expand("results/{dsid}/lp_substrate_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/lp_substrate_proteinregulon.rds", dsid=dsids)
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate meta substrate regulons
rule meta_substrate_regulon_generate:
    input:
        ref = "reference.rds",
        substrate_regulons = [],
        siteregulons = [expand("results/{dsid}/stdpi_substrate_siteregulon.rds", dsid=dsids), expand("results/{dsid}/hsm_substrate_siteregulon.rds", dsid=dsids), expand("results/{dsid}/lp_substrate_siteregulon.rds", dsid=dsids)],
        proteinregulons = [expand("results/{dsid}/stdpi_substrate_proteinregulon.rds", dsid=dsids), expand("results/{dsid}/hsm_substrate_proteinregulon.rds", dsid=dsids), expand("results/{dsid}/lp_substrate_proteinregulon.rds", dsid=dsids)]
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# prepare DPI activity regulon data
rule prepare_dpi_activity_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        meta_substrate_regulons = rules.stdpimeta_substrate_regulon_generate.output.meta_protein_regulons,
        fasta = "library.fasta"
    output:
        kinases = "results/{dsid}/prepare_dpi_activity_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_dpi_activity_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_dpi_activity_regulon/targets.txt",
        hsm_phosphointeractions = "results/{dsid}/prepare_dpi_activity_regulon/hsm_phosphointeractions.txt",
        pc_phosphointeractions = "results/{dsid}/prepare_dpi_activity_regulon/pc_phosphointeractions.txt",
        lp_phosphointeractions = "results/{dsid}/prepare_dpi_activity_regulon/lp_phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_dpi_activity_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_dpi_activity_regulon/matrix.txt"
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        hsm_threshold = 0,
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20
    threads: 4
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/prepare_activity_regulon.R"

# generate DPI activity regulon
rule dpi_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_dpi_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_dpi_activity_regulon.output.targets,
        matrix = rules.prepare_dpi_activity_regulon.output.matrix
    output:
        mit = "results/{dsid}/dpi_activity_regulon/fwer_computed.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule dpi_activity_regulon_bs:
    input:
        kinases = rules.prepare_dpi_activity_regulon.output.kinases,
        kinases_phosphatases = rules.prepare_dpi_activity_regulon.output.kinases_phosphatases,
        targets = rules.prepare_dpi_activity_regulon.output.targets,
        matrix = rules.prepare_dpi_activity_regulon.output.matrix,
        mit = rules.dpi_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/dpi_activity_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -a {input.kinases} -tg {input.targets} -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule dpi_activity_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/dpi_activity_regulon/{seed}/{seed}",
        peptides = rules.prepare_dpi_activity_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/dpi_activity_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/dpi_activity_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule dpi_activity_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/dpi_activity_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/dpi_activity_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule dpi_activity_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/dpi_activity_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/dpi_activity_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule dpi_activity_siteregulon_generate:
    input:
        network = rules.dpi_activity_siteregulon_consolidate.output.network,
        matrix = rules.prepare_dpi_activity_regulon.output.matrix,
        peptides = rules.prepare_dpi_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/dpi_activity_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule dpi_activity_proteinregulon_generate:
    input:
        network = rules.dpi_activity_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_dpi_activity_regulon.output.matrix,
        peptides = rules.prepare_dpi_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/dpi_activity_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule dpimeta_activity_regulon_generate:
    input:
        ref = rules.stdpimeta_substrate_regulon_generate.input.ref,
        substrate_regulons = rules.stdpimeta_substrate_regulon_generate.output.meta_protein_regulons,
        siteregulons = expand("results/{dsid}/dpi_activity_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/dpi_activity_proteinregulon.rds", dsid=dsids),
        fasta = "library.fasta"
    output:
        meta_redundantsite_regulons = "results/dpimeta_activity_redundantsite_regulon.rds",
        meta_site_regulons = "results/dpimeta_activity_site_regulon.rds",
        meta_protein_regulons = "results/dpimeta_activity_protein_regulon.rds",
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# prepare HSM activity regulon data
rule prepare_hsm_activity_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        meta_substrate_regulons = rules.hsmmeta_substrate_regulon_generate.output.meta_protein_regulons,
        fasta = "library.fasta"
    output:
        kinases = "results/{dsid}/prepare_hsm_activity_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_hsm_activity_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_hsm_activity_regulon/targets.txt",
        hsm_phosphointeractions = "results/{dsid}/prepare_hsm_activity_regulon/hsm_phosphointeractions.txt",
        pc_phosphointeractions = "results/{dsid}/prepare_hsm_activity_regulon/pc_phosphointeractions.txt",
        lp_phosphointeractions = "results/{dsid}/prepare_hsm_activity_regulon/lp_phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_hsm_activity_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_hsm_activity_regulon/matrix.txt"
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        hsm_threshold = 0,
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20
    threads: 4
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/prepare_activity_regulon.R"

# generate HSM activity regulon
rule hsm_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_hsm_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_hsm_activity_regulon.output.hsm_phosphointeractions,
        targets = rules.prepare_hsm_activity_regulon.output.targets,
        matrix = rules.prepare_hsm_activity_regulon.output.matrix
    output:
        mit = "results/{dsid}/hsm_activity_regulon/fwer_computed.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule hsm_activity_regulon_bs:
    input:
        kinases_phosphatases = rules.prepare_hsm_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_hsm_activity_regulon.output.hsm_phosphointeractions,
        targets = rules.prepare_hsm_activity_regulon.output.targets,
        matrix = rules.prepare_hsm_activity_regulon.output.matrix,
        mit = rules.hsm_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/hsm_activity_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule hsm_activity_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/hsm_activity_regulon/{seed}/{seed}",
        peptides = rules.prepare_hsm_activity_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/hsm_activity_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/hsm_activity_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule hsm_activity_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_activity_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/hsm_activity_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule hsm_activity_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/hsm_activity_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/hsm_activity_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule hsm_activity_siteregulon_generate:
    input:
        network = rules.hsm_activity_siteregulon_consolidate.output.network,
        matrix = rules.prepare_hsm_activity_regulon.output.matrix,
        peptides = rules.prepare_hsm_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_activity_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule hsm_activity_proteinregulon_generate:
    input:
        network = rules.hsm_activity_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_hsm_activity_regulon.output.matrix,
        peptides = rules.prepare_hsm_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/hsm_activity_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule hsmmeta_activity_regulon_generate:
    input:
        ref = rules.hsmmeta_substrate_regulon_generate.input.ref,
        substrate_regulons = rules.hsmmeta_substrate_regulon_generate.output.meta_protein_regulons,
        siteregulons = expand("results/{dsid}/hsm_activity_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/hsm_activity_proteinregulon.rds", dsid=dsids),
        fasta = "library.fasta"
    output:
        meta_redundantsite_regulons = "results/hsmmeta_activity_redundantsite_regulon.rds",
        meta_site_regulons = "results/hsmmeta_activity_site_regulon.rds",
        meta_protein_regulons = "results/hsmmeta_activity_protein_regulon.rds",
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# prepare LP activity regulon data
rule prepare_lp_activity_regulon:
    input:
        phospho = "{dsid}_phospho.rds",
        proteo = "{dsid}_proteo.rds",
        meta_substrate_regulons = rules.lpmeta_substrate_regulon_generate.output.meta_protein_regulons,
        fasta = "library.fasta"
    output:
        kinases = "results/{dsid}/prepare_lp_activity_regulon/kinases.txt",
        kinases_phosphatases = "results/{dsid}/prepare_lp_activity_regulon/kinases_phosphatases.txt",
        targets = "results/{dsid}/prepare_lp_activity_regulon/targets.txt",
        hsm_phosphointeractions = "results/{dsid}/prepare_lp_activity_regulon/hsm_phosphointeractions.txt",
        pc_phosphointeractions = "results/{dsid}/prepare_lp_activity_regulon/pc_phosphointeractions.txt",
        lp_phosphointeractions = "results/{dsid}/prepare_lp_activity_regulon/lp_phosphointeractions.txt",
        peptides = "results/{dsid}/prepare_lp_activity_regulon/peptides.txt",
        matrix = "results/{dsid}/prepare_lp_activity_regulon/matrix.txt"
    params:
        minimum_targets = 5,
        maximum_targets = 500,
        adaptive = True,
        fill = "rowmin",
        hsm_threshold = 0,
        ct_correction = True,
        ct_regulators_threshold = 0.05,
        ct_shadow_threshold = 0.05,
        ct_minimum_targets = 5,
        ct_penalty = 20
    threads: 4
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/prepare_activity_regulon.R"

# generate LP activity regulon
rule lp_activity_regulon_mit:
    input:
        kinases_phosphatases = rules.prepare_lp_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_lp_activity_regulon.output.lp_phosphointeractions,
        targets = rules.prepare_lp_activity_regulon.output.targets,
        matrix = rules.prepare_lp_activity_regulon.output.matrix
    output:
        mit = "results/{dsid}/lp_activity_regulon/fwer_computed.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} -o $(dirname {output}) -s 1 -t -j {threads} && touch {output}"

rule lp_activity_regulon_bs:
    input:
        kinases_phosphatases = rules.prepare_lp_activity_regulon.output.kinases_phosphatases,
        phosphointeractions = rules.prepare_lp_activity_regulon.output.lp_phosphointeractions,
        targets = rules.prepare_lp_activity_regulon.output.targets,
        matrix = rules.prepare_lp_activity_regulon.output.matrix,
        mit = rules.lp_activity_regulon_mit.output.mit
    output:
        iteration = temp("results/{dsid}/lp_activity_regulon/{seed}/{seed}")
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "cp $(dirname {input.mit})/fwer_*.txt $(dirname {output}) && "
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -e {input.matrix} -r {input.kinases_phosphatases} -i {input.phosphointeractions} -tg {input.targets} --noDPI -o $(dirname {output}) -s $(basename {output.iteration}) -j {threads} && touch {output.iteration}"

rule lp_activity_proteinregulon_patch:
    input:
        iteration = "results/{dsid}/lp_activity_regulon/{seed}/{seed}",
        peptides = rules.prepare_lp_activity_regulon.output.peptides
    output:
        sitenetwork = "results/{dsid}/lp_activity_siteregulon/bootstrapNetwork_{seed}.txt",
        proteinnetwork = "results/{dsid}/lp_activity_proteinregulon/bootstrapNetwork_{seed}.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/patch_regulon.R"

rule lp_activity_siteregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/lp_activity_siteregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/lp_activity_siteregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule lp_activity_proteinregulon_consolidate:
    input:
        iteration = expand("results/{{dsid}}/lp_activity_proteinregulon/bootstrapNetwork_{seed}.txt", seed=seed)
    output:
        network = "results/{dsid}/lp_activity_proteinregulon/network.txt"
    threads: 1
    container:
        "docker://ghcr.io/califano-lab/vespa.aracne:latest"
    singularity:
        "aracne.simg"
    shell:
        "java -Xmx14G -jar /aracne/dist/aracne.jar -ct 0 -o $(dirname {output}) -c -j {threads} && tar czvf $(dirname {output})/bootstrapNetwork.tgz $(dirname {output})/bootstrapNetwork_*.txt && rm $(dirname {output})/bootstrapNetwork_*.txt"

rule lp_activity_siteregulon_generate:
    input:
        network = rules.lp_activity_siteregulon_consolidate.output.network,
        matrix = rules.prepare_lp_activity_regulon.output.matrix,
        peptides = rules.prepare_lp_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/lp_activity_siteregulon.rds"
    params:
        proteinlevel = False
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule lp_activity_proteinregulon_generate:
    input:
        network = rules.lp_activity_proteinregulon_consolidate.output.network,
        matrix = rules.prepare_lp_activity_regulon.output.matrix,
        peptides = rules.prepare_lp_activity_regulon.output.peptides
    output:
        regulon = "results/{dsid}/lp_activity_proteinregulon.rds"
    params:
        proteinlevel = True
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_regulon.R"

rule lpmeta_activity_regulon_generate:
    input:
        ref = rules.lpmeta_substrate_regulon_generate.input.ref,
        substrate_regulons = rules.lpmeta_substrate_regulon_generate.output.meta_protein_regulons,
        siteregulons = expand("results/{dsid}/lp_activity_siteregulon.rds", dsid=dsids),
        proteinregulons = expand("results/{dsid}/lp_activity_proteinregulon.rds", dsid=dsids),
        fasta = "library.fasta"
    output:
        meta_redundantsite_regulons = "results/lpmeta_activity_redundantsite_regulon.rds",
        meta_site_regulons = "results/lpmeta_activity_site_regulon.rds",
        meta_protein_regulons = "results/lpmeta_activity_protein_regulon.rds",
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"

# generate meta activity regulons
rule meta_activity_regulon_generate:
    input:
        ref = rules.meta_substrate_regulon_generate.input.ref,
        substrate_regulons = rules.meta_substrate_regulon_generate.output.meta_protein_regulons,
        siteregulons = [expand("results/{dsid}/dpi_activity_siteregulon.rds", dsid=dsids), expand("results/{dsid}/hsm_activity_siteregulon.rds", dsid=dsids), expand("results/{dsid}/lp_activity_siteregulon.rds", dsid=dsids)],
        proteinregulons = [expand("results/{dsid}/dpi_activity_proteinregulon.rds", dsid=dsids), expand("results/{dsid}/hsm_activity_proteinregulon.rds", dsid=dsids), expand("results/{dsid}/lp_activity_proteinregulon.rds", dsid=dsids)],
        fasta = "library.fasta"
    output:
        meta_redundantsite_regulons = "results/meta_activity_redundantsite_regulon.rds",
        meta_site_regulons = "results/meta_activity_site_regulon.rds",
        meta_protein_regulons = "results/meta_activity_protein_regulon.rds",
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
    container:
        "docker://ghcr.io/califano-lab/vespa:latest"
    singularity:
        "vespa.simg"
    script:
        "scripts/generate_meta_regulon.R"
