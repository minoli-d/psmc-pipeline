import pandas as pd

# setup
BASEDIR = '/global/scratch/users/nicolas931010/sv_diversity'
SPECIES_TSV = '/global/scratch/users/minoli/psmc_pipeline/final_species_list.tsv'

# read
species_df = pd.read_csv(SPECIES_TSV, sep='\t')
species_info = species_df.set_index('scientific_name')
SPECIES_LIST = species_info.index.tolist()

# paths
def get_fasta(species):
    return BASEDIR + '/assemblies/' + species_info.loc[species, 'accession'] + '/' + species_info.loc[species, 'fasta']

def get_vcf(species):
    return BASEDIR + '/paftools/' + species_info.loc[species, 'biosample'] + '.for_mutyper.vcf'

def get_bed(species):
    return BASEDIR + '/paftools/' + species_info.loc[species, 'biosample'] + '.callable.bed'

def get_mutation_rate(species):
    return species_info.loc[species, 'mu']

def get_generation_time(species):
    return species_info.loc[species, 'gen']

def get_group_mu(species):
    return species_info.loc[species, 'group_mu']


# rules
rule all:
    input:
        expand("run_outputs/{species}/{species}.psmc", species=SPECIES_LIST),
        expand("run_outputs/{species}/{species}.psmcfa", species=SPECIES_LIST),
        expand("scaled_outputs/{species}_scaled.txt", species=SPECIES_LIST),
        expand("relative_outputs/{species}_scaled.txt", species=SPECIES_LIST)

rule filter_vcf:
    input:
        vcf = lambda wc: get_vcf(wc.species)
    output:
        snp = temp("run_outputs/{species}/{species}.snp")
    resources:
        threads=1,
        runtime="15m",
        mem_mb=2000
    conda:
        "envs/psmc.yaml"
    shell:
        """
        mkdir -p run_outputs/{wildcards.species}
        bcftools view -m2 -M2 -v snps {input.vcf} | ~/psmc/utils/vcf2snp.pl - > {output.snp}
        """

rule mutfa:
    input:
        ref = lambda wc: get_fasta(wc.species),
        snp = "run_outputs/{species}/{species}.snp"
    output:
        pseudo = temp("run_outputs/{species}/{species}.pseudo.fa")
    resources:
        threads=1,
        runtime="30m",
        mem_mb=2000
    conda:
        "envs/psmc.yaml"
    shell:
        """
        seqtk mutfa {input.ref} {input.snp} > {output.pseudo}
        """

rule mask:
    input:
        pseudo = "run_outputs/{species}/{species}.pseudo.fa",
        bed = lambda wc: get_bed(wc.species)
    output:
        masked = temp("run_outputs/{species}/{species}.masked.fa")
    resources:
        threads=1,
        runtime="30m",
        mem_mb=2000
    conda:
        "envs/psmc.yaml"
    shell:
        """
        seqtk seq -cM {input.bed} -l80 {input.pseudo} > {output.masked}
        """

rule fq2psmcfa:
    input:
        masked = "run_outputs/{species}/{species}.masked.fa"
    output:
        psmcfa = temp("run_outputs/{species}/{species}.psmcfa")
    resources:
        threads=1,
        mem_mb=1000,
        runtime="1h"    
    conda:
        "envs/psmc.yaml"
    shell:
        """
        ~/psmc/utils/fq2psmcfa {input.masked} > {output.psmcfa}
        """

rule run_psmc:
    input:
        psmcfa = "run_outputs/{species}/{species}.psmcfa"
    output:
        psmc = "run_outputs/{species}/{species}.psmc"
    resources:
        threads=4,
        runtime="5h",
        mem_mb=16000
    conda:
        "envs/psmc.yaml"
    shell:
        """
        psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}
        """

rule scale_psmc:
    input:
        psmc = "run_outputs/{species}/{species}.psmc"
    output:
        txt = "scaled_outputs/{species}_scaled.txt"
    resources:
        threads=1,
        runtime="15m",
        mem_mb=1000
    params:
        species = lambda wc: wc.species,
        mu = lambda wc: float(get_mutation_rate(wc.species)),
        gen = lambda wc: get_generation_time(wc.species)
    conda:
        "envs/psmc.yaml"
    shell:
        """
        mkdir -p scaled_outputs
        ~/psmc/utils/psmc_plot.pl -R -u {params.mu} -g {params.gen} {params.species} {input.psmc}
        mv {params.species}.0.txt {output.txt}
        rm *.gp
        rm *.par
        """


rule scale_psmc_relative:
    input:
        psmc = "run_outputs/{species}/{species}.psmc"
    output:
        txt = "relative_outputs/{species}_scaled.txt"
    resources:
        threads=1,
        runtime="15m",
        mem_mb=1000
    params:
        species = lambda wc: wc.species,
        mu = lambda wc: float(get_group_mu(wc.species)),
        gen = 1
    conda:
        "envs/psmc.yaml"
    shell:
        """
        mkdir -p relative_outputs
        ~/psmc/utils/psmc_plot.pl -R -u {params.mu} -g {params.gen} {params.species} {input.psmc}
        mv {params.species}.0.txt {output.txt}
        rm *.gp
        rm *.par
        """
