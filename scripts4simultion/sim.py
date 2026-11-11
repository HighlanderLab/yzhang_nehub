import msprime
import numpy as np
import tskit
import pandas as pd
import stat_plot
import concurrent.futures
import os
from pathlib import Path
import params
import demes
import io
import random

DEFAULTS = {
    'length': 80_000_000,  # Length of the genomic region in base pairs
    'rho': 9.5e-9, #21.5 Morgan (M) (0.95 cM/Mb) and the estimated female and male maps were 23.6 M (1.04 cM/Mb) and 19.5 M (0.86 cM/Mb), <Genetic variation in recombination rate in the pig>
    'mutation_rate': 6.3e-9, # Mutation rate: 4.32e-9 mutations per site per generation <Evolution of the germline mutation rate across vertebrates>
    'num_replicates': 10,
    'nsample': 200,
    'nsample1': 100, # one of population's sample size 
    'name': 'con_split_mig',
    'path': '../../data_simulation/CON_SPLIT_MIG_cor'
}
args = params.parse_args_with_defaults(DEFAULTS)

print("=== Parameters ===")
for key, value in vars(args).items():
    print(f"  {key}: {value}")
print("==================")

print(f"#### load {args.name}.yaml demography")
graph = demes.load(f"{args.name}.yaml")
demography = msprime.Demography.from_demes(graph)
# grouwth rate: r = -ln(Nt/N0)/t

def write_fasta(ts,outdir):
    """Write a random ancestral/reference FASTA.  Chromosome id = 'chr1'.

    Both the tskit, Relate, PSMC, and stairway modules need a FASTA whose
    sequence id matches the CHROM column in the VCF ('chr1').
    """
    print(f"  DEBUG outdir type: {type(outdir)}, value: {outdir}")
    print(f"  DEBUG ts type: {type(ts)}")
    print(f"  DEBUG sequence_length: {ts.sequence_length}")
    print(f"  DEBUG num_sites: {ts.num_sites}")
    print(f"  DEBUG num_mutations: {ts.num_mutations}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        s = [random.choices(
                ["A", "T", "C", "G"],
                [0.225, 0.225, 0.275, 0.275],
                k=1)[0]
             for i in range(int(ts.sequence_length))]
        print(f"  DEBUG sequence generated: {len(s)} bases")
    except Exception as e:
        print(f"  ERROR generating sequence: {e}")
        raise
    for i in ts.sites():
        s[int(i.position)]=i.ancestral_state
    # write reference genome, non-variant sites will be written randomly (ATCG)
    out = outdir / "simulation.fa"
    f=open(out,"w")
    f.write(f">chr1\n")
    for i in range(len(s)):
        if (i%60==0) and (i!=0):
            f.write("\n")
        f.write(s[i])
    f.write("\n")
    f.close()
    print(f"  FASTA: {outdir}/simulation.fa")
def write_vcf(mts, outdir):
    """Write bgzipped, tabix-indexed VCF. File: test_chr1.vcf.gz"""
    print("#### vcf ####")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw_vcf = outdir / "sim_chr1.vcf"
    print(f"  DEBUG raw_vcf: {raw_vcf}")
    try:
        buf = io.StringIO()
        mts.write_vcf(buf, contig_id="chr1",position_transform=lambda x: np.array(x) + 1)
        content = buf.getvalue()
    except Exception as e:
        print(f"  ERROR write_vcf failed: {e}")
        raise
    
    # check content is not empty
    if not content:
        raise ValueError("VCF content is empty — check mts has mutations")
    print(f"  VCF content: {len(content)} bytes, {mts.num_mutations} mutations")
    
    raw_vcf.write_text(content)
    # verify file was written
    if raw_vcf.stat().st_size == 0:
        raise ValueError(f"VCF file is empty: {raw_vcf}")
def write_poplabel(ts, outdir):
    """Write TSV population label using actual population names from ts."""
    rows = []
    print(f"  DEBUG num_individuals: {ts.num_individuals}")
    print(f"  DEBUG num_samples: {ts.num_samples}")
    for individual in ts.individuals():
        node = ts.node(individual.nodes[0])
        pop_name = ts.population(node.population).metadata.get("name", f"pop{node.population}")
        rows.append((
            f"tsk_{individual.id}",
            pop_name,
            pop_name,
            "NA"
        ))
    df = pd.DataFrame(rows, columns=["sample_id", "population", "group", "sex"])
    out = Path(outdir) / "popLabel.txt"
    df.to_csv(out, sep="\t", index=False)
    print(f"  popLabel: {out} ({len(df)} samples)")
def write_uniform_genetic_map(rho, length, chrom, outdir):
    positions = [0, int(length)]
    rates = [rho * 1e8, 0]  # cM/Mb
    distances = [0, rho * length * 100]  # cM

    df = pd.DataFrame({
        "pos": positions,
        "rate(cM/Mb)": rates,
        "Genetic_Map(cM)": distances
    })
    out = os.path.join(outdir, f"chr{chrom}.genetic_map.txt")
    df.to_csv(out, sep=" ", index=False)
def check_vcf_fasta_consistency(vcf_path, fasta_path):
    """Check VCF and FASTA are consistent."""
    print("#### checking VCF vs FASTA ####")
    
    # ── read FASTA ────────────────────────────────────────────────────────────
    fasta_seq = {}
    current_chr = None
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_chr = line[1:]
                fasta_seq[current_chr] = []
            else:
                fasta_seq[current_chr].extend(list(line))
    
    fasta_seq = {k: "".join(v) for k, v in fasta_seq.items()}
    print(f"  FASTA chromosomes: {list(fasta_seq.keys())}")
    for chrom, seq in fasta_seq.items():
        print(f"  FASTA {chrom}: {len(seq)} bp")

    # ── read VCF ──────────────────────────────────────────────────────────────
    mismatches = []
    matches = 0
    total = 0

    import gzip
    open_func = gzip.open if str(vcf_path).endswith(".gz") else open

    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos   = int(fields[1])    # 1-based
            ref   = fields[3]
            alt   = fields[4]

            # check chrom exists in FASTA
            if chrom not in fasta_seq:
                print(f"  ERROR: {chrom} not in FASTA")
                continue

            # check REF matches FASTA (convert 1-based to 0-based)
            fasta_base = fasta_seq[chrom][pos - 1].upper()
            total += 1

            if fasta_base != ref.upper():
                mismatches.append({
                    "chrom": chrom,
                    "pos":   pos,
                    "vcf_ref": ref,
                    "fasta_base": fasta_base,
                    "alt": alt
                })
            else:
                matches += 1

    # ── report ────────────────────────────────────────────────────────────────
    print(f"  Total variants: {total}")
    print(f"  REF matches:    {matches}")
    print(f"  REF mismatches: {len(mismatches)}")

    if mismatches:
        print(f"  First 5 mismatches:")
        for m in mismatches[:5]:
            print(f"    {m['chrom']}:{m['pos']} VCF_REF={m['vcf_ref']} FASTA={m['fasta_base']} ALT={m['alt']}")
    else:
        print("  ✅ All REF alleles match FASTA")

    return len(mismatches) == 0

def process_replicate(replicate_index, ts, mutation_rate,path):
    ts_mutated = msprime.sim_mutations(ts, rate=mutation_rate)
    Path(os.path.join(path,f"replicate{replicate_index}")).mkdir(parents=True, exist_ok=True)
    ts_mutated.dump(os.path.join(path,f"replicate{replicate_index}/sim_chr1.trees"))
    write_vcf(ts_mutated,os.path.join(path,f"replicate{replicate_index}"))
    write_fasta(ts_mutated,os.path.join(path,f"replicate{replicate_index}"))
    write_poplabel(ts_mutated,os.path.join(path,f"replicate{replicate_index}"))
    check_vcf_fasta_consistency(
        vcf_path   = os.path.join(path,f"replicate{replicate_index}/sim_chr1.vcf"),
        fasta_path = os.path.join(path,f"replicate{replicate_index}/simulation.fa")
    )


def Multi_thread(path, name, nsamp, nsample1,num_replicates, mutation_rate,demography,rho,length):
    # Simulate the tree sequence
    if "pop2" in [pop.name for pop in demography.populations]:
        nsample2=nsamp-nsample1
        sample_list={"pop1":nsample1, "pop2": nsample2} # Specify number of samples and ploidy
    else:
        sample_list={"pop1":nsample1}  # Specify number of samples and ploidy
    
    
    replicates = msprime.sim_ancestry(
        samples=sample_list,
        demography=demography,
        recombination_rate=rho,
        ploidy=2,
        # correct the LD bias of coalescence simulation
        model=[
            msprime.DiscreteTimeWrightFisher(duration=200),
            msprime.StandardCoalescent(),
        ],
        sequence_length=length,
        num_replicates=num_replicates,  # Set a random seed for reproducibility
    )
    Path(path).mkdir(parents=True, exist_ok=True) 
    write_uniform_genetic_map(rho, length, chrom=1, outdir=path)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_replicate, replicate_index, ts, mutation_rate,path)
            for replicate_index, ts in enumerate(replicates)
        ]
        
        concurrent.futures.wait(futures)

Multi_thread(args.path, args.name, args.nsample,args.nsample1, args.num_replicates, args.mutation_rate,demography,args.rho,args.length)
# stat_plot.Multi_population(args.path,args.name,args.nsample,args.nsample1, args.num_replicates, args.mutation_rate)



