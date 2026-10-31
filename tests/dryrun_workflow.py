#!/usr/bin/env python3
"""Optional real Snakemake EP/ETP DAG construction; no prediction commands."""
import argparse
import configparser
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def main(snakemake):
    results=[]
    with tempfile.TemporaryDirectory(prefix='gene_dag_') as directory:
        base=Path(directory).resolve()
        for mode in ['ep','etp']:
            work=base/mode; work.mkdir()
            genome=work/'genome.fa'; genome.write_text('>chr1\n'+'ACGTacgt'*1500+'\n')
            protein=work/'reference.faa'; protein.write_text('>reference1\nM'+'A'*100+'\n')
            cfg=configparser.ConfigParser(); cfg.read(ROOT/'config/braker.parameters.ini')
            cfg['paths']=dict(samples_file=str(work/'samples.csv'), augustus_config_path=str(work/'augustus_config'),
                              compleasm_download_path=str(work/'db'), busco_download_path=str(work/'db'))
            cfg['SLURM_ARGS']=dict(cpus_per_task='1',mem_of_node='120000',max_runtime='4320')
            image=work/'dummy.sif'; image.touch()
            cfg['containers']={k:str(image) for k in ['braker3_image','agat_image','busco_image']}
            ini=work/'config.ini'
            with ini.open('w') as stream: cfg.write(stream)
            fields='sample_name genome genome_masked protein_fasta bam_files fastq_r1 fastq_r2 sra_ids varus_genus varus_species isoseq_bam isoseq_fastq busco_lineage'.split()
            row=dict.fromkeys(fields,''); row.update(sample_name='example',genome=str(genome),genome_masked=str(genome),protein_fasta=str(protein),busco_lineage='arthropoda_odb12.2')
            if mode=='etp':
                for key in ['fastq_r1','fastq_r2']:
                    path=work/(key+'.fq'); path.write_text('@read1\nACGT\n+\nIIII\n'); row[key]=str(path)
            with (work/'samples.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerow(row)
            env={k:v for k,v in os.environ.items() if not k.startswith('BRAKER4_')}
            env.update(BRAKER4_CONFIG=str(ini),PYTHONNOUSERSITE='1')
            cmd=[snakemake,'--snakefile',str(ROOT/'workflow/Snakefile'),'--directory',str(work),'--cores','1','--dry-run']
            run=subprocess.run(cmd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
            if run.returncode:
                print(run.stdout); raise RuntimeError(mode+' DAG failed')
            results.append(dict(mode=mode,returncode=run.returncode,scope='DAG construction only, no external biological commands'))
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--snakemake',default='snakemake')
    main(parser.parse_args().snakemake)
