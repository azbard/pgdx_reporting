# %%
import subprocess as sp
import pandas as pd
from igv_reports.report import run as create_igv_html

# %%


def subprocess_cmd(command):
    "function to run command line code as part of Python script"
    process = sp.Popen(command, stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)


def make_igvreport(ttype, bam, bed1, bed2, hg19_fasta, html, output_path):
    "Read the genomic coordinates of the additional variants present in the script generated report."
    "The coordinates will be stored in separate lists"
    "output path like: f'../IGV/{case}-igv-report.html'"

    df_overall = pd.read_html(html)
    df = df_overall[0].copy()

    df = df.loc[(df["Type"] == "Mutation") & df["Genomic Location"].notna()]

    if not df.empty:

        df["Genomic Location"] = df["Genomic Location"].str.replace("C", "c")
        df[["chr", "coordinates"]] = df["Genomic Location"].str.split(":", expand=True)
        df[["start", "end"]] = df["coordinates"].str.split("-", expand=True)
        df["chr"] = df["chr"] + ".fa"
        df = df.rename(
            columns={
                "Amino Acid Change": "aminoacid",
                "Tier": "desc",
                "MAF (%)": "maf%",
            }
        )
        df["exon"] = "-"
        df["disease"] = ttype
        df.columns = df.columns.str.strip().str.lower()

        df = df[
            [
                "desc",
                "gene",
                "aminoacid",
                "consequence",
                "maf%",
                "exon",
                "disease",
                "chr",
                "start",
                "end",
            ]
        ]
        tsv = df.to_csv(sep="\t", index=False)
        # print(tsv)

    # "bed1 is the static bed file (Joe and Lauren's). bed2 is a temp text file where we append"
    # "all the regions and use this file as input to bamsnap command"
    pertinent_negatives = ""

    # "clearing the contents of the temp text file for use in next run"
    open(bed2, "w").close()

    with open(bed1, "r", encoding="utf8") as f, open(bed2, "a", encoding="utf8") as u:
        for i in f:
            pertinent_negatives += i
        # print(tsv + pertinent_negatives)
        u.write(tsv + pertinent_negatives)

    # "main igv-reports command" - may want to add info columns on end "chr start end"
    # command = f"create_report '{bed2}' '{hg19_fasta}' --tracks '{bam}' --begin 9 --end 10 --sequence 8 --info-columns desc gene aminoacid consequence maf% exon disease --output '{output_path}'"
    # subprocess_cmd(command)

    # call create_report not as a shell command but as a python function
    s = [
        bed2,
        hg19_fasta,
        "--tracks",
        bam,
        "--begin",
        "9",
        "--end",
        "10",
        "--sequence",
        "8",
        "--info-columns",
        "desc",
        "gene",
        "aminoacid",
        "consequence",
        "maf%",
        "exon",
        "disease",
        "--output",
        output_path,
    ]
    create_igv_html(s)

    # "deleting the contents of the temp text file for use in next run"
    open(bed2, "w").close()


if __name__ == "__main__":
    html = "CID21-11355_NA20-1452_B21-990 - PGDx Review - 20211211-155758.htm"
    # html = 'CID21-11219_NA21-4881_B21-990 - PGDx Review - 20211208-141837.htm'
    # bed1 = 'snapshot_hotspots.txt'
    bed1 = "regions_sample_fa_v2.txt"
    bed2 = "sample.txt"
    bam = "../../../igv-reports/bamfile.bam"
    hg19_fasta = "../../../igv-reports/hg19.fasta"
    # bam = 'test_alignments.bam'
    # vaf = "0.02"
    output_path = "outsnaps104.html"

    make_igvreport(bam, bed1, bed2, hg19_fasta, html, output_path)
