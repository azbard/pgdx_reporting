# Import dependencies
import os
import re
import time
import pandas as pd

from jinja2 import Environment, FileSystemLoader

from pgdx_reporting.reading import TMB_MSI, interp
from pgdx_reporting.igvreport import make_igvreport

# %%
def write_summary(interpreted_report_dict, report_dict):
    """
    Produces html string summary section
    from dictionary of dataframes
    """

    # Get lists of genes and findings
    f = interpreted_report_dict[0]

    # Initialize report section with title
    output_str = (
        "<p><b>SUMMARY:</b><br>Findings with Evidence of Clinical Significance: "
    )

    if len(f["significant_genes"]) == 0:
        output_str += "<b>None</b><br>"
    else:
        output_str += "<b>" + ", ".join(f["significant_genes"]) + "</b><br>"

    # Initialize report section with title

    output_str += (
        "Findings with Evidence of Clinical Significance in Other Tumor Types: "
    )

    if len(f["tier2_genes"]) == 0:
        output_str += "<b>None</b><br>"
    else:
        output_str += "<b>" + ", ".join(f["tier2_genes"]) + "</b><br>"

    # Get TMB result

    sample_sigs = TMB_MSI(report_dict)

    # Add TMB to report
    output_str += "Tumor Mutational Burden: <b>" + sample_sigs["TMB"] + " muts/Mb"

    if float(sample_sigs["TMB"]) >= 16:
        output_str += " (TMB-High)</b><br>"
    else:
        output_str += " (TMB-Low)</b><br>"

    # Get MSI result and add to report

    output_str += "Microsatellite Instability: "

    if sample_sigs["MSI"] == "MSS":
        output_str += "<b>MSS (Microsatellite Stable)</b><br>"
    else:
        output_str += "<b>MSI-High</b><br>"

    output_str += "Findings with Potential Clinical Significance: "

    # List genes with non-significant genes
    if len(f["other_genes"]) == 0:
        output_str += "None</p>"
    else:
        output_str += ", ".join(f["other_genes"]) + "</p><br>"

    return output_str + "<p><b>RESULTS:</b><br>"


def writeSNV_Indel(interpreted_report_dict):

    output_str = "<u>Sequence Variants:</u><br><br>"

    CNVtable = interpreted_report_dict[2].copy()
    CNVtable = CNVtable[CNVtable["Type"] == "Mutation"]

    if CNVtable.empty:
        output_str += "None<br>"
    else:

        for tier, intro in zip(
            ["Tier1", "Tier2", "Tier3"],
            [
                "Findings with Evidence of Clinical Significance: ",
                "Findings with Evidence of Clinical Significance in Other Tumor Types: ",
                "Findings with Potential Clinical Significance: ",
            ],
        ):

            tier_table = CNVtable[CNVtable["Tier"] == tier]

            if not tier_table.empty:
                output_str += (
                    "<b>"
                    + intro
                    + "<br>"
                    + "Gene&nbsp;&nbsp;Variant&nbsp;&nbsp;VAF(%)&nbsp;&nbsp;(Transcript)</b><br>"
                )
                tier_table = tier_table.sort_values("Gene").fillna("n/a")
                for _, row in tier_table.iterrows():
                    if row["Amino Acid Change"] == "n/a":
                        mutation = row["Consequence"]
                    elif row["Amino Acid Change"] == "MET Exon 14 Splice Site Event":
                        mutation = row["Amino Acid Change"]
                    else:
                        mutation = "p." + row["Amino Acid Change"]
                    output_str += "".join(
                        [
                            row["Gene"],
                            "&nbsp;&nbsp;",
                            mutation,
                            "&nbsp;&nbsp;",
                            row["MAF (%)"],
                            "%&nbsp;&nbsp;(",
                            row["Transcript"],
                            ")<br>",
                        ]
                    )
                output_str += "<br>"

    return output_str + "</p>"


def write_amp(interpreted_report_dict):

    CNVtable = interpreted_report_dict[2].copy()
    CNVtable = CNVtable[CNVtable["Type"] == "Amplification"]

    output_str = "<p><u>Copy Number Variants:</u><br>"

    if CNVtable.empty:
        output_str += "None<br>"
    else:
        for tier, intro in zip(
            ["Tier1", "Tier2", "Tier3"],
            [
                "Findings with Evidence of Clinical Significance: ",
                "Findings with Evidence of Clinical Significance in Other Tumor Types: ",
                "Findings with Potential Clinical Significance: ",
            ],
        ):

            tier_table = CNVtable[CNVtable["Tier"] == tier]

            filt = (tier_table["Status"] == "Copy Number Event Detected") | (
                tier_table["Status"] == "Amplified"
            )
            tier_table = tier_table[filt].sort_values("Gene")
            if not tier_table.empty:
                output_str += "<b>" + intro + "</b><br>"
                for _, row in tier_table.iterrows():
                    output_str += row["Gene"] + " Amplified<br>"

    return output_str + "</p>"


def write_trans(interpreted_report_dict):

    # Pulls out Sample Transolcation Analysis table and drop all NA values
    FusionTable = interpreted_report_dict[2].copy()
    FusionTable = FusionTable[FusionTable["Type"] == "Translocation"]

    output_str = "<p><u>Fusions and Translocations:</u><br>"

    if FusionTable.empty:
        output_str += "None<br>"
    else:
        for tier, intro in zip(
            ["Tier1", "Tier2", "Tier3"],
            [
                "Findings with Evidence of Clinical Significance: ",
                "Findings with Evidence of Clinical Significance in Other Tumor Types: ",
                "Findings with Potential Clinical Significance: ",
            ],
        ):

            tier_table = FusionTable[FusionTable["Tier"] == tier]

            if not tier_table.empty:

                output_str += "<b>" + intro + "</b><br>"

                if "Breakpoint" not in tier_table.columns:
                    tier_table["Breakpoint"] = (
                        tier_table["Gene"]
                        + "("
                        + tier_table["Genomic Location"]
                        + ")-"
                        + tier_table["Partner"]
                        + "("
                        + tier_table["Partner Genomic Location"]
                        + ")"
                    )

                tier_table.sort_values("Gene", inplace=True)

                for _, row in tier_table.iterrows():
                    # pulls out left fusion and then right fusion by group
                    chromLocation = re.search(
                        r"\((.*?)\).*?\((.*?)\)", row["Breakpoint"]
                    )

                    # extracts left breakpoint gene and right breakpoint gene by group
                    geneNames = re.search(
                        r"([A-Z0-9]+)\(.*\-(.+)\(.*?", row["Breakpoint"]
                    )

                    # Prints out in text format
                    output_str += "".join(
                        [
                            "A rearrangement was identified involving ",
                            geneNames.group(1),
                            " and ",
                            geneNames.group(2),
                            " with the following genomic breakpoints ",
                            row["Breakpoint"],
                            "<br>",
                        ]
                    )
                # SFA output for comparison:
                # A fusion transcript was identified involving SLC34A2 and ROS1 with the following genomic breakpoints (RET-Chr10:43611047-KIF5B-Chr10:32317346)

    return output_str + "</p>"


def write_qc(report_dict):
    """
    Takes dictionary from parsed PGDx csv and returns tuple:
    bool: pass/fail, str: qc section of report
    """
    output_str = "<p><u>Quality/Run Metrics</u><br>"

    # Check if overall pass/fail in case summary section
    output_str += "Overall: "
    summary_df = report_dict["Case Summary"]
    filt = summary_df["Metric"] == "Overall Case Pass/Fail"

    if (summary_df.loc[filt, "Value"] == "PASS").all():
        output_str += "<b>PASS</b><br>"
        success = True
    else:
        success = False
        output_str += "<b>FAIL</b><br>"

    # output_str += 'Targets with insufficient coverage: '

    # try to get list genes of indeterminate coverage, if fails
    # add 'none'

    return success, output_str + "</p>"


def write_interp_section(interpreted_report_dict):
    output_str = "<p><b>INTERPRETATION:</b><br>"

    f = interpreted_report_dict[0]["findings"]

    if len(f) == 0:
        output_str += "No Findings with Evidence of Clinical Significance"
    else:
        output_str += "Findings with Evidence of Clinical Significance:<br>"

    for finding in f:
        output_str += finding + "<br><br>"

    return output_str + "</p>"


def log(str, log_file):
    """
    takes string, returns it and adds it to a log file
    """
    log_file.write("\n" + str)
    return str


def write_report(
    report_dict,
    interpreted_report_dict,
    fda,
    compatible_versions,
    wiki_dir,
    pathologist_dir,
    batch_dir,
    igv_dir,
    req_dir,
    filename,
    log_file,
    ExonDF,
    mode,
):
    """
    Produces a report file in wiki_directory based on report_dict, fda validated
    targets, and compatible version of PGDx pipeline, and an excel report
    """
    # create a new file named by NA number date and time
    report_datetime = time.strftime("%Y%m%d-%H%M%S", time.localtime())

    case = report_dict["Case Summary"].set_index("Metric").at["Case Name", "Value"]

    # create report file
    file = open(
        os.path.join(wiki_dir, case + " - PGDx Report - " + report_datetime + ".txt"),
        "wb",
    )

    msg = log("Starting report: " + case + " - " + report_datetime, log_file,)
    if mode == "print":
        print(msg)
    elif mode == "yield":
        yield msg

    # Check Version Compatibility
    version = (
        report_dict["Case Summary"]
        .set_index("Metric")
        .loc[["Assay Name", "Assay Version", "Platform Version"]]
        .T.values[0]
        .tolist()
    )

    # If not compatible, add text warning to start of report
    if version in compatible_versions:
        version_message = ""
    else:
        version_message = (
            "<p><b><u>IMPORTANT:</u> Assay/Platform Name or "
            "Version does not match current compatibility and this "
            "report may not be accurate.</b><br> This [Assay Name, Assay "
            "Version, Platform Version]: {}<br>Currently compatible: {}</p>"
        ).format(version, compatible_versions)

    # Check QC pass
    qc = write_qc(report_dict)
    # If testing, add ttype and version numbers at top of report
    ttype_str = ""

    # If passes qc then write report else report = failure message
    if qc[0]:
        # Start writing report

        report_str = "".join(
            [
                version_message,
                ttype_str,
                write_summary(interpreted_report_dict, report_dict),
                writeSNV_Indel(interpreted_report_dict),
                "<br>",
                write_trans(interpreted_report_dict),
                "<br>",
                write_amp(interpreted_report_dict),
                "<br>",
                write_interp_section(interpreted_report_dict),
            ]
        )

        msg = log("Starting html report: " + case + " - " + report_datetime, log_file,)
        if mode == "print":
            print(msg)
        elif mode == "yield":
            yield msg

        html_filename = case + " - PGDx Review - " + report_datetime + ".htm"

        title = "CID#: {}&nbsp;&nbsp;&nbsp;&nbsp;NA#: {}&nbsp;&nbsp;&nbsp;&nbsp;Batch: {}".format(
            *case.split("_")
        )
        subtitle = "Report created: " + time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime()
        )

        cid_num = case.split("_")[0]

        wiki_link = (
            f"https://cid.mghpathology.org/mediawiki/index.php/{cid_num}#tab=Report"
        )

        table = interpreted_report_dict[1].to_html(
            index=False,
            justify="center",
            escape=False,
            table_id="dataframe1",
            classes="table table-sm table-striped table-hover table-responsive w-auto mx-auto",
        )

        igv_filename = f"{case}-igv-report-{report_datetime}.html"
        igv_link = f"../IGV/{igv_filename}"

        try:
            cnv_table = pd.read_csv(
                os.path.join(
                    batch_dir, "TextFiles", case + ".rawfoldchange_ETC-RUO.csv"
                )
            ).to_html(
                index=False,
                justify="center",
                escape=False,
                table_id="dataframe2",
                na_rep="n/a",
                classes="table table-sm table-striped table-hover table-responsive w-auto mx-auto text-center",
            )
        except:
            cnv_table = "No CNV table found."

        try:
            trans_table = pd.read_csv(
                os.path.join(
                    batch_dir, "TextFiles", case + ".translocation_ETC-RUO.csv"
                )
            ).to_html(
                index=False,
                justify="center",
                escape=False,
                na_rep="n/a",
                table_id="dataframe3",
                classes="table table-sm table-striped table-hover table-responsive w-auto mx-auto text-center",
            )
        except:
            trans_table = "No Translocation table found."

        ttype = (
            report_dict["Case Summary"]
            .set_index("Metric")
            .at["Details", "Value"]
            .lower()
            .title()
        )

        sample_sigs = TMB_MSI(report_dict)

        msi = "Microsatellite Analysis: "

        if sample_sigs["MSI"] == "MSI-H":
            msi += "MSI-H"
        else:
            msi += "MSS"

        tmb = f"TMB Muts/Mb: {sample_sigs['TMB']} "

        if float(sample_sigs["TMB"]) >= 16:
            tmb += "High (>=16)"
        else:
            tmb += "Low (<16)"

        loader = FileSystemLoader(req_dir)
        env = Environment(loader=loader)
        template = env.get_template("html_template.html")

        with open(os.path.join(pathologist_dir, html_filename), "w") as f:
            f.write(
                template.render(
                    title=title,
                    subtitle=subtitle,
                    table=table,
                    igv_link=igv_link,
                    cnv_table=cnv_table,
                    trans_table=trans_table,
                    ttype=ttype,
                    tmb=tmb,
                    msi=msi,
                    wiki_link=wiki_link,
                )
            )

        msg = log("Starting igv report: " + case + " - " + report_datetime, log_file,)
        if mode == "print":
            print(msg)
        elif mode == "yield":
            yield msg

        # get BAM
        stub_pat = re.compile(f"{case}.*\.bam$")

        bam = ""

        # search the BAM folder
        for dir in os.listdir(os.path.join(batch_dir, "BAMs")):
            if bool(stub_pat.match(dir)):
                bam = dir
                break

        bam = os.path.join(batch_dir, "BAMs", bam)
        bed1 = os.path.join(req_dir, "hotspots.txt")
        bed2 = os.path.join(req_dir, "bed2")
        hg19_fasta = os.path.join(req_dir, "hg19.fasta")
        html = os.path.join(pathologist_dir, html_filename)
        output_path = os.path.join(igv_dir, igv_filename)

        make_igvreport(ttype, bam, bed1, bed2, hg19_fasta, html, output_path)

    else:
        report_str = (
            "<p><b>RESULTS:</b><br>Targeted next-generation sequencing failed.</p>"
            "<p><b>INTERPRETATION:</b><br>There was insufficient tumor tissue/nucleic"
            " acid quality or quantity for analysis.</p>"
            "<p><b>COMMENT:</b><br>Please submit a request for testing"
            " on another specimen if clinically indicated.</p>"
        )

    file.write(report_str.encode(encoding="utf-8", errors="replace"))
    file.close()

    msg = log(f"Reporting {case} - {report_datetime} complete", log_file)
    if mode == "print":
        print(msg)
    elif mode == "yield":
        yield msg


# %%
