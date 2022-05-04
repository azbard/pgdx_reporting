import pandas as pd
from pandas.api.types import CategoricalDtype
import re
import os

# #Takes in files that have exon information from UCSC
# # (can take it regardless of genomic version) and returns a succinct dataframe
def createExonTable(file):
    # Reads the file and gives headers
    ExonDF = pd.read_csv(
        file,
        sep="\t",
        names=[
            "Chr",
            "ExonStart",
            "ExonEnd",
            "ccds",
            "unknown",
            "strand",
        ],
    )

    # splits columns to specify start and end location
    ExonDF[
        [
            "CCDS",
            "ExonRepeat",
            "Exon",
            "ExonRight",
            "Chrom",
            "Bases",
            "Direction",
        ]
    ] = ExonDF["ccds"].str.split("_", expand=True)

    # drops unnecessary columns
    doNotNeed = [
        "ccds",
        "unknown",
        "ExonRepeat",
        "ExonRight",
        "Bases",
        "Direction",
        "strand",
        "Chrom",
    ]

    ExonDF.drop(doNotNeed, axis=1, inplace=True)

    return ExonDF


def parse_csv(fullpath: str):
    """
    Parses tables in csv and returns dictionary where
    keys are table names, values are dataframes.
    """
    # Read in the file - we have to give column titles thru 'h'
    # because tables higher in the csv have fewer columns and
    # we need to make sure we import them all
    PGDx_xl_rprt_df = pd.read_csv(
        fullpath,
        names=["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        index_col=False,
    )

    # Create a column which will hold the sub-dataframes' title
    # fill with true if column 'a' contains something in square brackets
    PGDx_xl_rprt_df["title"] = PGDx_xl_rprt_df["a"].str.contains(
        "\[.+\]"
    )

    # replace title column with column 'a' if true and leave empty
    # otherwise and forward fill with title
    PGDx_xl_rprt_df["title"] = (
        PGDx_xl_rprt_df[["title", "a"]]
        .apply(lambda x: x[1] if x[0] == True else pd.NA, axis=1)
        .ffill()
    )

    # Group by 'title'
    grouped = PGDx_xl_rprt_df.groupby("title")

    # Initialize dictionary of tables
    report_dict = dict()

    for title in PGDx_xl_rprt_df.title.unique():
        # Select group of interest and drop title columns
        df = grouped.get_group(title).drop(columns=["title"])

        # Replace columns with second line
        df.columns = df.iloc[1, :]

        # Drop first two rows and any empty rows or columns
        df = (
            df.iloc[2:, :]
            .reset_index(drop=True)
            .rename_axis(None, axis=1)
            .dropna(axis=1, how="all")
            .dropna(axis=0, how="all")
        )

        # add to dictionary (stripped of square brackets)
        report_dict[title[1:-1]] = df

    return report_dict


def addExonData(report_dict, ExonDF):
    """
    Takes in the parsed report and Exon dataframe and merges them
    to add exon info to the report sequence variants table
    """
    # drops any null values in the table
    MutTable = report_dict["Sample Sequence Mutation Analysis"]

    # Splits column to seperate out Chromsome, Start, End of the mutation
    MutTable[["Chromsome", "Location"]] = MutTable[
        "Genomic Location"
    ].str.split(":", expand=True)

    MutTable[["MutStart", "MutEnd"]] = MutTable["Location"].str.split(
        "-", expand=True
    )
    MutTable.drop("Location", axis=1, inplace=True)

    # Merges with Exon table
    MergedDF = pd.merge(
        MutTable,
        ExonDF,
        how="left",
        left_on="Transcript",
        right_on="CCDS",
    )

    # checks to see if the mutation falls within the Exon or if not enough data
    FilteredDF = MergedDF[
        (pd.isnull(MergedDF["Chr"]))
        | (
            (
                pd.to_numeric(MergedDF["MutStart"])
                >= MergedDF["ExonStart"]
            )
            & (
                pd.to_numeric(MergedDF["MutEnd"])
                <= MergedDF["ExonEnd"]
            )
        )
    ].copy()

    # drops duplicate columns
    FilteredDF.drop(
        [
            "Chromsome",
            "MutStart",
            "MutEnd",
            "CCDS",
            "Chr",
            "ExonStart",
            "ExonEnd",
        ],
        axis=1,
        inplace=True,
    )

    return FilteredDF


def TMB_MSI(report_dict):
    """
    Takes parsed csv and produces dictionary
    of TMB and MSI results e.g. {'TMB': '4.6', 'MSI': 'MSS'}
    """
    if "Sample Genomic Signatures" in report_dict.keys():
        sample_sigs_df = report_dict[
            "Sample Genomic Signatures"
        ].set_index("Signature")

        output_dict = dict()

        # Different versions of the assay/platform have different
        # names for TMB so try to find the most complex first
        try:
            output_dict["TMB"] = sample_sigs_df.at[
                "TMB Muts/Mb (Sequenced)", "Status/Score"
            ]
        except:
            output_dict["TMB"] = sample_sigs_df.at[
                "TMB Muts/Mb", "Status/Score"
            ]

        output_dict["MSI"] = sample_sigs_df.at[
            "Microsatellite Analysis", "Status/Score"
        ]
    else:
        output_dict = {"TMB": None, "MSI": None}
    return output_dict


# helper function
def astype_inplace(df: pd.DataFrame, dct: dict):
    """
    takes dataframe and dict of column names: desired dtypes
    returns dataframe with columns in desired format
    """
    df.loc[:, list(dct.keys())] = df.astype(dct)[list(dct.keys())]
    return df


def split_mut_pattern_if_std(x: str):
    """
    Looks for standard missense mutation pattern and returns a dict of parts if found
    """
    # pattern for a standard muation e.g. V600E i.e. no del, ins, etc.
    std_mut_pattern = re.compile(
        r"^(?P<first>[A-Z])(?P<codon>[0-9]{1,4})(?P<last>[A-Z,\*\?])$"
    )

    codon_pattern = re.search(std_mut_pattern, str(x))

    if bool(codon_pattern):
        # return int(codon_pattern.group(2))
        return codon_pattern.groupdict()
    else:
        return {"first": pd.NA, "codon": pd.NA, "last": pd.NA}


def check_K3326_and_up(x):
    """
    Takes a row in the mutation table and checks if its mutation is of the form
    K and then a number >= 3326, if so returns false, else true
    """
    parsed_mutation = split_mut_pattern_if_std(x["Amino Acid Change"])
    if pd.isnull(parsed_mutation["first"]):
        out = True
    elif (parsed_mutation["first"] == "K") & (
        int(parsed_mutation["codon"]) >= 3326
    ):
        out = False
    else:
        out = True

    return out


def add_link_columns(df0, batch_dir, filename):
    """
    Takes in the dataframe generated as an output from the interp fn and adds
    columns with links
    """
    df = df0.reset_index(drop=True).copy()

    # find bam file by stipping PGDx report file of ending
    file_stub = filename.replace("CCR_ETC-RUO.csv", "")

    # add .bam on end
    stub_pat = re.compile(f"{file_stub}.*\.bam$")

    bam_file = ""

    # search the BAM folder
    for dir in os.listdir(os.path.join(batch_dir, "BAMs")):
        if bool(stub_pat.match(dir)):
            bam_file = dir
            break

    macfile = f"/Volumes/PGDx_v1$/ElioConnect_Output/{os.path.basename(batch_dir)}/BAMs/{bam_file}"
    windowsfile = f"\\\\smbgpb\pgdx_v1\ElioConnect_Output\{os.path.basename(batch_dir)}\BAMs\{bam_file}"

    # Splitting columns by delimiter to obtain Chromosome, start and end point
    df[["Chromosome", "Start-End Point"]] = df[
        "Genomic Location"
    ].str.split(":", expand=True)
    df[["Start Point", "End Point"]] = df[
        "Start-End Point"
    ].str.split("-", expand=True)

    # Filling NaN with empty string for string concatention later
    df.fillna(" ", inplace=True)

    # Making empty list to be filled
    google_list = []
    clinvar_list = []
    gnomad_list = []
    mycancegenome_list = []
    WindowsIGV_list = []
    MacIGV_list = []

    for i in range(df.shape[0]):
        gene_name = df["Gene"].loc[i]
        variant = df["Amino Acid Change"].loc[i]
        chr_name = df["Chromosome"].loc[i]
        start = df["Start Point"].loc[i]
        end = df["End Point"].loc[i]
        type_var = df["Type"].loc[i]

        # Populating the Google, ClinVar, gnomAD and Mycancergenome columns
        clinvar_value = " "
        mycancegenome_value = " "
        gnomad_value = " "
        WindowsIGV_value = " "
        MacIGV_value = " "

        if df.at[i, "Type"] == "MSI" or df.at[i, "Type"] == "TMB-H":
            google_value = " "
        elif type_var == "Translocation":
            google_value = f'<a href="https://www.google.com/search?q={gene_name}+fusion" target="_blank" rel="noopener noreferrer">Google</a>'
        elif type_var == "Amplification":
            google_value = f'<a href="https://www.google.com/search?q={gene_name}+amplification" target="_blank" rel="noopener noreferrer">Google</a>'
            WindowsIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={windowsfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Win</a>'
            MacIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={macfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Mac</a>'
        elif variant == " ":
            google_value = f'<a href="https://www.google.com/search?q={gene_name}+{type_var}" target="_blank" rel="noopener noreferrer">Google</a>'
            clinvar_value = f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term={gene_name}[sym]{type_var}" target="_blank" rel="noopener noreferrer">ClinVar</a>'
            mycancegenome_value = f'<a href="https://www.mycancergenome.org/content/alteration/{gene_name}-{type_var}" target="_blank" rel="noopener noreferrer">MCG</a>'
            gnomad_value = f'<a href="https://gnomad.broadinstitute.org/region/{chr_name}-{start}-{end}" target="_blank" rel="noopener noreferrer">gnomAD</a>'
            WindowsIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={windowsfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Win</a>'
            MacIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={macfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Mac</a>'
        else:
            google_value = f'<a href="https://www.google.com/search?q={gene_name}+{variant}" target="_blank" rel="noopener noreferrer">Google</a>'
            clinvar_value = f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term={gene_name}[sym]{variant}" target="_blank" rel="noopener noreferrer">ClinVar</a>'
            mycancegenome_value = f'<a href="https://www.mycancergenome.org/content/alteration/{gene_name}-{variant}" target="_blank" rel="noopener noreferrer">MCG</a>'
            gnomad_value = f'<a href="https://gnomad.broadinstitute.org/region/{chr_name}-{start}-{end}" target="_blank" rel="noopener noreferrer">gnomAD</a>'
            WindowsIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={windowsfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Win</a>'
            MacIGV_value = f'<a href="http://localhost:60151/load?merge=false&amp;genome=hg19&amp;file={macfile}&amp;locus={chr_name}%3A{str(int(start) - 1)}" >IGV_Mac</a>'
        google_list.append(google_value)
        clinvar_list.append(clinvar_value)
        gnomad_list.append(gnomad_value)
        mycancegenome_list.append(mycancegenome_value)
        WindowsIGV_list.append(WindowsIGV_value)
        MacIGV_list.append(MacIGV_value)

    df["Google"] = google_list
    df["ClinVar"] = clinvar_list
    df["gnomAD"] = gnomad_list
    df["MCG"] = mycancegenome_list
    df["IGV: Win"] = WindowsIGV_list
    df["IGV: Mac"] = MacIGV_list

    if (
        ("Tier" in df.columns)
        and ("Type" in df.columns)
        and ("Gene" in df.columns)
    ):
        for col in ["Gene", "Type", "Tier"]:
            mid = df[col]
            df.drop(labels=[col], axis=1, inplace=True)
            df.insert(0, col, mid)

    amps = df["Type"] == "Amplification"
    if not df[amps].empty:
        df.loc[amps, "Consequence"] = "Amplified"

    trans = df["Type"] == "Translocation"
    if not df[trans].empty:
        df.loc[trans, "Genomic Location"] = df[trans]["Breakpoint"]
        df.loc[trans, "Gene"] = (
            df[trans]["Partner"] + "::" + df[trans]["Gene"]
        )

    max_columns_to_drop = set(
        [
            "Exon",
            "Mutant",
            "Reference",
            "Transcript",
            "Chromosome",
            "Start-End Point",
            "Start Point",
            "End Point",
            "Status",
            "Breakpoint",
            "Partner",
        ]
    )

    columns_to_drop = list(set(df.columns) & max_columns_to_drop)

    if len(columns_to_drop) != 0:
        df.drop(labels=columns_to_drop, axis=1, inplace=True)

    return df


def METexon14skipping(report_dict):
    """
    Takes report_dict and returns dataframe either empty
    or with only MET Exon14 skipping rows
    """
    # drops any null values in the table
    MutTable = report_dict["Sample Sequence Mutation Analysis"]
    MutTable = MutTable[MutTable["Gene"].notna()].copy()

    # pulls all MET mutations from table
    METtable = MutTable.loc[
        (MutTable["Gene"] == "MET")
        # these requirements discarded per LR 10/21/21
        # & (MutTable['Consequence'] == 'Splice site donor')
        # & (MutTable['Amino Acid Change'].isnull())
    ].copy()

    if METtable.empty:
        return METtable
    else:
        # Splits column to seperate out Chromsome, Start, End of the mutation
        METtable[["Chromsome", "Location"]] = METtable[
            "Genomic Location"
        ].str.split(":", expand=True)
        METtable[["MutStart", "MutEnd"]] = METtable[
            "Location"
        ].str.split("-", expand=True)
        METtable.drop("Location", axis=1, inplace=True)

        # changes columns to integers instead of strings
        astype_inplace(METtable, {"MutStart": int, "MutEnd": int})

        # Iterates through all MET mutations to find FDA approved ones
        for index, row in METtable.iterrows():
            # Defines the ranges of Mutations as two positions outside of MET exon14 > 7:116411903-116412043
            MutRange = range((row["MutStart"]), (row["MutEnd"] + 1))
            #  Note: this range changed 10/21/21
            LowerRange = range(116411900, 116411904)
            UpperRange = range(116412043, 116412047)

            # Creates a list of those that overlap between mutants and +/- 2 positions of Exon 14
            WithinLower = list(set(MutRange) & set(LowerRange))
            WithinUpper = list(set(MutRange) & set(UpperRange))

            # Checks if anything was in those lists and marks them on the dataframe
            if not WithinLower and not WithinUpper:
                METtable.at[index, "FDA"] = False
            else:
                METtable.at[index, "FDA"] = True

        # Makes table of FDA approved targets previously identified and drops all columns used in calculations
        FDAapprovedDF = METtable.loc[METtable["FDA"] == True].copy()
        FDAapprovedDF.drop(
            ["FDA", "Chromsome", "MutStart", "MutEnd"],
            axis=1,
            inplace=True,
        )

        return FDAapprovedDF


def interp(report_dict, fda, batch_dir, ExonDF, filename):
    """
    Takes parsed PGDx file and list of fda validated findings
    and returns a tuple:
        1) dictionary of lists of strings. (output_dict)
            Keys are:
            - significant genes
            - second tier genes
            - other genes
            - findings
        2) an dataframe for excel file creation of 3)
            but with links added and fewer columns
        3) a dataframe of all mutations with
            - a tier [Tier1-Tier3] and
            - a type [mutation, translocation, amplification, MSI, TMB]
        4) a copy of the report_dict but
            - with tiers and,
            - if Tier ==1, a tier_variant_id for later linking to a comment

    """
    # fda_version = fda.at[0, "version"]

    # Intitialize
    output_dict = {
        "significant_genes": [],
        "tier2_genes": [],
        "other_genes": [],
        "findings": [],
    }

    output_excel = {
        "significant_genes": (pd.DataFrame({"Tier": ["Tier1"]})),
        "tier2_genes": (pd.DataFrame({"Tier": ["Tier2"]})),
        "other_genes": (pd.DataFrame({"Tier": ["Tier3"]})),
    }

    # # create a dict to fill with tier info
    tiered_dict = report_dict.copy()

    # empty result tables in the new dict
    for tbl in [
        "Sample Genomic Signatures",
        "Sample Sequence Mutation Analysis",
        "Sample Amplification Analysis",
        "Sample Translocation Analysis",
    ]:
        if tbl in tiered_dict.keys():
            tiered_dict[tbl] = tiered_dict[tbl][0:0]
            len_cols = len(tiered_dict[tbl].columns)
            tiered_dict[tbl].insert(len_cols, "tier", [])
            tiered_dict[tbl].insert(len_cols + 1, "tiered_var_id", [])

    # Get MSI results add to findings if MSI

    sample_sigs = TMB_MSI(report_dict)
    sigs = report_dict["Sample Genomic Signatures"]
    filt = sigs["Signature"] == "Microsatellite Analysis"
    row = sigs.loc[filt].loc[0].to_list()

    if sample_sigs["MSI"] == "MSI-H":
        fda_by_type = fda.set_index("mut_type", drop=True)
        fda_comment = fda_by_type.at["MSI", "Interpretation"]
        fda_code = fda_by_type.at["MSI", "tiered_var_id"]
        output_dict["findings"] += [fda_comment]
        msi = "MSI"
        output_excel["significant_genes"] = output_excel[
            "significant_genes"
        ].append({"Type": "MSI"}, ignore_index=True)
        row.extend(["1", str(fda_code)])
    else:
        msi = "MSS"
        row.extend(["0", "0"])

    tiered_dict["Sample Genomic Signatures"].loc[0] = row

    filt = sigs["Signature"].str.startswith("TMB")
    rows = sigs.loc[filt]

    if float(sample_sigs["TMB"]) >= 16:
        fda_by_type = fda.set_index("mut_type", drop=True)
        fda_comment = fda_by_type.at["TMB", "Interpretation"]
        fda_code = fda_by_type.at["TMB", "tiered_var_id"]
        output_dict["findings"] += [fda_comment]
        tmb = "H"
        output_excel["significant_genes"] = output_excel[
            "significant_genes"
        ].append({"Type": "TMB-H"}, ignore_index=True)
        tier_info = pd.DataFrame(
            {"tier": ["0", "1"], "tiered_var_id": ["0", fda_code]},
            index=[1, 2],
        )
    else:
        tmb = "L"
        tier_info = pd.DataFrame(
            {"tier": ["0", "0"], "tiered_var_id": ["0", "0"]},
            index=[1, 2],
        )

    to_add = pd.concat(
        [
            rows,
            tier_info,
        ],
        axis=1,
    )
    tiered_dict["Sample Genomic Signatures"] = tiered_dict[
        "Sample Genomic Signatures"
    ].append(
        to_add,
        ignore_index=True,
    )

    # Get ttype
    ttype = (
        report_dict["Case Summary"]
        .set_index("Metric")
        .at["Details", "Value"]
        .upper()
    )

    # Initalize
    list_o_lungs = [
        "LUNG NSCLC",
        "LUNG NSCLC-AD",
        "LUNG NSCLC-LCNEC",
        "LUNG NSCLC-SQ",
    ]
    fda_filt = (fda.ttype == ttype) | (fda.ttype == "Solid Tumor")
    fda_in = fda[fda_filt]
    fda_out = fda[~fda_filt]

    # Start working through muations
    if not report_dict["Sample Sequence Mutation Analysis"].empty:

        # add exons to muation table (only works on EGFR)
        mut_table = addExonData(report_dict, ExonDF)

        # Add a Type columns (will be used for excel output)
        mut_table.insert(0, "Type", "Mutation")

        # Special Case for Lung NSCLC Splice site events on Exon 14
        filt = mut_table.Gene == "MET"

        # If there's a MET splice site event
        if not mut_table.loc[filt].empty:
            # check if it's exon 14
            met_ex_14_sk = METexon14skipping(report_dict)
            met_ex_14_sk.insert(0, "Type", "Mutation")

            # if so and ttype is Lung NSCLC, add to significant genes and findings
            if ttype in list_o_lungs and not met_ex_14_sk.empty:
                forexcel = met_ex_14_sk.copy()
                forexcel[
                    "Amino Acid Change"
                ] = "MET Exon 14 Splice Site Event"
                fortiered = met_ex_14_sk.copy()

                output_dict["significant_genes"] += ["MET"]

                met_filt = (
                    (fda["ttype"] == ttype)
                    & (fda["mut_type"] == "Splice Site Event")
                    & (fda["gene"] == "MET")
                )

                fda_comment, fda_code = fda.loc[
                    met_filt, ["Interpretation", "tiered_var_id"]
                ].values[0]
                fda_code = str(fda_code)

                output_dict["findings"] += [fda_comment]

                output_excel["significant_genes"] = output_excel[
                    "significant_genes"
                ].append(forexcel, ignore_index=True)

                len_rows, len_cols = fortiered.shape
                fortiered.insert(len_cols, "tier", len_rows * ["1"])
                fortiered.insert(
                    len_cols + 1,
                    "tiered_var_id",
                    len_rows * [fda_code],
                )
                tiered_dict[
                    "Sample Sequence Mutation Analysis"
                ] = tiered_dict[
                    "Sample Sequence Mutation Analysis"
                ].append(
                    fortiered, ignore_index=True
                )

                #  and remove those rows from further consideration
                mut_table = pd.concat(
                    [mut_table, met_ex_14_sk]
                ).drop_duplicates(keep=False)

            # if so and other ttype, just add to tier2 genes
            elif ttype not in list_o_lungs and not met_ex_14_sk.empty:
                forexcel = met_ex_14_sk.copy()
                forexcel[
                    "Amino Acid Change"
                ] = "MET Exon 14 Splice Site Event"
                output_dict["tier2_genes"] += ["MET"]
                output_excel["tier2_genes"] = output_excel[
                    "tier2_genes"
                ].append(forexcel, ignore_index=True)
                fortiered = met_ex_14_sk.copy()
                len_rows, len_cols = fortiered.shape
                fortiered.insert(len_cols, "tier", len_rows * ["2"])
                fortiered.insert(
                    len_cols + 1, "tiered_var_id", len_rows * ["0"]
                )

                tiered_dict[
                    "Sample Sequence Mutation Analysis"
                ] = tiered_dict[
                    "Sample Sequence Mutation Analysis"
                ].append(
                    fortiered, ignore_index=True
                )
                #         #  and remove those rows from further consideration
                mut_table = pd.concat(
                    [mut_table, met_ex_14_sk]
                ).drop_duplicates(keep=False)
    # #############################################################################

    #     # Check for main results table aginst fda validated and fda validated for other
    #  ttypes

    for correct_ttype, fda_table, key in zip(
        [True, False],
        [fda_in, fda_out],
        ["significant_genes", "tier2_genes"],
    ):
        if not report_dict["Sample Sequence Mutation Analysis"].empty:

            for _, row in mut_table.iterrows():
                #     # Approach here is to mask the fda table in multiple ways
                #     # 1) does the regular expression match amino acid change?
                re_mask = fda_table.re.apply(
                    lambda x: bool(
                        re.search(
                            str(x), str(row["Amino Acid Change"])
                        )
                    )
                )

                #     # 2) if there's an exon requirement, does it match?
                exon_mask = fda_table.exon.apply(
                    lambda x: True
                    if (pd.isnull(x) or (x == float(row.Exon)))
                    else False
                )

                #     # 3) Gene name must match
                #     # and 4) must be 'Mutation' (as opposed to Translocation, etc)
                gen_mut_mask = (fda_table.gene == row.Gene) & (
                    fda_table.mut_type == "Mutation"
                )

                #     # 5) if there's a Consequence requirement, does it match?
                conseq_mask = fda_table.consequence.apply(
                    lambda x: True
                    if (pd.isnull(x) or (x == row.Consequence))
                    else False
                )

                #   # 6) if there's a codon K3326 (currently 10/26/21 the only one) requirement, does it match?
                #  NOTE: If additional codon requirements are added, we will need to add additional functions
                #  to check

                codon_mask = fda_table.codon.apply(
                    lambda x: True
                    if (pd.isnull(x) or check_K3326_and_up(row))
                    else False
                )

                #   Gather all masks
                mask = (
                    re_mask
                    & exon_mask
                    & gen_mut_mask
                    & conseq_mask
                    & codon_mask
                )

                #     # If there are no matching fda_table approved mutations add the gene
                #     # to the 'other_genes' list
                if fda_table[mask].empty:
                    if not correct_ttype:
                        output_dict["other_genes"] += [row.Gene]

                        output_excel["other_genes"] = output_excel[
                            "other_genes"
                        ].append(row, ignore_index=True)

                        tiered_row = row.append(
                            pd.Series(
                                ["3", "0"],
                                index=["tier", "tiered_var_id"],
                            )
                        )
                        tiered_dict[
                            "Sample Sequence Mutation Analysis"
                        ] = tiered_dict[
                            "Sample Sequence Mutation Analysis"
                        ].append(
                            tiered_row, ignore_index=True
                        )

                # if we get a match, add it to the appropriate genes list and if correct ttype,
                # and add gene and description to the 'findings' list
                else:
                    output_dict[key] += [row.Gene]
                    output_excel[key] = output_excel[key].append(
                        row, ignore_index=True
                    )

                    if correct_ttype:
                        if row["Amino Acid Change"] == "n/a":
                            for _, fda_row in fda_table.loc[
                                mask
                            ].iterrows():
                                output_dict["findings"] += [
                                    row.Gene
                                    + " "
                                    + row.Consequence
                                    + ": "
                                    + fda_row.Interpretation
                                ]
                        else:
                            mask3 = (
                                fda_table.description.apply(
                                    lambda x: x[0:3]
                                )
                                == row["Amino Acid Change"][0:3]
                            )
                            if not fda_table[mask3].empty:
                                for _, fda_row in fda_table.loc[
                                    mask
                                ].iterrows():
                                    output_dict["findings"] += [
                                        row.Gene
                                        + " "
                                        + row["Amino Acid Change"]
                                        + ": "
                                        + fda_row.Interpretation
                                    ]
                            else:
                                for _, fda_row in fda_table.loc[
                                    mask
                                ].iterrows():
                                    output_dict["findings"] += [
                                        fda_row.gene
                                        + " "
                                        + fda_row.description
                                        + ": "
                                        + fda_row.Interpretation
                                    ]
                        fda_comment, fda_code = fda_table.loc[
                            mask, ["Interpretation", "tiered_var_id"]
                        ].values[0]
                        fda_code = str(fda_code)
                        tiered_row = row.append(
                            pd.Series(
                                ["1", fda_code],
                                index=["tier", "tiered_var_id"],
                            )
                        )
                    else:
                        tiered_row = row.append(
                            pd.Series(
                                ["2", "0"],
                                index=["tier", "tiered_var_id"],
                            )
                        )

                    tiered_dict[
                        "Sample Sequence Mutation Analysis"
                    ] = tiered_dict[
                        "Sample Sequence Mutation Analysis"
                    ].append(
                        tiered_row, ignore_index=True
                    )

        # Trans and fusions
        if not report_dict["Sample Translocation Analysis"].empty:

            mask = fda_table.mut_type == "Translocation"
            fda_subset = fda_table[mask]

            if (
                "Type"
                not in report_dict[
                    "Sample Translocation Analysis"
                ].columns
            ):
                report_dict["Sample Translocation Analysis"].insert(
                    0, "Type", "Translocation"
                )

            subset = report_dict["Sample Translocation Analysis"]

            # for translocation, either gene or partner just needs to be in fda list
            mask2 = (subset["Gene"].isin(fda_subset.gene)) | (
                subset["Partner"].isin(fda_subset.gene)
            )

            to_add = subset[mask2]

            if not to_add.empty:
                add_list = (
                    to_add[["Gene", "Partner"]]
                    .apply(
                        lambda x: x[1] + "::" + x[0] + " " + "Fusion",
                        axis=1,
                    )
                    .tolist()
                )
                output_dict[key] += add_list
                output_excel[key] = output_excel[key].append(
                    to_add, ignore_index=True
                )
                fortiered = to_add.copy()
                len_rows, len_cols = fortiered.shape

                if correct_ttype:
                    for _, row in to_add.iterrows():
                        for _, fda_row in fda_subset.iterrows():
                            if (row["Gene"] == fda_row["gene"]) or (
                                row["Partner"] == fda_row["gene"]
                            ):
                                output_dict["findings"].append(
                                    f"{row['Partner']}::{row['Gene']} Fusion: {fda_row['Interpretation']}"
                                )

                    fortiered.insert(
                        len_cols, "tier", len_rows * ["1"]
                    )
                    fortiered.insert(
                        len_cols + 1,
                        "tiered_var_id",
                        len_rows * ["0"],
                    )
                    for _, row in fortiered.iterrows():
                        for _, fda_row in fda_subset.iterrows():
                            if (row["Gene"] == fda_row["gene"]) or (
                                row["Partner"] == fda_row["gene"]
                            ):
                                row["tiered_var_id"] = fda_row[
                                    "tiered_var_id"
                                ]
                    tiered_dict[
                        "Sample Translocation Analysis"
                    ] = tiered_dict[
                        "Sample Translocation Analysis"
                    ].append(
                        fortiered, ignore_index=True
                    )

                if key == "tier2_genes":
                    fortiered.insert(
                        len_cols, "tier", len_rows * ["2"]
                    )
                    fortiered.insert(
                        len_cols + 1,
                        "tiered_var_id",
                        len_rows * ["0"],
                    )
                    tiered_dict[
                        "Sample Translocation Analysis"
                    ] = tiered_dict[
                        "Sample Translocation Analysis"
                    ].append(
                        fortiered, ignore_index=True
                    )

            elif not subset[~mask2].empty:
                output_dict["other_genes"] += (
                    subset[~mask2][["Gene", "Partner"]]
                    .apply(
                        lambda x: x[1] + "::" + x[0] + " " + "Fusion",
                        axis=1,
                    )
                    .tolist()
                )
                output_excel["other_genes"] = output_excel[
                    "other_genes"
                ].append(subset[~mask2], ignore_index=True)
                fortiered = subset[~mask2].copy()
                len_rows, len_cols = fortiered.shape
                fortiered.insert(len_cols, "tier", len_rows * ["3"])
                fortiered.insert(
                    len_cols + 1, "tiered_var_id", len_rows * ["0"]
                )
                tiered_dict[
                    "Sample Translocation Analysis"
                ] = tiered_dict[
                    "Sample Translocation Analysis"
                ].append(
                    fortiered, ignore_index=True
                )

        if not report_dict["Sample Amplification Analysis"].empty:

            mask = fda_table.mut_type == "Amplification"
            fda_subset = fda_table[mask]

            if (
                "Type"
                not in report_dict[
                    "Sample Amplification Analysis"
                ].columns
            ):
                report_dict["Sample Amplification Analysis"].insert(
                    0, "Type", "Amplification"
                )

            subset = report_dict["Sample Amplification Analysis"]

            # for amplification, we need to see Status == "Copy Number Event Detected" or "amplified"
            subset = subset[
                subset["Status"].isin(
                    ["Copy Number Event Detected", "Amplified"]
                )
            ]

            mask2 = subset["Gene"].isin(fda_subset.gene)
            to_add = subset[mask2]
            fortiered = tiered_dict[
                "Sample Amplification Analysis"
            ].copy()

            if not to_add.empty:
                output_dict[key] += (
                    to_add["Gene"]
                    .map(lambda x: x + " " + "Amplification")
                    .tolist()
                )
                output_excel[key] = output_excel[key].append(
                    to_add, ignore_index=True
                )
            if correct_ttype:
                mask_amp = fda_subset.gene.isin(subset["Gene"])
                fda_subdf = fda_subset[mask_amp][
                    ["gene", "Interpretation"]
                ].apply(
                    lambda x: x["gene"]
                    + " Amplification: "
                    + x["Interpretation"],
                    axis=1,
                )
                if len(fda_subdf) != 0:
                    output_dict["findings"] += fda_subdf.tolist()
                for _, row in to_add.iterrows():
                    for _, fda_row in fda_subset.iterrows():
                        if row["Gene"] == fda_row["gene"]:
                            fortiered = fortiered.append(
                                row.append(
                                    pd.Series(
                                        {
                                            "tier": "1",
                                            "tiered_var_id": fda_row[
                                                "tiered_var_id"
                                            ],
                                        }
                                    )
                                ),
                                ignore_index=True,
                            )
            if key == "tier2_genes":
                for _, row in to_add.iterrows():
                    for _, fda_row in fda_subset.iterrows():
                        if row["Gene"] == fda_row["gene"]:
                            fortiered = fortiered.append(
                                row.append(
                                    pd.Series(
                                        {
                                            "tier": "2",
                                            "tiered_var_id": "0",
                                        }
                                    )
                                ),
                                ignore_index=True,
                            )
            else:
                output_dict["other_genes"] += subset["Gene"][
                    ~mask2
                ].tolist()
                output_excel["other_genes"] = output_excel[
                    "other_genes"
                ].append(subset[~mask2], ignore_index=True)
                to_add = subset[~mask2].copy()
                len_rows, len_cols = to_add.shape
                to_add.insert(len_cols, "tier", len_rows * ["3"])
                to_add.insert(
                    len_cols + 1, "tiered_var_id", len_rows * ["0"]
                )
                fortiered = fortiered.append(
                    to_add, ignore_index=True
                )

            tiered_dict[
                "Sample Amplification Analysis"
            ] = tiered_dict["Sample Amplification Analysis"].append(
                fortiered, ignore_index=True
            )

    # De-duplicate, remove items in higher tier lits from lower tiers
    # and sort lists
    output_dict = {k: set(v) for k, v in output_dict.items()}

    output_dict["other_genes"] = (
        output_dict["other_genes"]
        - output_dict["tier2_genes"]
        - output_dict["significant_genes"]
    )

    output_dict["tier2_genes"] = (
        output_dict["tier2_genes"] - output_dict["significant_genes"]
    )

    output_dict = {k: sorted(list(v)) for k, v in output_dict.items()}

    # Sort output for excel

    output_excel = (
        output_excel["significant_genes"]
        .append(output_excel["tier2_genes"], ignore_index=True)
        .append(output_excel["other_genes"], ignore_index=True)
    )

    output_excel["Tier"] = output_excel["Tier"].ffill()

    if "Type" in output_excel.columns:
        output_excel.dropna(subset=["Type"], inplace=True)

    output_excel.sort_index(axis=1, inplace=True)

    # allow for custom sort order
    output_excel["Type"] = pd.Categorical(
        output_excel["Type"],
        [
            "MSI",
            "TMB-H",
            "Mutation",
            "Translocation",
            "Amplification",
            " ",
            "n/a",
        ],
    )

    # remove lower tier duplicates
    important_columns = output_excel.columns.tolist()
    important_columns.remove("Tier")

    output_excel = output_excel.sort_values(
        ["Tier", "Type", "Gene"]
    ).drop_duplicates(important_columns, keep="first")
    # remove lower tier dupes from tiered dfs
    for tbl in [
        "Sample Genomic Signatures",
        "Sample Sequence Mutation Analysis",
        "Sample Amplification Analysis",
        "Sample Translocation Analysis",
    ]:
        important_columns = tiered_dict[tbl].columns.tolist()
        important_columns.remove("tier")
        important_columns.remove("tiered_var_id")
        cat_tier_order = CategoricalDtype(
            ["1", "2", "3", "0"], ordered=True
        )
        tiered_dict[tbl]["tier"] = tiered_dict[tbl]["tier"].astype(
            cat_tier_order
        )
        tiered_dict[tbl] = (
            tiered_dict[tbl]
            .sort_values("tier")
            .drop_duplicates(important_columns, keep="first")
        )
        for col in ["type", "exon", "Type", "Exon"]:
            if col in tiered_dict[tbl].columns:
                tiered_dict[tbl] = tiered_dict[tbl].drop(
                    columns=[col]
                )

    # add tmb and msi to significant genes
    if msi == "MSI":
        output_dict["significant_genes"].append("MSI-High")
    if tmb == "H":
        output_dict["significant_genes"].append("TMB-High")

    return (
        output_dict,
        add_link_columns(output_excel, batch_dir, filename),
        output_excel,
        tiered_dict,
    )


if __name__ == "__main__":
    report_dict = parse_csv(
        "/mnt/c/Users/adzab/OneDrive - Mass General Brigham/Operations/PGDx/pgdx-reporting/ElioConnect_Output/B22-4_ETCR-CFX-545448/Reports/CID21-13393_NA21-5481_B22-4.CCR_ETC-RUO.csv"
    )
    fda_file = pd.read_excel(
        "/mnt/c/Users/adzab/OneDrive - Mass General Brigham/Operations/PGDx/pgdx-web-app/requirements/fda_comment.xlsx",
        engine="openpyxl",
    )
    batch_dir = "/mnt/c/Users/adzab/OneDrive - Mass General Brigham/Operations/PGDx/pgdx-reporting/ElioConnect_Output/B22-4_ETCR-CFX-545448"
    ExonDF = createExonTable(
        "/mnt/c/Users/adzab/OneDrive - Mass General Brigham/Operations/PGDx/pgdx-web-app/requirements/refGene-CDS_HG19"
    )
    filename = "testing.html"
    interped = interp(
        report_dict, fda_file, batch_dir, ExonDF, filename
    )
    for tbl in [
        "Sample Genomic Signatures",
        "Sample Sequence Mutation Analysis",
        "Sample Amplification Analysis",
        "Sample Translocation Analysis",
    ]:
        print(interped[3][tbl])
