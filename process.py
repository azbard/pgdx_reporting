# %%
import os

# %%
def progress(count, total, status=""):
    """
    A simple progress bar generator
    """
    bar_len = 9
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = f"{'=' * filled_len}{'-' * (bar_len - filled_len)}"

    # sys.stdout.write(f"[{bar}] {percents}% ... {status}\r")
    # sys.stdout.flush()
    return f"[{bar}] {percents}% ... {status}\r"


# %%
def pgdx_process(batch_dir, req_dir, mode="print"):

    input_directory = os.path.join(batch_dir, "Reports")
    igv_dir = os.path.join(batch_dir, "IGV")
    pathologist_dir = os.path.join(batch_dir, "All_Pathologist")
    wiki_dir = os.path.join(batch_dir, "Wiki")

    HappyHippo = (
        "\n"
        "    ________________________________________\n"
        "    / I'm getting to work!                   \\\n"
        "    \ Hope you have great a day :)           /\n"
        "    ----------------------------------------\n"
        "            \  ^__^\n"
        "            \  (oo)\_______\n"
        "               (__)\       )\/\\\n"
        "                    ||----w |\n"
        "                    ||     ||\n"
    )

    # print(HappyHippo)
    print(HappyHippo)

    total = 4
    i = 0

    i += 1
    if mode == "print":
        print(progress(i, total, status="Importing Dependencies"))
    # elif mode == "yield":
    #     yield progress(i, total, status="Importing Dependencies")
    import pandas as pd
    import time

    i += 1
    if mode == "print":
        print(progress(i, total, status="Importing Custom Code "))
    # elif mode == "yield":
    #     yield progress(i, total, status="Importing Custom Code ")
    from pgdx_reporting import reading

    i += 1
    if mode == "print":
        print(progress(i, total, status="Importing Custom Code "))
    # elif mode == "yield":
    #     yield progress(i, total, status="Importing Custom Code ")
    from pgdx_reporting import writing
    from pgdx_reporting import UploadToSql

    #
    # Initialize Directories
    # Plan is to run in the network root
    # then define the batch folder and add adresses from there.

    i += 1
    if mode == "print":
        print(progress(i, total, status="Setup Complete        \n"))
    # elif mode == "yield":
    #     yield progress(i, total, status="Setup Complete        \n")
    # if any required directories don't exist then create them
    new_dirs = [wiki_dir, pathologist_dir, igv_dir]

    for new_dir in new_dirs:
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)

    # Import Data Requirements

    # List of compatible versions [Assay Name, Assay Version, Platform Version ]
    compatible_versions = [
        ["PGDx elio tissue complete (03)", "2.0.0.4", "v1.4.1-10"],
        ["PGDx elio tissue complete - RUO (07)", "3.4.4.1", "v1.3.3-16"],
        ["PGDx elio tissue complete - RUO (07)", "4.0.0.4", "v1.4.1-10"],
        ["PGDx elio tissue complete - RUO (07)", "4.1.1.0", "v1.5.2-7"],
    ]

    # List of fda validated gene mutations/alterations
    fda = pd.read_excel(os.path.join(req_dir, "FDA_comment.xlsx"), engine="openpyxl")
    # Work with latest version
    fda = fda.sort_values("version").drop_duplicates(
        subset=[
            "ttype",
            "mut_type",
            "gene",
            "description",
            "exon",
            "re",
            "consequence",
            "codon",
        ],
        keep="last",
    )

    # List of Exon positions
    ExonDF = reading.createExonTable(os.path.join(req_dir, "refGene-CDS_HG19"))

    with open(os.path.join(req_dir, "html_template.html")) as f:
        html_template = f.read()

    # Initialize
    counter = 0
    success_counter = 0
    problem_files = str()

    # Create/open log file (auto-create if a new year)
    log_file = open(os.path.join(req_dir, "log" + time.strftime("%Y") + ".txt"), "a")

    # Long Batch Folder
    writing.log(f"Processing folder {os.path.basename(batch_dir)}", log_file)

    # count files to process:
    file_list = os.listdir(input_directory)
    n_files = len(file_list)

    # Loop through input directory
    for filename in file_list:

        counter += 1

        try:
            opening_msg = (
                f'\n{time.strftime("%Y%m%d-%H%M%S", time.localtime())}'
                f" - {str(counter)}/{str(n_files)}"
                f"\nOpening File: {filename}"
            )
            if mode == "print":
                print(writing.log(opening_msg, log_file,))
            # elif mode == "yield":
            #     yield writing.log(
            #         opening_msg, log_file,
            #     )

            fullpath = os.path.join(input_directory, filename)

            #     Read in data from the file
            report_dict = reading.parse_csv(fullpath)

            summary_df = report_dict["Case Summary"]
            filt = summary_df["Metric"] == "Overall Case Pass/Fail"

            if (summary_df.loc[filt, "Value"] == "PASS").all():

                interpreted_report_dict = reading.interp(
                    report_dict, fda, batch_dir, ExonDF, filename
                )
                #     Write the report to the report directory
                writing.write_report(
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
                )

                success_counter += 1

                progress_msg = (
                    f"\SQL Upload {filename} -"
                    f'{time.strftime("%Y%m%d-%H%M%S", time.localtime())}'
                )

                if mode == "print":
                    print(writing.log(progress_msg, log_file,))
                # elif mode == "yield":
                #     yield writing.log(
                #         progress_msg, log_file,
                #     )

                case = (
                    report_dict["Case Summary"]
                    .set_index("Metric")
                    .at["Case Name", "Value"]
                )

                tableDicOut = UploadToSql.tableCreation(
                    interpreted_report_dict[3],
                    fullpath,
                    os.path.join(
                        batch_dir, "TextFiles", case + ".translocation_ETC-RUO.csv"
                    ),
                    os.path.join(
                        batch_dir, "TextFiles", case + ".rawfoldchange_ETC-RUO.csv"
                    ),
                    os.path.join(
                        batch_dir, "TextFiles", case + ".samplehotspots_ETC-RUO.csv"
                    ),
                    os.path.join(req_dir, "FDA_comment.xlsx"),
                    mode,
                )

                # Write the output tables into created folder as .csv
                UploadToSql.writeOutputTables(fullpath, tableDicOut)
                # uploads the csv to sql
                UploadToSql.uploadToSql(
                    tableDicOut, os.path.join(req_dir, "FDA_comment.xlsx"), mode
                )
            else:
                writing.fail(
                    report_dict, wiki_dir, log_file, mode="print",
                )
                success_counter += 1

        except OSError as err:
            msg = "************ OS error: {0} ************".format(err)
            writing.log(msg, log_file)
            problem_files += filename + "\n" + msg + "\n\n"

    final_msg = "" if len(problem_files) == 0 else "Problem files:\n\n" + problem_files

    final_final_msg = (
        f"\nSuccessfully processed "
        f"{success_counter} out of {counter} files.\n"
        f"{final_msg}"
    )

    if mode == "print":
        print(writing.log(final_final_msg, log_file,))
    # elif mode == "yield":
    #     yield writing.log(
    #         final_final_msg, log_file,
    #     )

    log_file.close()
    complete_reported_file = os.path.join(batch_dir, "COMPLETE-REPORTED")                                                                                 ─╯
    with open(complete_reported_file, "w") as fp: 
        pass  


# %%
def pgdx_main(batch_dir, req_dir):
    complete_file = os.path.join(batch_dir, "COMPLETE")
    complete_reported_file = os.path.join(batch_dir, "COMPLETE-REPORTED")
    if os.path.exists(complete_file):
        print(f"hello {batch_dir}")
        if os.path.exists(complete_reported_file):
            print(f"{batch_dir} Already Complete")
        else:
            try:
                print(f"trying {batch_dir}")
                pgdx_process(batch_dir, req_dir)
                with open(complete_reported_file, "w") as fp:
                    pass
                print("great success")
            except:
                print("processing failed")
                # email me


if __name__ == "__main__":
    pgdx_process(
        "/mnt/pgdx_v1/ElioConnect_Output/B21-1146_ETCR-SPF-66845",
        "/mnt/pgdx_v1/.Bioinformatics/requirements",
    )
