import pandas as pd
import sqlalchemy
from sqlalchemy import create_engine, text
import os
import re
import sys
import warnings
import numpy
from psycopg2.extensions import register_adapter, AsIs

from pgdx_reporting.reading import interp

# These lines necessary for data translation
# between pandas and postgresql
def addapt_numpy_float64(numpy_float64):
    return AsIs(numpy_float64)


def addapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)


register_adapter(numpy.float64, addapt_numpy_float64)
register_adapter(numpy.int64, addapt_numpy_int64)

warnings.filterwarnings("ignore")

# %%
def caseInfo(filename):
    # Remove the junk at the end of the file name
    real_filename = os.path.basename(filename)
    file_stub = real_filename.split(".")
    id_area = file_stub[0]
    # Extract out the CID, NA, and Batch to distinct Variables
    id_list = ["-".join(id.split("-")[:2]) for id in id_area.split("_")]

    CID, NA, Batch = id_list
    return CID, NA, Batch


# %%
def makeColumnsLowercase(df):
    # Makes every column name lowercase
    df.columns = map(str.lower, df.columns)

    # Replaces all spaces with underscores in the column names
    df.columns = df.columns.str.replace(" ", "_")

    return df


def addIdentifiers(df, CID, NA, Batch):
    # Adds in the column and corresponding information for identifiers
    # starting with CID, NA, and then Batch
    df.insert(0, "cid_number", CID)
    df.insert(1, "na", NA)
    df.insert(2, "batch", Batch)

    return df


def smallTransform(df, filename):
    # Transpose to the table to make it 1 by X
    df_t = df.transpose()

    # Move the first column to the name and drop the first column
    df_t = df_t.rename(columns=df_t.iloc[0]).iloc[1:, :]

    # Remove the Index
    df_t.reset_index(drop=True, inplace=True)

    CID, NA, Batch = caseInfo(filename)

    df_t = addIdentifiers(df_t, CID, NA, Batch)

    df_t = df_t.iloc[:1]

    return df_t


def bigTransform(df, filename):
    # Transforms the first two columns
    dfLeft = smallTransform(df.iloc[:, 0:2], filename)

    # Converts the first column into a list
    firstColumnList = df.iloc[:, 0].tolist()

    # Creates and Index to join on
    dfLeft["fakeIndex"] = "2"

    # Checks the first two columns
    firstTwoColumnNames = df.iloc[:, :2].columns

    for column in df:
        # Checks  it is not the first two columns, to start transforming the rest of the data
        if column not in firstTwoColumnNames:

            # Combines the old column name with the first column list to create the new column names
            newColumnNames = [s + " " + column for s in firstColumnList]

            # Transposes the column into a row
            newDF = df[column].to_frame().transpose()

            # Renames the column names
            newDF.columns = newColumnNames

            # Creates same index on transitionary df
            newDF["fakeIndex"] = "2"

            # Adds new DF to our existing df
            dfLeft = dfLeft.merge(newDF, on="fakeIndex", how="outer")

    # Drops fake index
    dfLeft.drop(["fakeIndex"], axis=1, inplace=True)

    return dfLeft


def tableCreation(
    bigDic,
    reportPath,
    translocationPath,
    rawfoldchangePath,
    samplehotspotPath,
    pgdx_fda_comment_path,
    mode="print",
):

    # Takes the CID, NA, and Batch from the reportPath
    CID, NA, Batch = caseInfo(reportPath)

    # ---------------------- Sample Hotspot Table --------------------
    # As this table does not interact with the main report, this is generated first and seperately
    # Reads in the csv
    samplehotspotTable = pd.read_csv(samplehotspotPath)

    # Adds in identifers and transforms the columns to be lowercase
    samplehotspotTable = makeColumnsLowercase(
        addIdentifiers(samplehotspotTable, CID, NA, Batch)
    )

    # Renames start and end to have location prefix for clarity
    samplehotspotTable.rename(
        {"start": "location_start", "end": "location_end"}, axis="columns", inplace=True
    )

    # defines the columns that will be used as the primary keys
    qc_hotspots_primaryKeys = [
        "cid_number",
        "na",
        "batch",
        "chrom",
        "location_start",
        "location_end",
        "mutant",
        "exon",
        "amino_acid",
    ]

    # since the primary keys can be null, this concatenates them into one column using | as the deliminator
    samplehotspotTable["primary_key"] = samplehotspotTable[
        qc_hotspots_primaryKeys
    ].apply(lambda row: "|".join(row.values.astype(str)), axis=1)

    # ------------------------------------------------------------------

    # This is the list of all the tables in the report that contain information relating to the batch
    batchList = [
        "Case Summary",
        "Run Information",
        "Run Sequencer Version Information",
        "Run Quality Metrics",
        "Run Unassigned Indices",
        "Run Consumables",
        "External Control Information",
    ]

    # This is a list of all the tables in the report that contain sample information that can be joined together
    sampleList = [
        "Sample Information",
        "Sample Quality Metrics",
        "Sample Contamination Analysis",
        "Sample Genomic Signatures",
    ]

    # Creates a table of just the CID, NA, and Batch information for everything to be joined onto
    batchTable = pd.DataFrame({"cid_number": [CID], "na": [NA], "batch": [Batch]})

    # Creates a table of just the CID, NA, and Batch information for everything to be joined onto
    sampleTable = pd.DataFrame({"cid_number": [CID], "na": [NA], "batch": [Batch]})

    for key in bigDic:

        # Quickly renames a column if it encounters that specific table
        if key == "Sample Contamination Analysis":
            bigDic["Sample Contamination Analysis"].rename(
                {"Pass/Fail": "Sample Contamination Pass/Fail"},
                axis="columns",
                inplace=True,
            )

        # Checks to see if the it is a small or large column and joins it on the initalized dataframe above

        # ---------------- Sample Table --------------

        # First it starts with a simple dataframe with 2 columns
        if key in batchList and len(bigDic[key].columns) == 2:

            # Creates a temporary dataframe which is transposed
            smallMergeDF = smallTransform(bigDic[key], reportPath)
            # Takes that transformed table and joins it on the intialized dataframe
            batchTable = batchTable.merge(
                smallMergeDF,
                how="outer",
                left_on=["cid_number", "na", "batch"],
                right_on=["cid_number", "na", "batch"],
                suffixes=("", "_y"),
            )

        # Now it repeats the same concept with the larger tables
        # Probably a cleaner way to do this, maybe discuss with Adam later?
        if key in batchList and len(bigDic[key].columns) > 2:

            # Makes a temporary dataframe that is transposed and has all the additional columns
            # transposed with more intuitive namining scheme
            bigMergeDF = bigTransform(bigDic[key], reportPath)

            # Joins the transformed tables onto the initialized dataframe
            batchTable = batchTable.merge(
                bigMergeDF,
                how="outer",
                left_on=["cid_number", "na", "batch"],
                right_on=["cid_number", "na", "batch"],
                suffixes=("", "_y"),
            )

        # ---------------- Batch Table --------------

        if key in sampleList and len(bigDic[key].columns) == 2:

            # Creates a temporary dataframe which is transposed
            smallMergeDF = smallTransform(bigDic[key], reportPath)
            # Takes that transformed table and joins it on the intialized dataframe
            sampleTable = sampleTable.merge(
                smallMergeDF,
                how="outer",
                left_on=["cid_number", "na", "batch"],
                right_on=["cid_number", "na", "batch"],
                suffixes=("", "_y"),
            )

        # Now it repeats the same concept with the larger tables
        # Probably a cleaner way to do this, maybe discuss with Adam later?
        if key in sampleList and len(bigDic[key].columns) > 2:

            # Makes a temporary dataframe that is transposed and has all the additional columns
            # transposed with more intuitive namining scheme
            bigMergeDF = bigTransform(bigDic[key], reportPath)

            # Joins the transformed tables onto the initialized dataframe
            sampleTable = sampleTable.merge(
                bigMergeDF,
                how="outer",
                left_on=["cid_number", "na", "batch"],
                right_on=["cid_number", "na", "batch"],
                suffixes=("", "_y"),
            )

        # ---------------- Translocation Table --------------

        if key == "Sample Translocation Analysis":
            # Checks if the dataframe is empty
            if not bigDic[key].empty:

                # Creates a dataframe of the amplification table
                translocationTable = bigDic[key]

                # Identifies the supplemental table
                translocationETC = pd.read_csv(translocationPath)

                # Merges the amplification table onto the supplemental table
                # as that table contains informations on all the genes

                translocationTable_n = translocationTable.merge(
                    translocationETC,
                    how="left",
                    left_on=["Gene", "Partner"],
                    right_on=["Gene", "Partner Gene"],
                ).dropna()

                translocationTable_r = (
                    translocationTable.merge(
                        translocationETC,
                        how="left",
                        left_on=["Gene", "Partner"],
                        right_on=["Partner Gene", "Gene"],
                    )
                    .drop(columns=["Partner Gene"])
                    .rename({"Gene_x": "Gene", "Gene_y": "Partner Gene"}, axis=1)
                    .dropna()
                )

                translocationTable = pd.concat(
                    [translocationTable_n, translocationTable_r]
                ).drop_duplicates()

                # Drops redundant columns
                translocationTable.drop(["Partner", "Breakpoint"], axis=1, inplace=True)

                # Adds in identifiers to the table
                translocationTable = addIdentifiers(translocationTable, CID, NA, Batch)

                # Makes the table's columns lowercase
                translocationTable = makeColumnsLowercase(translocationTable)

            else:
                translocationTable = pd.DataFrame()

        # ---------------- Amplification Table --------------

        if key == "Sample Amplification Analysis":

            # puts the amplification table into a new dataframe
            ampAnalysis = bigDic[key]

            # reads in the supplemental rawfoldchange table
            cnvTable = pd.read_csv(rawfoldchangePath)

            # Merges on the two tables together
            try:
                cnvTable = cnvTable.merge(ampAnalysis, how="left", on="Gene")
            except:
                msg = "Amp table may be empty!"
                if mode == "print":
                    print(msg)
                # elif mode == "yield":
                #     yield msg
            # status is a defined word in sql so had to be renamed
            cnvTable.rename(
                {"Status": "alteration_status"}, axis="columns", inplace=True
            )

            # Adds in identifiers to the table and makes it lowercase
            cnvTable = makeColumnsLowercase(addIdentifiers(cnvTable, CID, NA, Batch))

        # ---------------- Mutation Table --------------

        if key == "Sample Sequence Mutation Analysis":
            # Puts the information into a new dataframe
            mutTable = bigDic[key]

            # Alters the dataframe to make the columns lowercase and add identifiers
            mutTable = makeColumnsLowercase(addIdentifiers(mutTable, CID, NA, Batch))

            # defines the columns that will be used as the primary keys
            qc_hotspots_primaryKeys = [
                "cid_number",
                "na",
                "batch",
                "gene",
                "genomic_location",
                "reference",
                "mutant",
            ]

            # since the primary keys can be null, this concatenates them into one column using | as the deliminator
            mutTable["primary_key"] = mutTable[qc_hotspots_primaryKeys].apply(
                lambda row: "|".join(row.values.astype(str)), axis=1
            )

    # This drops all the additional columns in batch table
    # They were joined so that if there was a duplicate it would have a specific suffix
    # This looks for that suffix and drops it
    batchTable.drop(
        batchTable.filter(regex="_y$").columns.tolist(), axis=1, inplace=True
    )

    # Makes every column name lowercase in the batch table
    # Replaces all spaces with underscores in the column names
    batchTable = makeColumnsLowercase(batchTable)

    # drops the sample specific identifiers from the batch table
    batchTable.drop(["cid_number", "na", "total"], axis=1, inplace=True)

    # Reads in the most up to date FDA comment version and adds it to the batch table
    pgdx_fda_commentDF = pd.read_excel(pgdx_fda_comment_path, engine="openpyxl")
    batchTable["fda_comment_version"] = pgdx_fda_commentDF["version"].max()

    # Repeats for sample table
    sampleTable.drop(
        sampleTable.filter(regex="_y$").columns.tolist(), axis=1, inplace=True
    )
    sampleTable = makeColumnsLowercase(sampleTable)
    # Creates a dictionary where the keys are the table names, and the values are the corresponding DF
    tableDic = {
        "pgdx_batch": batchTable,
        "pgdx_case_msi": sampleTable,
        "pgdx_translocations": translocationTable,
        "pgdx_amplifications": cnvTable,
        "pgdx_mutations": mutTable,
        "pgdx_qc_hotspots": samplehotspotTable,
    }
    return tableDic


# This function writes out the tables
def writeOutputTables(reportPath, tableDicOut, output_dir):

    # Generates the CID, NA, Batch as variables to use in naming
    CID, NA, Batch = caseInfo(reportPath)

    newFolderPath = output_dir

    # Checks to see if the folder to momentarily place all of the files
    # If it doesn't exist, it makes a folder called 'sqlUpload'
    if not os.path.exists(newFolderPath):
        os.makedirs(newFolderPath)

    # Generates the filename and path for the batch
    batchFilename = os.path.join(newFolderPath + Batch + "_batchTableUpload.csv")

    # loops through dictionary of output tables
    for key in tableDicOut:

        # Creates a path and filename for each sample and its output table
        sampleFilename = os.path.join(
            newFolderPath + CID + "_" + NA + "_" + Batch + "_" + key + "_Upload.csv"
        )

        # For batches it checks if it exists
        if key == "pgdx_batch" and not os.path.exists(batchFilename):

            # If the batch file does not exist, it will create one
            with open(batchFilename, "x") as batchFile:

                # Transforms the dataframe to a string
                contents = tableDicOut[key].to_csv(index=False)

                # Writes the string out to the created file
                batchFile.write(contents)

        # Checks to see if the dataframe is empty as translocation tables can be empty
        # Also checks to see it is not already created
        # Lastly checks that it is not the batch table
        if (
            not tableDicOut[key].empty
            and key != "pgdx_batch"
            and not os.path.exists(sampleFilename)
        ):

            # if those criteria are met, it will create the csv
            with open(sampleFilename, "x") as sampleFile:

                # Converts dataframe to a string
                contents = tableDicOut[key].to_csv(index=False)

                # Writes the string to created file
                sampleFile.write(contents)


# Function generates SQL connection
def connectToDatabase(mode="print"):
    """Connect to the PostgreSQL database server"""

    try:
        # connect to the PostgreSQL server
        msg = "Connecting to the PostgreSQL database..."
        if mode == "print":
            print(msg)
        # elif mode == "yield":
        #     yield msg
        engine = create_engine(
            "postgresql://sv739:cpath%40mgh@172.27.81.127:31644/cider"
        )

    except (Exception, sqlalchemy.exc.DatabaseError) as error:
        msg = error
        if mode == "print":
            print(msg)
        # elif mode == "yield":
        #     yield msg

        sys.exit(1)

    msg = "Connection successful"
    if mode == "print":
        print(msg)
    # elif mode == "yield":
    #     yield msg

    return engine


# Function loops through the csv's created and uploads them to the sql database
def uploadToSql(tableDicOut, pgdx_fda_comment_path, mode="print"):

    # Establishes the connections to the SQL database
    engine = connectToDatabase()

    for filetype in tableDicOut:

        if filetype == "pgdx_batch":
            checkFDA(
                tableDicOut["pgdx_batch"], pgdx_fda_comment_path, engine,
            )

        # if filetype == "pgdx_translocations":
        #     msg = tableDicOut[filetype]
        #     if mode == "print":
        #         print(msg)
        # elif mode == "yield":
        #     yield msg

        try:
            if not tableDicOut[filetype].empty:
                # Uses the connection previously created
                with engine.begin() as connection:

                    # Finds the corresponding file in the dictionary and appends it to the sql database
                    # table with the same name
                    tableDicOut[filetype].to_sql(
                        filetype,
                        con=connection,
                        if_exists="append",
                        index=False,
                        schema="pgdx",
                    )
                    msg = filetype + " successfully uploaded to sql"
                    if mode == "print":
                        print(msg)
                    # elif mode == "yield":
                    #     yield msg

            else:
                msg = filetype + " empty - not uploaded to sql"
                if mode == "print":
                    print(msg)
                # elif mode == "yield":
                #     yield msg

        except (Exception, sqlalchemy.exc.DBAPIError) as error:
            msg = filetype + " failed to upload to sql"
            if mode == "print":
                print(msg)
            # elif mode == "yield":
            #     yield msg

            if "duplicate key value violates unique constraint" in str(error):
                msg = "Already Uploaded!"
                if mode == "print":
                    print(msg)
                # elif mode == "yield":
                #     yield msg

            else:
                msg = error
                if mode == "print":
                    print(msg)
                # elif mode == "yield":
                #     yield msg


# Checks to see if the fda table is up to date in sql and updates if not
def checkFDA(df, pgdx_fda_comment_path, engine, mode="print"):

    # Retrieving data to see what the most recent version uploaded to the database is
    stmt = text("SELECT max(version) from pgdx.pgdx_fda_comment")
    with engine.connect() as conn2:
        fda_comment_table_version = conn2.execute(stmt)

    # Checks to see if the database's most recent version matches the most recent FDA table version
    if not fda_comment_table_version.first()[0] == df["fda_comment_version"][0]:

        # if not matches will transform
        pgdx_fda_commentDF = pd.read_excel(pgdx_fda_comment_path, engine="openpyxl")

        try:
            # Uses the connection previously created
            with engine.begin() as connection:

                # Finds the corresponding file in the dictionary and appends it to the sql database
                # table with the same name
                pgdx_fda_commentDF.to_sql(
                    "pgdx_fda_comment",
                    con=connection,
                    if_exists="replace",
                    index=False,
                    schema="pgdx",
                )

        except (Exception, sqlalchemy.exc.DBAPIError) as error:
            msg = "Upload of FDA comment table failed"
            if mode == "print":
                print(msg)
            # elif mode == "yield":
            #     yield msg

            msg = error
            if mode == "print":
                print(msg)
            # elif mode == "yield":
            #     yield msg

        msg = "FDA comment table uploaded successfully"
        if mode == "print":
            print(msg)
        # elif mode == "yield":
        #     yield msg

    else:
        msg = "FDA version up to date"
        if mode == "print":
            print(msg)
        # elif mode == "yield":
        #     yield msg


if __name__ == "__main__":

    from pgdx_reporting.reading import parse_csv, interp, createExonTable

    batch_dir = os.path.abspath(
        os.path.join(
            "..", "..", "pgdx-reporting", "ElioConnect_Output", "B22-4_ETCR-CFX-545448"
        )
    )
    req_dir = os.path.abspath(os.path.join("..", "requirements"))
    fda_comment_path = os.path.join(req_dir, "FDA_Comment.xlsx")
    fda = pd.read_excel(fda_comment_path, engine="openpyxl")
    ExonDF = createExonTable(os.path.join(req_dir, "refGene-CDS_HG19"))
    filename = "CID21-13393_NA21-5481_B22-4.CCR_ETC-RUO.csv"
    reportPath = os.path.join(batch_dir, "Reports", filename)
    csv_dir = os.path.join(os.path.basename(os.path.abspath(req_dir)), "sqlUpload"))

    tranlocationPath = os.path.join(
        batch_dir, "TextFiles", "CID21-13393_NA21-5481_B22-4.translocation_ETC-RUO.csv"
    )
    rawfoldchangePath = os.path.join(
        batch_dir, "TextFiles", "CID21-13393_NA21-5481_B22-4.rawfoldchange_ETC-RUO.csv"
    )
    samplehotspotPath = os.path.join(
        batch_dir, "TextFiles", "CID21-13393_NA21-5481_B22-4.samplehotspots_ETC-RUO.csv"
    )

    pre_bigDic = parse_csv(reportPath)
    bigDic = interp(pre_bigDic, fda, batch_dir, ExonDF, filename)[3]

    case = pre_bigDic["Case Summary"].set_index("Metric").at["Case Name", "Value"]

    # Creates the dictionary with all of the dataframes to ouput
    tableDicOut = tableCreation(
        bigDic,
        reportPath,
        tranlocationPath,
        rawfoldchangePath,
        samplehotspotPath,
        fda_comment_path,
    )

    # Write the output tables into created folder as .csv
    writeOutputTables(reportPath, tableDicOut, csv_dir)

    # uploads the csv to sql
    uploadToSql(tableDicOut, fda_comment_path)
