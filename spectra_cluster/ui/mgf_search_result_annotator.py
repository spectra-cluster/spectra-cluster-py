"""
The mgf_search_result_annotator embeds identification data in MGF files to be
processed by the spectra_cluster algorithm tool suite.

This tool adds search results to an MGF file by adding the
identified peptide sequence as the SEQ= field to MGF files. This identification
data is picked up by the spectra-cluster tools and added to the .clustering output
files. Modification data is currently omitted (this is not the case in the
internal PRIDE Cluster pipeline).

The spectra-cluster tools (ie. the spectra-cluster-cli tool) expect identification
data to be embedded in the processed MGF files. Even though this method is
unorthodox it significantly simplifies the development of clustering tools as these
do not have to worry about the used search engine or search result formats.
Additionally, when building the PRIDE Cluster resource we have to rely on this
technique since the identification and spectrum data is exported from the PRIDE
Archive database.

Usage:
    mgf_search_result_annotator.py --input=<spectra.mgf> --search=<search_result.mzid> --output=<annotated_spectra.mgf>
                                   [--format=<MSGF+>] [--fdr=<0.01>] [--decoy_string=<REVERSED>]
    mgf_search_result_annotator.py (--help | --version)

Options:
    -i, --input=<spectra.mgf>             The original MGF file to use as input.
    -s, --search=<search_result.mzid>     The path to the search result. **Note**: The search
                                          must have been  performed on the input MGF file directly.
                                          Otherwise, the matching between identification data and
                                          spectra may go wrong.
    -o, --output=<annotated_spectra.mgf>  Path to where the annotated MGF file should be written to.
    -f, --format <MSGF+>                  The format of the search results. Possible options
                                          are "MSGF+", "MSAmanda", "Scaffold", "XTandem".
                                          [default: "MSGF+"]
    -d, --fdr=<0.01>                      Define the FDR by which the input search results are
                                          filtered. If the FDR is set to '2' for Scaffold output,
                                          the original cut-off is used. [default: 0.01]
    --decoy_string=<REVERSED>             The string to use to identify decoy proteins.
                                          [default: REVERSED]
    -h, --help                            Print this help message.
    -v, --version                         Print the current version.
"""

from docopt import docopt
import sys
import csv
import operator
import os
from pyteomics import mzid


csv.field_size_limit(sys.maxsize)


def parse_msgfplus(filename, fdr):
    """
    Extracts the PSMs from a MSGF+ search result file.

    :param filename: Filename of the MSGF+ search result file (only text file supported)
    :param fdr: Target FDR as fractional (ie. 0.01 for 1%)
    :return: A list of PSM objects
    """
    msgfplus_results = list()
    with open(filename, newline="") as result_file:
        msgfplus_result_reader = csv.DictReader(result_file, delimiter="\t", quoting=csv.QUOTE_NONE)

        for msgfplus_psm in msgfplus_result_reader:
            # ignore all PSMs below the set FDR
            if float(msgfplus_psm['PepQValue']) > fdr:
                continue

            psm = Psm(int(msgfplus_psm["SpecID"][6:]), msgfplus_psm["Peptide"])
            msgfplus_results.append(psm)

    return msgfplus_results


def parse_msamanda(filename, fdr, mgf_filename=None):
    """
    Extracts identification data from a MSAmanda search result file.

    :param filename: Filename of the MSAmanda search output (text file format)
    :param fdr: The target FDR as fractional (ie. 0.01 for 1%)
    :param mgf_filename: If set, only PSMs that were retrieved from the set search
                         result file will be taken into consideration
    :return: A list of PSM objects
    """
    msamanda_results = list()
    short_mgf_filename = mgf_filename

    # remove any path information from the filename
    if "/" in short_mgf_filename:
        short_mgf_filename = short_mgf_filename[short_mgf_filename.rfind("/") + 1:]

    with open(filename, newline="") as result_file:
        msamanda_result_reader = csv.DictReader(filter(lambda row: row[0] != '#', result_file), delimiter="\t", )

        for msamanda_psm in msamanda_result_reader:
            # ignore incorrectly formatted lines
            if msamanda_psm['Rank'] is None:
                continue

            # ignore all non-rank 1 PSMs
            if int(msamanda_psm['Rank']) > 1:
                continue

            # ignore all PSMs that are coming from a different file
            if msamanda_psm['Filename'] != short_mgf_filename:
                continue

            # label decoy hits
            if "REVERSED_" in msamanda_psm["Protein Accessions"]:
                msamanda_psm["Decoy"] = True
            else:
                msamanda_psm["Decoy"] = False

            # fix the MsAmanda Score
            if "," in msamanda_psm['Amanda Score']:
                msamanda_psm['Amanda Score'] = float(msamanda_psm['Amanda Score'].replace(",", "."))
            else:
                msamanda_psm['Amanda Score'] = float(msamanda_psm['Amanda Score'])

            # simply save the result to enable FDR filtering
            msamanda_results.append(msamanda_psm)

        # sort according to score
        msamanda_results.sort(key=operator.itemgetter('Amanda Score'), reverse=True)

        # filter decoys
        filtered_psms = list()
        n_target = 0
        n_decoy = 0

        for msamanda_psm in msamanda_results:
            if msamanda_psm["Decoy"]:
                n_decoy += 1
            else:
                n_target += 1

            current_fdr = n_decoy * 2 / (n_target + n_decoy)

            if current_fdr <= fdr:
                psm_object = Psm(Psm.MISSING_INDEX, msamanda_psm["Sequence"], msamanda_psm['Title'])
                filtered_psms.append(psm_object)
            else:
                break

        return filtered_psms


def parser_mzident(filename, score_field, title_field=None,
                   fdr=0.01, larger_score_is_better=False, decoy_string="DECOY"):
    """
    A general parsing function for mzIdentML files.

    Several exporters of mzIdentML do not report the correct spectrum indexes. X!Tandem, for example,
    uses the spectrum's title as "id" instead of the correct "index=N" format for MGF files. Therefore,
    it is possible to supply the index_field and title_field separately. Later, missing indexes will be
    resolved through the titles.

    :param filename: The path to the mzIdentML file
    :param score_field: The name of the score's field (**Important**: do not supply the accession
                        but only the name)
    :param title_field: The name of the field supplying the spectrum's title (in SpectrumIdentificationResult).
    :param fdr: Target FDR (default 0.01). If set to "2" the original cut-off is used.
    :param larger_score_is_better: Logical indicating whether better scores mean a more reliable
                                   result. Default is False as most search engines report
                                   probabilities
    :param decoy_string: String used to identify decoy proteins.
    :return: A list of PSM objects
    """
    mzid_psms = list()

    # load all PSMs from the file
    with mzid.read(filename) as object_reader:
        with mzid.read(filename) as reader:
            for spec_ref in reader:
                for spec_ident in spec_ref["SpectrumIdentificationItem"]:
                    # filter based on original FDR if set right away
                    if fdr == 2 and not spec_ident["passThreshold"]:
                        continue

                    # only use rank 1 ids
                    if spec_ident["rank"] > 1:
                        continue

                    if score_field not in spec_ident:
                        raise Exception("Failed to find supplied score field '" + score_field +
                                        "' in mzIdentML file.")
                    if title_field is not None and title_field not in spec_ref:
                        raise Exception("Failed to find supplied title field '" + title_field +
                                        "' in mzIdentML file.")

                    mzid_psm = dict()

                    mzid_psm["score"] = spec_ident[score_field]

                    # the index should be used as id
                    if spec_ref["spectrumID"][:6] == "index=":
                        mzid_psm["index"] = int(spec_ref["spectrumID"][6:])
                    else:
                        mzid_psm["index"] = Psm.MISSING_INDEX

                    # spectrum title is optional in mzIdentML
                    if title_field is not None:
                        mzid_psm["title"] = spec_ref[title_field]
                    elif "spectrum title" in spec_ref:
                        mzid_psm["title"] = spec_ref["spectrum title"]

                    # get the sequence in an mzIdentML "secure" way
                    peptide_evidence_id = spec_ident["PeptideEvidenceRef"][0]["peptideEvidence_ref"]
                    peptide_evidence = object_reader.get_by_id(peptide_evidence_id)
                    peptide_id = peptide_evidence["peptide_ref"]
                    mzid_psm["sequence"] = object_reader.get_by_id(peptide_id)["PeptideSequence"]
                    # TODO: PTMs are stored in peptide["Modification"]

                    # get the protein details to detect decoy hits
                    protein = object_reader.get_by_id(peptide_evidence["dBSequence_ref"])
                    protein_accession = protein["accession"]
                    # add the name if present just to be sure to catch any mention of the decoy status
                    if "name" in protein:
                        protein_accession += protein["name"]
                    if "protein description" in protein:
                        protein_accession += protein["protein description"]

                    mzid_psm["is_decoy"] = decoy_string in protein_accession

                    mzid_psms.append(mzid_psm)

    # sort the psms based on probability
    mzid_psms.sort(key=operator.itemgetter('score'), reverse=larger_score_is_better)

    # filter decoys
    filtered_psms = list()
    n_target = 0
    n_decoy = 0

    for mzid_psm in mzid_psms:
        # only filter if the FDR wasn't set to 2
        if fdr != 2:
            if mzid_psm["is_decoy"]:
                n_decoy += 1
            else:
                n_target += 1

            current_fdr = n_decoy * 2 / (n_target + n_decoy)

            if current_fdr > fdr:
                break

        # convert the psm
        if not mzid_psm["is_decoy"]:
            filtered_psms.append(Psm(mzid_psm["index"], mzid_psm["sequence"], mzid_psm["title"]))

    return filtered_psms


def parse_xtandem_mzident(filename, fdr, decoy_string="REVERSED"):
    """
    Parses X!Tandem output in the mzIdentML format.

    :param filename: The path to the mzIdentML file.
    :param fdr: The target fdr (ie. 0.01 for 1%)
    :param decoy_string: The decoy string to use to identify decoy proteins.
    :return: A list of PSM objects.
    """

    psms = parser_mzident(filename, score_field="X\\!Tandem:expect",
                          title_field="spectrumID", decoy_string=decoy_string,
                          fdr=fdr)

    # fix the incorrect titles
    for psm in psms:
        if "RTINSECONDS=" in psm.get_title():
            offset = psm.title.find("RTINSECONDS=")
            psm.title = psm.title[:offset].strip()

    return psms


def parse_scaffold(filename, fdr, decoy_string="REVERSED"):
    """
    Parses identification data from a Scaffold mzid output file

    :param filename: The filename of the mzIdentML file
    :param fdr: Target FDR as fractional (ie. 0.01 for 1%). If set to '2'
                the original cut-off is used.
    :param decoy_string: The string to identify decoy proteins.
    :return: A list of PSM objects
    """
    return parser_mzident(filename=filename,
                          score_field="Scaffold:Peptide Probability",
                          fdr=fdr,
                          decoy_string=decoy_string)


class Psm:
    """
    Represent a peptide spectrum match extracted from the search results.

    :ivar index: The 0-based index of the spectrum
    :ivar sequence: The peptide sequence
    :ivar title: The spectrum's "title"
    """
    MISSING_INDEX = -1

    def __init__(self, index, sequence, title=None):
        # 0-based index of the spectrum
        self.index = index
        self.sequence = sequence
        self.title = title

    def get_index(self):
        return self.index

    def set_index(self, index):
        self.index = index

    def get_sequence(self):
        return self.sequence

    def get_title(self):
        return self.title

    def set_title(self, title):
        self.title = title


def create_title_to_index_dict(mgf_filename):
    """
    Creates a dictionary with the spectrum's title as key and the 0-based index as value.

    :param mgf_filename: The MGF file to parse.
    :return: A dict with the spectra' titles as keys and their 0-based index as value.
    """
    title_to_index = dict()
    with open(mgf_filename, "r") as mgf_file:
        current_spec_index = 0
        for line in mgf_file:
            # ignore non-title lines
            if line.strip()[0:6] != "TITLE=":
                continue

            title = line[6:].strip()
            title_to_index[title] = current_spec_index
            current_spec_index += 1

    return title_to_index


def fix_missing_index(search_results, title_to_index, echo=False):
    """
    Some search engines do not report the spectrum's index
    in the search results but only the titles. This function uses
    the title_to_index dictionary to fix these missing index.

    :param search_results: A list of PSM objects.
    :param title_to_index: A dict with the spectra's titles as key and their index as value.
    :param echo: Boolean indicating whether a summary should be printed to stdout (Default false)
    :return: A float representing the fraction of matched missing indexes
    """
    fixed_indexes = 0
    has_missing_index = 0

    for psm in search_results:
        if psm.get_index() == Psm.MISSING_INDEX:
            has_missing_index += 1

            if psm.get_title() in title_to_index:
                psm.set_index(title_to_index[psm.get_title()])
                fixed_indexes += 1
            else:
                print("  Failed to find index for PSM with title '" + psm.get_title() + "'")

    if echo:
        print("Matched " + str(fixed_indexes) + "/" + str(has_missing_index) + " missing indices")

    if has_missing_index > 0:
        return fixed_indexes / has_missing_index
    else:
        # if no missing indexes are present, simply return 1 to
        # indicate that everything is matched
        return 1


def write_annotated_mgf(input_mgf, sequence_dictionary, output_mgf):
    """
    Creates an annotated MGF file based on the input MGF file and the sequence
    dictionary.

    :param input_mgf: The source MGF file.
    :param sequence_dictionary: A dict with the 0-based spectrum index as key and
                                the peptide sequence as value.
    :param output_mgf: The path to write the newly created MGF file to.
    """
    with open(output_mgf, "w") as output_file:
        with open(input_mgf, "r") as mgf_file:

            current_spec_index = 0

            for line in mgf_file:
                # simply copy the content
                output_file.write(line)
                # add the sequence information after the TITLE=
                if (line[0:6] == "TITLE=") & (current_spec_index in sequence_dictionary.keys()):
                    output_file.write("SEQ=" + sequence_dictionary[current_spec_index] + "\n")

                if line[0:8] == "END IONS":
                    current_spec_index += 1


def main():
    arguments = docopt(__doc__, version='mgf_search_result_annotator.py 1.0 BETA')

    input_file = arguments["--input"]
    search_file = arguments["--search"]
    output_file = arguments["--output"]

    # make sure input and search file exist
    if not os.path.isfile(input_file):
        print("Error: Cannot find peak list file '" + input_file + "'")
        sys.exit(1)
    if not os.path.isfile(search_file):
        print("Error: Cannot find search result file '" + search_file)
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(output_file):
        print("Error: Output file '" + output_file + "' already exists")
        sys.exit(1)

    # Create the title => index map
    title_to_index = create_title_to_index_dict(input_file)

    # Process the search result file
    search_results = list()
    search_format = arguments["--format"]

    fdr = float(arguments["--fdr"])

    if search_format.lower() == "msgf++":
        search_results = parse_msgfplus(search_file, fdr)
    elif search_format.lower() == "msamanda":
        search_results = parse_msamanda(search_file, fdr, input_file)
    elif search_format.lower() == "scaffold":
        search_results = parse_scaffold(filename=search_file, fdr=fdr,
                                        decoy_string=arguments["--decoy_string"])
    elif search_format.lower() == "xtandem":
        search_results = parse_xtandem_mzident(filename=search_file, fdr=fdr,
                                               decoy_string=arguments["--decoy_string"])
    else:
        print("Error: Unknown search engine result format set. "
              "Allowed values are: [MSGF+, MSAmanda, Scaffold, XTandem]")

    print("Extracted " + str(len(search_results)) + " PSMs @ FDR = " + str(fdr))

    # fix missing indices
    rel_matched = fix_missing_index(search_results, title_to_index, True)

    if rel_matched != 1:
        print("Error: Failed to match missing indexes.")
        sys.exit(1)

    # convert the search results into a dictionary
    sequence_dictionary = dict()

    for psm in search_results:
        # if there's a duplicate, use the lexical lower version
        if psm.get_index() in sequence_dictionary.keys():
            if sequence_dictionary[psm.get_index()] < psm.get_sequence():
                continue

        # use the spectrum's index as the key
        sequence_dictionary[psm.get_index()] = psm.get_sequence()

    # Create the annotated MGF file
    write_annotated_mgf(input_file, sequence_dictionary, output_file)

    print("Annotated MGF file successfully written to " + output_file)

if __name__ == "__main__":
    main()
