import os
import argparse
import logging

from glypy import Composition
from gagdecon.builder import ChainCompositionBuilder, MolecularComposition
from gagdecon.deconvoluter import deconvolute, read_peaklist, pick_peaks_from_file
from gagdecon.reporter import CSVOutputWriter, rich_output

if rich_output:
    from gagdecon.reporter import RichOutputWriter

logger = logging.getLogger("gagdecon")


class GAGDeconvoluter(object):
    def __init__(self, gag_type, length_range, has_anhydromanose=False, losses=None, composition_rules=None,
                 pick_peaks=False):
        if composition_rules is None:
            composition_rules = {}
        self.gag_type = gag_type
        self.length_range = length_range
        self.has_anhydromanose = has_anhydromanose
        self.losses = losses
        self.composition_rules = composition_rules
        self.pick_peaks = pick_peaks

    def deconvolute_file(self, peaklist_path, charge_range=(-1, -10)):
        if self.pick_peaks:
            peaklist = pick_peaks_from_file(peaklist_path)
        else:
            peaklist = read_peaklist(peaklist_path)
        return self.deconvolute_peaklist(peaklist, charge_range)

    def deconvolute_peaklist(self, peaklist, charge_range=(-1, -10)):
        database_builder = ChainCompositionBuilder(
            self.gag_type, self.length_range, self.has_anhydromanose,
            self.losses, self.composition_rules)
        database = list(database_builder)
        return deconvolute(peaklist, database, charge_range)

    def _get_component_names(self):
        database_builder = ChainCompositionBuilder(
            self.gag_type, self.length_range, self.has_anhydromanose,
            self.losses, self.composition_rules)
        return database_builder.component_names()


def run(peaklist_path, gag_type, length_range, has_anhydromanose=False, losses=None,
        composition_rules=None, max_charge=-10, output_path=None, output_format=('csv',),
        pick_peaks=False):
    dec = GAGDeconvoluter(gag_type, length_range, has_anhydromanose, losses, composition_rules)

    logger.info("Deconvoluting and Matching %s", peaklist_path)
    results = dec.deconvolute_file(peaklist_path, charge_range=(-1, max_charge))

    if output_format == "html" and not rich_output:
        output_format = "csv"
        logger.error("HTML Output Not Supported: Missing Dependencies")

    if output_path is None:
        path = os.path.splitext(peaklist_path)[0]
        output_path = path + ".gagdecon."

    results = sorted(results, key=lambda x: x["score"])
    for output_format_ in output_format:
        output_path_ = output_path + output_format_
        if output_format_ == 'csv':
            logger.info("Writing CSV")
            fields = dec._get_component_names()
            CSVOutputWriter(results, output_path_, fields).write()
        elif output_format_ == 'html':
            logger.info("Writing HTML")
            RichOutputWriter(results, losses, output_path_).write()
        elif output_format_ == 'pickle':
            import pickle
            logger.info("Writing Pickle")
            pickle.dump(results, open(output_path_ + 'pkl', 'wb'))
    return results


app = argparse.ArgumentParser("gagdecon")
app.add_argument("-g", "--gag-type", choices=('hs', 'cs'), help='The type of GAG chain to produce')
# app.add_argument("-r", "--reduced", default=None, help="The reducing end composition, e.g. H2")
app.add_argument("-l", "--losses", action='append',
                 help="Any chemical loss compositions, e.g. H2O or NH. May specify more than once.")
app.add_argument("-a", "--has-anhydromanose", action='store_true', required=False, default=False,
                 help='Does the composition contain a reduced anhydromanose residue at the reducing end')
app.add_argument("-d", "--chain-length-range", nargs=2, type=int,
                 help="The range of chain lengths to consider, given as two numbers low and high")
app.add_argument("-c", "--max-charge", type=int, help='The maximum negative charge to search for')
app.add_argument("-o", "--output-path", required=False, default=None, help='The path to write all output to')
app.add_argument("-f", "--output-format", required=False, action='append', default=[],
                 choices=('csv', 'html', 'pickle'), help="The format to write output to. May specify more than once.")

app.add_argument("-p", "--pick-peaks", required=False, default=False, help="Should peak picking be performed")

app.add_argument("peaklist_path", help='Path to the peaklist file')


def main():
    logging.basicConfig(level='INFO', format="%(name)s: %(message)s")
    args = app.parse_args()
    molecular_composition_losses = [None]  # The None case for the loss-less case
    for loss in args.losses:
        logger.info("Converting loss %s -> %s", loss, Composition(loss))
        molecular_composition_losses.append(MolecularComposition(loss, Composition(loss)))
    length_range = sorted(map(int, args.chain_length_range))
    max_charge = -abs(args.max_charge)
    has_anhydromanose = bool(args.has_anhydromanose)
    gag_type = args.gag_type

    logger.info("GAG Chain Range: %d-%d" % tuple(length_range))

    # reducing_end_type = args.reduced
    # if reducing_end_type:
    #     reducing_end_type = Composition(reducing_end_type)

    output_path = args.output_path
    output_format = args.output_format
    if not output_format:
        output_format = ['csv']

    pick_peaks = args.pick_peaks

    run(
        args.peaklist_path, gag_type, length_range, has_anhydromanose, molecular_composition_losses,
        None, max_charge, output_path, output_format, pick_peaks=pick_peaks)


if __name__ == '__main__':
    main()
