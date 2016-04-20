import glypy
from glypy.structure.monosaccharide import ReducedEnd
from glypy.composition import glycan_composition
from glypy.composition.glycan_composition import MolecularComposition, SubstituentResidue
from collections import OrderedDict


hexnac = glypy.MonosaccharideResidue.from_iupac_lite("HexNAc")
hexose = glypy.MonosaccharideResidue.from_iupac_lite("Hex")
hexa = glypy.MonosaccharideResidue.from_iupac_lite("HexA")
hexn = glypy.MonosaccharideResidue.from_iupac_lite("HexN")

water = MolecularComposition("water", glypy.Composition("H2O"))
sulfate = SubstituentResidue.from_iupac_lite("sulfate")
acetyl = SubstituentResidue.from_iupac_lite("acetyl")

_anhydroman = glypy.monosaccharides.Man
_anhydro = glypy.Substituent("anhydro")
_anhydroman.add_substituent(_anhydro, position=2, parent_loss=glypy.Composition("H"))
anhydroman = glycan_composition.MonosaccharideResidue.from_monosaccharide(_anhydroman)

nhloss = MolecularComposition("nhloss", glypy.Composition("NH"))


class GlycosaminoglycanComposition(object):
    def __init__(self, composition, chain_length):
        self.glycan_composition = composition
        self.composition = composition.total_composition()
        self.chain_length = chain_length

    def __iter__(self):
        for key in self.composition:
            yield key, self.composition[key]

    def __getitem__(self, key):
        try:
            x = self.composition[key]
        except (TypeError, KeyError):
            x = 0
        if x != 0:
            return x
        else:
            return self.glycan_composition[key]

    def mass(self):
        return self.glycan_composition.mass()

    def __repr__(self):
        return "GlycosaminoglycanComposition(%r, %d)" % (self.glycan_composition, self.chain_length)

    def packed(self):
        pack = "[%(HexNAc)d,%(HexA)d,%(HexN)d,%(acetyl)d,%(sulfate)d]" % self.glycan_composition
        if self.glycan_composition[anhydroman]:
            pack += "+%dAnhydroMan" % self.glycan_composition[anhydroman]

        keys = set(self.glycan_composition.keys())
        leftovers = keys - {hexn, hexa, hexnac, acetyl, sulfate, anhydroman}
        parts = [
            "+%d%s" % (self.glycan_composition[part], part.name) if self.glycan_composition[part] > 0 else
            "-%d%s" % (self.glycan_composition[part], part.name)
            for part in leftovers]
        return ''.join(([pack] + parts))


class CombinatorialGenerator(object):
    def __init__(self, hexn, hexa, hexnac, sulfate=0, acetyl=0, anhydroman=0,
                 losses=None, chain_length=0, sulfation_ratio=(0.75, 1.5),
                 reduced=False):
        if losses is None:
            losses = [None]
        if anhydroman:
            reduced = True
        self.hexn = hexn
        self.hexa = hexa
        self.hexnac = hexnac
        self.sulfate = sulfate
        self.acetyl = acetyl
        self.anhydroman = anhydroman
        self.losses = losses
        self.chain_length = chain_length
        self.sulfation_ratio = sulfation_ratio
        self.reduced = reduced

    def generate(self):
        min_sulfate, max_sulfate = [int(s * self.chain_length) for s in self.sulfation_ratio]
        for n_acetyl in range(min(self.acetyl, self.hexn) + 1):
            for n_sulfate in range(min_sulfate, min(max_sulfate, self.sulfate + 1)):
                for loss in self.losses:
                    composition = glypy.GlycanComposition()
                    composition[hexn] = self.hexn
                    composition[hexa] = self.hexa
                    composition[hexnac] = self.hexnac
                    composition[sulfate] = n_sulfate
                    composition[acetyl] = n_acetyl
                    composition[anhydroman] = self.anhydroman
                    if self.reduced:
                        composition.reducing_end = ReducedEnd()
                    if loss is not None:
                        composition[loss] = -1
                    yield GlycosaminoglycanComposition(composition, self.chain_length)

    __iter__ = generate


class ChainCompositionBuilder(object):
    def __init__(self, gag_type, length_range, anhydroman=0, losses=None, composition_rules=None):
        if composition_rules is None:
            composition_rules = {}
        if losses is None:
            losses = [None]
        self.gag_type = gag_type
        self.length_range = length_range
        self.anhydroman = anhydroman
        self.losses = losses
        self.composition_rules = composition_rules

    def component_names(self):
        constants = OrderedDict({
            hexa: "HexA",
            hexn: "HexN",
            hexnac: "HexNAc",
            sulfate: "SO3",
            acetyl: "Ac"
        })
        if self.anhydroman:
            constants[anhydroman] = "AnhydroMan"

        for loss in self.losses:
            if loss is None:
                continue
            name = loss.name
            constants[loss] = name

        return constants

    def infer_components(self, chain_length):
        if self.gag_type == 'hs':
            hexnac = 0
            hexn = chain_length / 2 - \
                self.anhydroman if chain_length % 2 == 0 else (chain_length + 1) / 2 - self.anhydroman
        else:
            hexn = 0
            hexnac = chain_length / 2 - \
                self.anhydroman if chain_length % 2 == 0 else (chain_length + 1) / 2 - self.anhydroman
        hexa = chain_length / 2
        acetyl = min(3, hexn)
        sulfate = hexn * 3 + hexnac * 2 + hexa + 2 * self.anhydroman

        components = {
            "sulfate": sulfate, "hexa": hexa, "hexnac": hexnac,
            "hexn": hexn, "acetyl": acetyl, "anhydroman": self.anhydroman
        }
        components.update(self.composition_rules)
        return components

    def generate(self):
        for chain_length in range(*sorted(self.length_range)):
            components = self.infer_components(chain_length)
            combin = CombinatorialGenerator(chain_length=chain_length, losses=self.losses, **components)
            for gag in combin:
                yield gag

    __iter__ = generate
