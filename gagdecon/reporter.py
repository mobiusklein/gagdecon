import os
import urllib
import csv
import brainpy
from collections import OrderedDict

try:
    from cStringIO import StringIO
except:
    from io import StringIO

try:
    rich_output = True
    import jinja2
    from matplotlib import rcParams as mpl_params
    from matplotlib import pyplot as plt
    from matplotlib.axes import Axes
    from ms_peak_picker.utils import draw_peaklist, draw_raw

    mpl_params.update({
        'figure.facecolor': 'white',
        'figure.edgecolor': 'white',
        'figure.figsize': (5, 3.5),
        'font.size': 10,
        # 72 dpi matches SVG
        # this only affects PNG export, as SVG has no dpi setting
        'savefig.dpi': 72,
        # 10pt still needs a little more room on the xlabel:
        'figure.subplot.bottom': .125})

    def png_plot(figure, **kwargs):
        buffer = render_plot(figure, format='png', **kwargs)
        return "<img src='data:image/png;base64,%s'>" % urllib.quote(buffer.getvalue().encode("base64"))

    def svg_plot(figure, **kwargs):
        buffer = render_plot(figure, format='svg', **kwargs)
        return buffer.getvalue()

    def render_plot(figure, **kwargs):
        if isinstance(figure, Axes):
            figure = figure.get_figure()
        if "height" in kwargs:
            figure.set_figheight(kwargs["height"])
        if "width" in kwargs:
            figure.set_figwidth(kwargs['width'])
        if kwargs.get("bbox_inches") != 'tight' or kwargs.get("patchless"):
            figure.patch.set_visible(False)
            figure.axes[0].patch.set_visible(False)
        buffer = StringIO()
        figure.savefig(buffer, **kwargs)
        plt.close(figure)
        return buffer

    def envelope_to_vector(envelope):
        mzs = []
        intensities = []
        for mz, intensity in envelope:
            mzs.append(mz - .000001)
            intensities.append(0.)
            mzs.append(mz)
            intensities.append(intensity)
            mzs.append(mz + .000001)
            intensities.append(0.)
        return (mzs), (intensities)

    def draw_envelope(env, ax=None):
        ax = draw_raw(*envelope_to_vector(env), ax=ax)
        lo, hi = ax.get_xlim()
        lo -= 0.5
        hi += 0.5
        ax.set_xlim(lo, hi)
        return ax

    def draw_tid(composition, charge, ax=None):
        tid = brainpy.isotopic_variants(composition, charge=charge)
        if ax is None:
            fig, ax = plt.subplots(1)
        ax = draw_peaklist(tid, ax=ax)
        lo, hi = ax.get_xlim()
        lo -= 0.5
        hi += 0.5
        ax.set_xlim(lo, hi)
        return ax

    class RichOutputWriter(object):
        def __init__(self, results, losses=None, path='report.html'):
            if losses is None:
                losses = []
            losses = [loss for loss in losses if loss is not None]
            self.results = results
            self.losses = losses
            self.path = path

        def write(self):
            loader = jinja2.FileSystemLoader(os.path.dirname(__file__))
            jenv = jinja2.Environment(loader=loader)
            jenv.filters['png_plot'] = png_plot
            jenv.globals.update(draw_tid=draw_tid, draw_envelope=draw_envelope)

            template = jenv.get_template("template.jinja2")
            g = template.stream(results=self.results, losses=self.losses)
            f = open(self.path, 'wb')
            i = 0
            for chunk in g:
                i += 1
                f.write(chunk)
                if i % 100:
                    f.flush()
            f.close()
            return self.path

    class FullSpectrumPlot(object):
        def __init__(self, results, peaklist, path='spectrum.pdf'):
            self.results = results
            self.peaklist = peaklist
            self.path = path
            self.ax = None

        def _draw_whole_spectrum(self):
            ax = draw_peaklist(self.peaklist, lw=0.25, alpha=0.5)
            self.ax = ax
            return ax

        def _highlight_matched_envelopes(self):
            ax = self.ax
            for result in self.results:
                score = result['score']
                a = max(1 - score, 0.1)
                draw_peaklist(result.fit.experimental, ax=ax, color='red', alpha=a, lw=0.5)

        def _ghost_in_theoretical_envelopes(self):
            ax = self.ax
            for result in self.results:
                draw_peaklist(result.fit.theoretical, ax=ax, color='green', alpha=0.5, lw=0.5)

        def write(self):

            ax = self._draw_whole_spectrum()
            self._highlight_matched_envelopes()
            self._ghost_in_theoretical_envelopes()

            max_yticks = 500
            xloc = plt.MaxNLocator(max_yticks)
            ax.xaxis.set_major_locator(xloc)

            fig = ax.get_figure()
            fig.set_figwidth(120 * 5)
            fig.savefig(self.path, bbox_inches="tight")
            plt.close(fig)

except ImportError:
    rich_output = False


class CSVOutputWriter(object):
    def __init__(self, results, path, fieldmap=None):
        self.results = results
        self.path = path
        self.fieldmap = OrderedDict()
        self.fieldmap['mz'] = 'm/z'
        self.fieldmap["charge"] = 'charge'
        self.fieldmap["score"] = "score"
        self.fieldmap['chain_length'] = "chain length"
        self.fieldmap['intensity'] = 'intensity'
        self.fieldmap["ppm_error"] = 'ppm error'
        self.fieldmap.update(fieldmap)

    def write(self):
        writer = csv.DictWriter(open(self.path, 'wb'), fieldnames=self.fieldmap.values())
        writer.writeheader()
        results = self.results
        for obj in results:
            writer.writerow({v: obj[k] for k, v in self.fieldmap.items()})
