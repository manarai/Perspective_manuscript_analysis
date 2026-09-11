"""Build Figure 2 for the Perspective revision.

Panels:
  a  Twelve computations placed on (Q, omega). y-axis is symlog with
     linthresh = 0.02, so NF-kB (+0.009), Neubert-Caswell (+0.012),
     TCR (+0.075) and the heterogeneous repressilator (+0.095) sit
     visibly clear of the zero boundary. Horizontal bars = range of
     omega across the three per-state rescalings in
     results/sensitivity.csv. All twelve pairs of symmetric /
     heterogeneous variants (repressilator, toggle, relay) are joined
     by a connector. Region labels sit in the correct corners and the
     background is not shaded — the two boundary lines (omega = 0
     theorem boundary, Q = 1 stated convention) carry the classification.
     Single ink for markers so the reader does not read region
     membership out of colour.
  b  Six amplifying models paired with their spectrum-matched normal
     surrogates (nu = 0 by construction). y-axis baseline is 1.00
     exactly, so bar heights encode excess gain honestly. Each bar is
     annotated with G_max (3 decimals when excess < 0.1) and with the
     Henrici index nu of A, so the reader can see the TCR / Neubert-
     Caswell ν-similarity with widely different G_max.
  c  ||exp(At)||_2 vs normalised time for one representative of each
     of the four regions, annotated with G_max. tau = -1/alpha where
     alpha = max_i Re(lambda_i(A)) is the spectral abscissa.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
from models import all_models


ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)

INK = "#1e1e1e"
INK_LIGHT = "#666666"
BAR_A = "#c0392b"
BAR_NULL = "#b0b0b0"

DISPLAY = {
    "tcr":               "TCR",
    "nfkb":              "NF-κB",
    "murphy_miller":     "Murphy–Miller",
    "predator_prey":     "Neubert–Caswell",
    "repressilator":     "Repressilator (sym)",
    "repressilator_het": "Repressilator (het)",
    "goodwin":           "Goodwin",
    "selkov":            "Selkov",
    "toggle":            "Toggle (sym)",
    "toggle_het":        "Toggle (het)",
    "relay":             "Relay (sym)",
    "relay_het":         "Relay (het)",
}

PAIRS = [
    ("repressilator", "repressilator_het"),
    ("toggle", "toggle_het"),
    ("relay",  "relay_het"),
]

# omega symlog parameters
LINTHRESH = 0.02
LINSCALE  = 1.2


def load_descriptors():
    rows = {}
    with (RESULTS / "descriptors.csv").open() as f:
        for r in csv.DictReader(f):
            key = (r["model"], r["variant"])
            for k in ("alpha", "Q", "omega", "G_max", "nu", "r_PR"):
                r[k] = float(r[k])
            rows[key] = r
    return rows


def load_sensitivity():
    """Read the distribution over 500 random per-state metric probes.

    Returns per-model, per-descriptor: baseline, p05, p95, and sign_flip
    (fraction of draws whose sign flips relative to baseline, or None
    when the descriptor's sign is not meaningful).
    """
    ranges = {}
    with (RESULTS / "sensitivity_dist.csv").open() as f:
        for r in csv.DictReader(f):
            sf = r["frac_sign_flip"]
            ranges.setdefault(r["model"], {})[r["descriptor"]] = {
                "baseline": float(r["baseline"]),
                "p05": float(r["p05"]),
                "p95": float(r["p95"]),
                "sign_flip": float(sf) if sf != "" else None,
            }
    return ranges


def panel_a(ax, desc, sens):
    xlim = (-0.4, 4.9)
    ylim = (-5.0, 5.0)

    ax.set_xlim(*xlim)
    ax.set_yscale("symlog", linthresh=LINTHRESH, linscale=LINSCALE)
    ax.set_ylim(*ylim)

    ax.axhline(0.0, color=INK, lw=0.9, zorder=1)
    ax.axvline(1.0, color=INK_LIGHT, lw=0.8, ls="--", zorder=1)

    # Region labels in the correct corners.
    region_labels = [
        (xlim[0] + 0.10, 4.6,   "Amplifying,\nnon-oscillatory", "left",  "top"),
        (xlim[1] - 0.10, 4.6,   "Amplifying &\noscillating",    "right", "top"),
        (xlim[1] - 0.10, -4.6,  "Oscillating,\nnot amplifying", "right", "bottom"),
        (xlim[0] + 0.10, -4.6,  "Neither",                      "left",  "bottom"),
    ]
    for x, y, txt, ha, va in region_labels:
        ax.text(x, y, txt, ha=ha, va=va, fontsize=8.8, color=INK,
                weight="semibold", zorder=2, linespacing=1.15)

    lookup = {k: desc[(k, "A")] for k in DISPLAY}

    # ω sign-flippers under 500 random per-state metric probes — flagged
    # amber when the sign of ω flips in >=50 % of draws (i.e., the
    # assignment is decided by the metric as often as by the operator).
    SIGN_FLIP_COLOR = "#d97706"
    SIGN_FLIP_THRESHOLD = 0.50
    sign_flip = set()
    for key in lookup:
        sf = sens[key]["omega"]["sign_flip"]
        if sf is not None and sf >= SIGN_FLIP_THRESHOLD:
            sign_flip.add(key)

    # Label positions. All Q = 0 models pile up on the y-axis; symlog
    # spreads the near-zero values so we can label at their actual y
    # with a leader.
    LABEL_POS = {
        # Q = 0 stack.
        "tcr":               ( 0.14,  0.28,  "left"),
        "predator_prey":     ( 0.14,  0.0055,"left"),
        "relay_het":         ( 0.14, -0.10,  "left"),
        "relay":             ( 0.14, -0.35,  "left"),
        "toggle_het":        ( 0.14, -0.95,  "left"),
        "toggle":            ( 0.14, -2.20,  "left"),
        # Interior points — right cluster spread to avoid overlap with the
        # symmetric-repressilator range bar at Q ≈ 4.0.
        "murphy_miller":     ( 3.05,  3.70,  "left"),
        "selkov":            ( 3.05,  0.35,  "left"),
        "goodwin":           ( 2.25, -0.85,  "left"),
        "nfkb":              ( 3.40,  0.030, "right"),
        "repressilator":     ( 3.05, -1.15,  "left"),
        "repressilator_het": ( 4.70,  0.28,  "right"),
    }

    # Point + range: use error-bar caps so the marker is unambiguously
    # the baseline value and the bar reads as its 5-95 % percentile range
    # over 500 random per-state metric probes.
    for key, r in lookup.items():
        s = sens[key]["omega"]
        omin, omax = s["p05"], s["p95"]
        Q, w = r["Q"], r["omega"]
        omin_c = max(omin, ylim[0])
        omax_c = min(omax, ylim[1])
        is_flip = key in sign_flip
        ecolor  = SIGN_FLIP_COLOR if is_flip else INK
        elw     = 1.6 if is_flip else 1.2
        capsize = 4.0 if is_flip else 3.0

        ax.errorbar([Q], [w],
                    yerr=[[max(0.0, w - omin_c)], [max(0.0, omax_c - w)]],
                    fmt="none", ecolor=ecolor, elinewidth=elw,
                    capsize=capsize, capthick=elw, alpha=0.85, zorder=3)
        ax.scatter([Q], [w], s=48,
                   c=(SIGN_FLIP_COLOR if is_flip else INK),
                   edgecolor="white", linewidth=0.9, zorder=4)

        lx, ly, ha = LABEL_POS[key]
        leader_gap = 0.03 if ha == "left" else -0.03
        ax.plot([Q, lx - leader_gap], [w, ly],
                color=INK, lw=0.5, alpha=0.45, zorder=3)
        label_color = SIGN_FLIP_COLOR if is_flip else INK
        ax.text(lx, ly, DISPLAY[key], fontsize=7.8,
                color=label_color, ha=ha, va="center", zorder=5)

    # Sign-flip flag — placed in the empty upper-mid region of the panel
    # where there is no data and no region label, so it reads as a key
    # rather than as a point annotation.
    ax.text(1.20, 2.0,
            "amber  =  ω sign-flips in ≥ 50 % of 500 metric probes",
            ha="left", va="center",
            fontsize=8.0, color=SIGN_FLIP_COLOR, style="italic",
            zorder=6)

    # Custom y-ticks.
    yticks = [-3, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 3]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{v:g}" for v in yticks])
    ax.yaxis.set_minor_locator(mpl.ticker.SymmetricalLogLocator(
        base=10.0, linthresh=LINTHRESH, subs=np.arange(2, 10)))

    ax.set_xlabel(r"$Q\ =\ |\mathrm{Im}\,\lambda_1|\,/\,|\mathrm{Re}\,\lambda_1|$")
    ax.set_ylabel(r"$\omega\ =\ \lambda_{\max}((A + A^{\top})/2)$")
    ax.set_title("a  Descriptor plane over twelve computations",
                 loc="left", fontsize=10.5, weight="semibold")
    ax.tick_params(direction="out", length=3, which="major")
    ax.tick_params(direction="out", length=1.5, which="minor")
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def _fmt_gmax(g):
    excess = g - 1.0
    return f"{g:.3f}" if abs(excess) < 0.1 else f"{g:.2f}"


def panel_b(ax, desc):
    # Neubert-Caswell and TCR are placed adjacent so the reader can
    # read the nu-similarity / excess-gain-mismatch off the bar tops.
    amplifiers = ["predator_prey", "tcr",
                  "repressilator_het", "nfkb",
                  "selkov", "murphy_miller"]
    x = np.arange(len(amplifiers))
    width = 0.42
    g_A    = [desc[(k, "A")]["G_max"]      for k in amplifiers]
    nu_A   = [desc[(k, "A")]["nu"]         for k in amplifiers]
    excess_A = [g - 1.0 for g in g_A]

    # Amplifier excess-gain bars (from 0).
    bars_A = ax.bar(x, excess_A, width, color=BAR_A,
                    label=r"$A$: excess gain $G_{\max} - 1$")

    # Surrogate markers: short horizontal segments sitting at excess = 0
    # for each amplifier position, colour-coded distinctly.
    for xi in x:
        ax.hlines(0.0, xi - width/2, xi + width/2,
                  color="#4a4a4a", lw=3.2, zorder=5,
                  clip_on=False)
    surrogate_proxy = Line2D([0], [0], color="#4a4a4a", lw=3.2,
                             label=(r"$A_{\rm null}$ "
                                    r"(spectrum-matched, $\nu = 0$):"
                                    r"  excess = 0"))

    # Annotate G_max and nu above each red bar. clip_on=False so the
    # leftmost (Neubert-Caswell) label does not clip against the axes.
    ymax = max(excess_A)
    for i, (xi, e, g, nu) in enumerate(zip(x, excess_A, g_A, nu_A)):
        ha = "left" if i == 0 else "center"
        xoff = -0.15 if i == 0 else 0.0
        ax.text(xi + xoff, e + ymax * 0.02,
                f"$G_{{\\max}}$ = {_fmt_gmax(g)}\n$\\nu$ = {nu:.2f}",
                ha=ha, va="bottom", fontsize=7.6, color=BAR_A,
                linespacing=1.15, clip_on=False)

    ax.axhline(0, color=INK, lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY[k] for k in amplifiers], rotation=25,
                       ha="right", fontsize=8.5)
    ax.set_ylabel(r"Excess gain,  $G_{\max} - 1$")
    ax.set_ylim(0, ymax * 1.28)
    ax.legend(handles=[bars_A, surrogate_proxy],
              loc="upper left", frameon=False, fontsize=8.2)
    ax.set_title("b  Spectrum-matched surrogates remove transient gain",
                 loc="left", fontsize=10.5, weight="semibold")
    ax.tick_params(direction="out", length=3)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def panel_c(ax):
    ms = {m.key: m for m in all_models()}
    reps = [
        ("murphy_miller", "Murphy–Miller (amplifying & oscillating)", "#c0392b"),
        ("tcr",           "TCR (amplifying, non-oscillatory)",         "#e59866"),
        ("goodwin",       "Goodwin (oscillating, not amplifying)",     "#2874a6"),
        ("toggle",        "Toggle (neither)",                          "#616a6b"),
    ]

    ts_norm = np.linspace(0.0, 6.0, 500)

    handles = []
    for key, label, color in reps:
        A = ms[key].A
        alpha = float(np.max(np.linalg.eigvals(A).real))
        tau = -1.0 / alpha if alpha < 0 else 1.0
        gains = np.array([np.linalg.norm(expm(A * (t * tau)), 2)
                          for t in ts_norm])
        gmax = gains.max()
        zorder = 4 if key == "murphy_miller" else 3
        (line,) = ax.plot(ts_norm, gains, color=color, lw=1.8,
                          alpha=0.95, zorder=zorder,
                          label=f"{label}   $G_{{\\max}}$ = {gmax:.2f}")
        handles.append(line)

    ax.axhline(1.0, color=INK_LIGHT, lw=0.8, ls="--")
    ax.set_xlabel(r"$t\ /\ \tau$"
                  "     "
                  r"($\tau = -1/\alpha$;  $\alpha = \max_i \mathrm{Re}\,\lambda_i(A)$)")
    ax.set_ylabel(r"$\|e^{At}\|_2$")
    ax.set_title("c  Transient gain trajectories", loc="left",
                 fontsize=10.5, weight="semibold")
    ax.set_xlim(0, 6)
    ax.set_ylim(0.0, 3.4)
    ax.legend(handles=handles, loc="upper right", frameon=False,
              fontsize=8.2, handlelength=1.6, borderpad=0.3,
              labelspacing=0.35)
    ax.tick_params(direction="out", length=3)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def build():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "figure.dpi": 130,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    desc = load_descriptors()
    sens = load_sensitivity()

    fig = plt.figure(figsize=(11.0, 8.6))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.28, 1.0],
                  hspace=0.44, wspace=0.28,
                  left=0.08, right=0.98, top=0.94, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    panel_a(ax_a, desc, sens)
    panel_b(ax_b, desc)
    panel_c(ax_c)

    for ext in ("pdf", "png"):
        out = FIGDIR / f"fig2.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"wrote {out.relative_to(ROOT)}")
    plt.close(fig)


if __name__ == "__main__":
    build()
