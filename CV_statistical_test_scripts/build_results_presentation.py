#!/usr/bin/env python3
"""Build the PowerPoint summary of the pillar-angle tests."""

from __future__ import annotations

from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

from project_paths import presentation_dir, project_root

ASSETS = Path(
    "/Users/jiacheng/.cursor/projects/Users-jiacheng-Desktop-AngleMeasurementsFromKira/assets"
)
NAVY = RGBColor(0x1F, 0x4E, 0x79)
INK = RGBColor(0x1A, 0x1A, 0x1A)
MUTED = RGBColor(0x4A, 0x4A, 0x4A)
RULE = RGBColor(0xD0, 0xD7, 0xE2)


def _set_run(run, text: str, size: int, bold: bool = False, color: RGBColor = INK, font: str = "Calibri") -> None:
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = color
    run.font.name = font


def add_title(slide, text: str) -> None:
    box = slide.shapes.add_textbox(Inches(0.45), Inches(0.22), Inches(12.4), Inches(0.55))
    frame = box.text_frame
    frame.word_wrap = True
    paragraph = frame.paragraphs[0]
    _set_run(paragraph.add_run(), text, 28, bold=True, color=NAVY)
    line = slide.shapes.add_shape(
        1,  # MSO_SHAPE.RECTANGLE
        Inches(0.45),
        Inches(0.84),
        Inches(12.4),
        Emu(18000),
    )
    line.fill.solid()
    line.fill.fore_color.rgb = NAVY
    line.line.fill.background()


def add_text(slide, left, top, width, height, blocks: list[tuple[str, int, bool]]) -> None:
    box = slide.shapes.add_textbox(left, top, width, height)
    frame = box.text_frame
    frame.word_wrap = True
    frame.auto_size = None
    frame.margin_left = Inches(0.04)
    frame.margin_right = Inches(0.04)
    for index, (text, size, bold) in enumerate(blocks):
        paragraph = frame.paragraphs[0] if index == 0 else frame.add_paragraph()
        paragraph.alignment = PP_ALIGN.LEFT
        paragraph.space_after = Pt(8)
        paragraph.level = 0
        color = MUTED if size <= 13 else INK
        _set_run(paragraph.add_run(), text, size, bold=bold, color=color)
    frame.auto_size = None


def fit_picture(slide, path: Path, left, top, max_width, max_height):
    with Image.open(path) as image:
        width_px, height_px = image.size
    aspect = width_px / height_px
    box_aspect = max_width / max_height
    if aspect > box_aspect:
        width = max_width
        height = int(max_width / aspect)
    else:
        height = max_height
        width = int(max_height * aspect)
    picture_left = left + (max_width - width) // 2
    picture_top = top + (max_height - height) // 2
    slide.shapes.add_picture(str(path), picture_left, picture_top, width, height)


def blank(prs: Presentation):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    background = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, prs.slide_height)
    background.fill.solid()
    background.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    background.line.fill.background()
    slide.shapes._spTree.remove(background._element)
    slide.shapes._spTree.insert(2, background._element)
    return slide


def build(image_dir: Path, output: Path) -> None:
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    prs.core_properties.title = "Pillar angle statistical tests"
    prs.core_properties.subject = "p-values for the pillar-angle figures"

    slide = blank(prs)
    add_title(slide, "Pillar angle statistical tests")
    add_text(
        slide,
        Inches(0.55),
        Inches(1.3),
        Inches(12.1),
        Inches(5.4),
        [
            ("Each slide is one figure, with the p-value from the test that matches it.", 22, False),
            ("Coefficient of variation (CV) uses the sample standard deviation divided by the mean.", 20, False),
            ("A before/after CV test resamples embryos and keeps each embryo’s two angles together. The statistic is the percent change, (CV after − CV before) / CV before. There are 10,000 replicates and the seed is 42.", 20, False),
            ("The wls versus tyrosinase test uses 100,000 replicates and seed 42. The volume-bin slope uses 10,000 replicates and seed 42.", 20, False),
            ("An interval that contains 0, or a p-value above 0.05, means that result is compatible with no change.", 20, False),
        ],
    )

    slide = blank(prs)
    add_title(slide, "wls versus tyrosinase: CV by region")
    fit_picture(slide, image_dir / "wls_pillarAngle.png", Inches(0.45), Inches(1.05), Inches(12.4), Inches(2.35))
    add_text(
        slide,
        Inches(0.55),
        Inches(3.55),
        Inches(12.1),
        Inches(3.6),
        [
            ("Question: is the pillar-angle CV larger in wls than in tyrosinase?", 20, True),
            ("Anterior: wls CV 0.043 versus 0.033. Difference +0.010. Bootstrap p = 0.2836.", 18, False),
            ("Posterior: wls CV 0.031 versus 0.019. Difference +0.012. Bootstrap p = 0.1493.", 18, False),
            ("Ventral: wls CV 0.035 versus 0.018. Difference +0.017. Bootstrap p = 0.1691.", 18, False),
            ("Holm adjustment across the three regions gives p = 0.4479 for every region. The higher wls CV in each region is compatible with equal CVs.", 18, False),
        ],
    )

    slide = blank(prs)
    add_title(slide, "CV versus lumen volume")
    fit_picture(
        slide,
        image_dir / "volume_cv_posterior.png",
        Inches(0.35),
        Inches(1.05),
        Inches(6.3),
        Inches(6.1),
    )
    add_text(
        slide,
        Inches(6.8),
        Inches(1.15),
        Inches(6.1),
        Inches(5.9),
        [
            ("Weighted slope of CV against the 90 pL bin center. The same bins are used for every aspect. The screenshot is the posterior panel.", 16, False),
            ("Anterior, 119 embryos: slope −0.0065 per 90 pL. 95% interval −0.0088 to −0.0035. p < 0.0002.", 18, True),
            ("Posterior, 102 embryos: slope −0.0081 per 90 pL. 95% interval −0.0106 to −0.0039. p < 0.0002.", 18, True),
            ("Ventral, 51 embryos: slope −0.0258 per 90 pL. 95% interval −0.0324 to −0.0137. p = 0.0002.", 18, True),
            ("Every 95% interval excludes 0. Volume and age move together, so the slope is the association along the time course.", 15, False),
        ],
    )

    slide = blank(prs)
    add_title(slide, "PBS and HAase, 0.5 h")
    fit_picture(slide, image_dir / "haase.png", Inches(0.4), Inches(1.05), Inches(5.3), Inches(6.05))
    add_text(
        slide,
        Inches(5.9),
        Inches(1.1),
        Inches(7.0),
        Inches(6.0),
        [
            ("Angle, PBS, n = 10: mean change +4.0°. Wilcoxon p = 0.08398. Paired t p = 0.1095.", 17, False),
            ("Angle, HAase, n = 9: mean change +10.8°. Wilcoxon p = 0.01953. Paired t p = 0.01178.", 17, False),
            ("CV, PBS: 0.062 to 0.060, percent change −4.0%. 95% interval −28% to +33%. p = 0.7168.", 17, False),
            ("CV, HAase: 0.057 to 0.014, percent change −75%. 95% interval −88% to −48%. p = 0.0006.", 17, True),
            ("HAase minus PBS, percent CV change: −71 percentage points. 95% interval −112 to −37. p = 0.0028.", 17, True),
            ("One HAase embryo on 5 March has no after image and is omitted. The printed figure percents use the rounded CV labels.", 14, False),
        ],
    )

    slide = blank(prs)
    add_title(slide, "DMSO and forskolin, 1 h")
    fit_picture(slide, image_dir / "forskolin.png", Inches(0.35), Inches(1.05), Inches(5.5), Inches(6.05))
    add_text(
        slide,
        Inches(6.05),
        Inches(1.1),
        Inches(6.85),
        Inches(6.0),
        [
            ("The bracket on the figure is the angle comparison. Forskolin minus DMSO is +8.1°. Equal-variance t, df = 18, p = 0.01416. Exact Mann–Whitney p = 0.01150.", 17, True),
            ("DMSO angle, n = 10: mean change +2.7°. Wilcoxon p = 0.2754. Paired t p = 0.1684.", 17, False),
            ("Forskolin angle, n = 10: mean change +10.8°. Wilcoxon p = 0.001953. Paired t p = 0.001509.", 17, False),
            ("CV, DMSO: −34%. 95% interval −68% to −10%. p = 0.0104.", 17, False),
            ("CV, forskolin: −44%. 95% interval −51% to −22%. p = 0.0018.", 17, False),
            ("Forskolin minus DMSO, percent CV change: −9.7 percentage points. 95% interval −33 to +31. p = 0.6292.", 17, True),
        ],
    )

    slide = blank(prs)
    add_title(slide, "Anterior ablation")
    fit_picture(slide, image_dir / "anterior_ablation.png", Inches(0.35), Inches(1.0), Inches(7.15), Inches(3.85))
    add_text(
        slide,
        Inches(0.5),
        Inches(4.95),
        Inches(12.3),
        Inches(2.35),
        [
            ("Plotted embryos, n = 10. Angle mean change −11.3°. Paired t p = 0.0003595. Wilcoxon p = 0.001953. Lumen volume paired t p = 0.0005289.", 17, False),
            ("CV of those 10 embryos: 0.076 to 0.093, percent change +22%. 95% interval −5% to +111%. p = 0.1434.", 18, True),
            ("The printed labels 0.070 to 0.092 use 13 embryos. That set changes by +31%. 95% interval +1% to +106%. p = 0.0402.", 17, False),
        ],
    )

    slide = blank(prs)
    add_title(slide, "Posterior ablation")
    fit_picture(slide, image_dir / "posterior_ablation.png", Inches(0.35), Inches(1.0), Inches(7.15), Inches(4.15))
    add_text(
        slide,
        Inches(7.6),
        Inches(1.15),
        Inches(5.35),
        Inches(5.9),
        [
            ("Figure set, n = 10. Angle mean change −13.1°. Paired t p < 0.0001. Wilcoxon p = 0.001953.", 16, False),
            ("All 11 Prism pairs: angle paired t p < 0.0001. Lumen volume paired t p = 0.0003105.", 16, False),
            ("CV of the 10: 0.092 to 0.112, percent change +21%. 95% interval −0.6% to +50%. p = 0.0524.", 17, True),
            ("All 11 Prism pairs: percent CV change +24%. 95% interval −0.8% to +59%. p = 0.0548.", 16, False),
            ("Rows 7–11 have no embryo file name in the source folder. Row 10 is left out of the n = 10 CV.", 14, False),
        ],
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    prs.save(output)


def main() -> None:
    image_dir = project_root() / "source_figures"
    output = presentation_dir() / "pillar_angle_statistical_tests.pptx"
    build(image_dir, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
