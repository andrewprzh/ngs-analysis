#!/usr/bin/env python

import os
import sys
import cairo
import math

GENOME_START_POS = 1000
GENOME_END_POS = 100000
TOTAL_GENES = 10
SPACE_BETWEEN_LEVELS = 10
BASE_LEVEL = 3

CANVAS_WIDTH = min(2000, GENOME_END_POS - GENOME_START_POS)
CANVAS_HEIGHT = max(300, SPACE_BETWEEN_LEVELS * (TOTAL_GENES + BASE_LEVEL))

XMARGIN = 50
YMARGIN = 50

WIDTH = CANVAS_WIDTH + 2 * XMARGIN
HEIGHT = CANVAS_HEIGHT + 2 * YMARGIN

XSCALE = float(CANVAS_WIDTH) / float(GENOME_END_POS - GENOME_START_POS)
YSCALE = float(CANVAS_HEIGHT) / float(TOTAL_GENES + BASE_LEVEL)

CODON_SIZE = 3


def draw_element(contex, genome_coords, level):
    x1 = XMARGIN + genome_coords[0] * XSCALE
    x2 = XMARGIN + genome_coords[0] * XSCALE
    y = YMARGIN + level * YSCALE
    context.move_to(x1, y)
    context.line_to(x2, y)
    context.stroke()


def draw_gene(contex, coords):
    context.set_source_rgb(0.2, 0.2, 0.2)
    context.set_line_width(2)
    draw_element(contex, coords, 0)


def draw_transcript(contex, coords, transcript_num):
    context.set_source_rgb(0.5, 0.5, 0.5)
    context.set_line_width(3)
    draw_element(contex, coords, transcript_num + BASE_LEVEL)


def draw_exon(contex, coords, transcript_num):
    context.set_source_rgb(0.1, 0.1, 0.9)
    context.set_line_width(3)
    draw_element(contex, coords, transcript_num + BASE_LEVEL)


def draw_codon(contex, genome_coords, level):
    x = XMARGIN + (genome_coords[0] + genome_coords[1]) * XSCALE / 2
    y = YMARGIN + level * YSCALE
    contex.arc(x, y, CODON_SIZE, 0, 2 * math.pi)
    context.fill()
    context.stroke()


#transcript_num = 0 for reference drawing
def draw_start_codon(contex, coords, transcript_num = 0):
    context.set_source_rgb(0.1, 0.9, 0.1)
    draw_codon(contex, coords, transcript_num + BASE_LEVEL)


#transcript_num = 0 for reference drawing
def draw_stop_codon(contex, coords, transcript_num = 0):
    context.set_source_rgb(0.8, 0.1, 0.1)
    draw_codon(contex, coords, transcript_num + BASE_LEVEL)

def draw_coords(context):
    pass

def set_legend(context, gene_id, chromosome, coordinates):
    pass


def init(transcript_num, gene_coords):
    GENOME_START_POS = gene_coords[0]
    GENOME_END_POS = gene_coords[1]
    TOTAL_GENES = transcript_num

    CANVAS_WIDTH = min(5000, GENOME_END_POS - GENOME_START_POS)
    CANVAS_HEIGHT = SPACE_BETWEEN_LEVELS * (TOTAL_GENES + BASE_LEVEL)

    WIDTH = CANVAS_WIDTH + 2 * XMARGIN
    HEIGHT = CANVAS_HEIGHT + 2 * YMARGIN

    XSCALE = float(CANVAS_WIDTH) / float(GENOME_END_POS - GENOME_START_POS)
    YSCALE = float(CANVAS_HEIGHT) / float(TOTAL_GENES + BASE_LEVEL)


def init_context(fname):
    surface = cairo.SVGSurface(fname + ".svg", WIDTH, HEIGHT)
    context = cairo.Context(surface)
    context.rectangle(0, 0, WIDTH, HEIGHT)
    context.set_source_rgb(1, 1, 1)
    context.fill()
    context.stroke()
    return context




