#!/usr/bin/env python

import os
import sys
import cairo
import math

class GeneDrawer:
    SPACE_BETWEEN_LEVELS = 10
    LOW_SPACE_BETWEEN_LEVELS = 8
    BASE_LEVEL = 3
    XMARGIN = 100
    YMARGIN = 50
    CODON_SIZE = 3
    LABEL_X_MARGIN = 10
    LABEL_Y_MARGIN = 5

    def __init__(self, transcript_num, contig_num, gene_coords, fname):
        self.GENOME_START_POS = gene_coords[0]
        self.GENOME_END_POS = gene_coords[1]
        self.TOTAL_GENES = transcript_num
        self.TOTAL_CONTIGS = contig_num

        self.CANVAS_WIDTH = min(2000, self.GENOME_END_POS - self.GENOME_START_POS)
        self.GENE_CANVAS_HEIGHT = max(150, self.SPACE_BETWEEN_LEVELS * (self.TOTAL_GENES + self.BASE_LEVEL))
        self.CONTIG_CANVAS_HEIGHT = max(150, self.LOW_SPACE_BETWEEN_LEVELS * (self.TOTAL_CONTIGS + self.BASE_LEVEL))
        self.CANVAS_HEIGHT = self.CONTIG_CANVAS_HEIGHT + self.GENE_CANVAS_HEIGHT
        self.GENE_LEVEL = self.GENE_CANVAS_HEIGHT

        self.WIDTH = self.CANVAS_WIDTH + 2 * self.XMARGIN
        self.HEIGHT = self.CANVAS_HEIGHT + 2 * self.YMARGIN

        self.XSCALE = float(self.CANVAS_WIDTH) / float(self.GENOME_END_POS - self.GENOME_START_POS)
        self.GENE_YSCALE = float(self.GENE_CANVAS_HEIGHT) / float(self.TOTAL_GENES + self.BASE_LEVEL)
        self.CONTIG_YSCALE = float(self.CONTIG_CANVAS_HEIGHT) / float(self.TOTAL_CONTIGS + self.BASE_LEVEL)

        self.surface = cairo.SVGSurface(fname + ".svg", self.WIDTH, self.HEIGHT)
        self.context = cairo.Context(self.surface)
        self.context.rectangle(0, 0, self.WIDTH, self.HEIGHT)
        self.context.set_source_rgb(1, 1, 1)
        self.context.fill()
        self.context.stroke()


    def draw_element(self, genome_coords, level):
        x1 = self.XMARGIN + (genome_coords[0] - self.GENOME_START_POS)  * self.XSCALE
        x2 = self.XMARGIN + (genome_coords[1] - self.GENOME_START_POS)  * self.XSCALE
        y = 0
        if level < 0:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.GENE_YSCALE
        else:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.CONTIG_YSCALE
        self.context.move_to(x1, y)
        self.context.line_to(x2, y)
        self.context.stroke()
        #print(x1,x2,y)


    def label_alignment(self, coords, level, label):
        self.context.set_source_rgb(0.1, 0.1, 0.1)
        x = self.XMARGIN + (genome_coords[1] - self.GENOME_START_POS)  * self.XSCALE + LABEL_X_MARGIN
        y = 0
        if level < 0:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.GENE_YSCALE
        else:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.CONTIG_YSCALE
        self.context.move_to(x, y)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.set_font_size(10)
        self.context.show_text(label)
        self.context.stroke()


    def label_genome(self, coords, gene_id, chr_id):
        self.context.set_source_rgb(0.1, 0.1, 0.1)
        x = self.XMARGIN + ((genome_coords[1] + genome_coords[0]) / 2 - self.GENOME_START_POS)  * self.XSCALE + LABEL_X_MARGIN
        y = - LABEL_Y_MARGIN
        self.context.move_to(x, y)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.set_font_size(12)
        self.context.show_text(chr_id + ": " + str(genome_coords[0]) + "-" + str(genome_coords[1]) + "," + )
        self.context.stroke()

    def draw_codon(self, genome_coords, level):
        x = self.XMARGIN + (genome_coords[0] + genome_coords[1] - 2 * self.GENOME_START_POS) * self.XSCALE / 2
        y = 0
        if level < 0:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.GENE_YSCALE
        else:
            y = self.YMARGIN + self.GENE_LEVEL + level * self.CONTIG_YSCALE

        self.context.arc(x, y, self.CODON_SIZE, 0, 2 * math.pi)
        self.context.fill()
        self.context.stroke()


    def draw_gene(self, coords):
        self.context.set_source_rgb(0.2, 0.2, 0.2)
        self.context.set_line_width(3)
        self.draw_element(coords, 0)


    def draw_transcript(self, coords, transcript_num):
        self.context.set_source_rgb(0.5, 0.5, 0.5)
        self.context.set_line_width(1)
        self.draw_element(coords, - transcript_num - self.BASE_LEVEL)


    def draw_exon(self, coords, transcript_num):
        self.context.set_source_rgb(0.1, 0.1, 0.9)
        self.context.set_line_width(3)
        self.draw_element(coords, - transcript_num - self.BASE_LEVEL)

    def draw_contig(self, coords, contig_num):
        self.context.set_source_rgb(0.5, 0.5, 0.5)
        self.context.set_line_width(1)
        self.draw_element(coords, contig_num + self.BASE_LEVEL)

    def draw_block(self, coords, contig_num, coverage = 1):
        self.context.set_source_rgb(0.1, 0.1, 0.9)
        self.context.set_line_width(min(6, 3 + math.log(coverage, 4)))
        self.draw_element(coords, contig_num + self.BASE_LEVEL)


    #level_num = transcript / contig num, level_sign = 0 for gene leve, -1 for transcript, +1 for contig
    def draw_start_codon(self, coords, level_num = 0, level_sign = 0):
        self.context.set_source_rgb(0.1, 0.9, 0.1)
        self.draw_codon(coords, level_sign * (level_num + self.BASE_LEVEL))

    #level_num = transcript / contig num, level_sign = 0 for gene leve, -1 for transcript, +1 for contig
    def draw_stop_codon(self, coords, level_num = 0, level_sign = 0):
        self.context.set_source_rgb(0.8, 0.1, 0.1)
        self.draw_codon(coords, level_sign * (level_num + self.BASE_LEVEL))



    def set_legend(self, gene_id, chromosome, coordinates):
        pass







