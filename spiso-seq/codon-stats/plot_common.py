#!/usr/bin/env python

import os
import sys
import cairo
import math

class GeneDrawer:
    SPACE_BETWEEN_LEVELS = 10
    BASE_LEVEL = 3
    XMARGIN = 50
    YMARGIN = 50
    CODON_SIZE = 3

    def __init__(self, transcript_num, gene_coords, fname):
        self.GENOME_START_POS = gene_coords[0]
        self.GENOME_END_POS = gene_coords[1]
        self.TOTAL_GENES = transcript_num

        self.CANVAS_WIDTH = min(2000, self.GENOME_END_POS - self.GENOME_START_POS)
        self.CANVAS_HEIGHT = max(300, self.SPACE_BETWEEN_LEVELS * (self.TOTAL_GENES + self.BASE_LEVEL))

        self.WIDTH = self.CANVAS_WIDTH + 2 * self.XMARGIN
        self.HEIGHT = self.CANVAS_HEIGHT + 2 * self.YMARGIN

        self.XSCALE = float(self.CANVAS_WIDTH) / float(self.GENOME_END_POS - self.GENOME_START_POS)
        self.YSCALE = float(self.CANVAS_HEIGHT) / float(self.TOTAL_GENES + self.BASE_LEVEL)

        self.surface = cairo.SVGSurface(fname + ".svg", self.WIDTH, self.HEIGHT)
        self.context = cairo.Context(self.surface)
        self.context.rectangle(0, 0, self.WIDTH, self.HEIGHT)
        self.context.set_source_rgb(1, 1, 1)
        self.context.fill()
        self.context.stroke()


    def draw_element(self, genome_coords, level):
        x1 = self.XMARGIN + (genome_coords[0] - self.GENOME_START_POS)  * self.XSCALE
        x2 = self.XMARGIN + (genome_coords[1] - self.GENOME_START_POS)  * self.XSCALE
        y = self.YMARGIN + level * self.YSCALE
        self.context.move_to(x1, y)
        self.context.line_to(x2, y)
        self.context.stroke()
        #print(x1,x2,y)


    def draw_codon(self, genome_coords, level):
        x = self.XMARGIN + (genome_coords[0] + genome_coords[1] - 2 * self.GENOME_START_POS) * self.XSCALE / 2
        y = self.YMARGIN + level * self.YSCALE
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
        self.draw_element(coords, transcript_num + self.BASE_LEVEL)


    def draw_exon(self, coords, transcript_num):
        self.context.set_source_rgb(0.1, 0.1, 0.9)
        self.context.set_line_width(3)
        self.draw_element(coords, transcript_num + self.BASE_LEVEL)



    #transcript_num = 0 for reference drawing
    def draw_start_codon(self, coords, transcript_num = -1):
        self.context.set_source_rgb(0.1, 0.9, 0.1)
        self.draw_codon(coords, transcript_num + (1 if transcript_num == -1 else self.BASE_LEVEL))


    #transcript_num = -self.BASE_LEVEL for reference drawing
    def draw_stop_codon(self, coords, transcript_num = -1):
        self.context.set_source_rgb(0.8, 0.1, 0.1)
        self.draw_codon(coords, transcript_num + (1 if transcript_num == -1 else self.BASE_LEVEL))

    def draw_coords(self):
        pass

    def set_legend(self, gene_id, chromosome, coordinates):
        pass







