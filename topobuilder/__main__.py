# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-28 12:35:27
# @Last modified by:   bonet
# @Last modified time: 02-Sep-2019

import argparse
import os
import sys
import shutil
import getpass
from multiprocessing import Pool
import signal

import topoIO
import utils

from .form.FormFabric import FormFabric
from .interfaces.frontend import serve as front_serve
from .interfaces.backend import serve as back_serve

def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Plot Rosetta's fragments as an image")

    parser.add_argument('-input', dest='input', action='store', help='JSON input')
    parser.add_argument('-shape', dest='shape', action='store_true',
                        help='Print only a shape of the final sketch', default=False)
    parser.add_argument('-hurry', dest='hurry', action='store_true',
                        help='Print only folds fullfilling all constraints', default=False)
    parser.add_argument('-user', dest='user', action='store', help='cluster username',
                        default = getpass.getuser())
    parser.add_argument('-show',  dest='show',  action='store_true',
                        help='Show mode (requires keypress)', default=False)

    parser.add_argument('-X', '-interactive', dest='interactive', action='store_true',
                        help='Open in interactive mode', default=False)

    options = parser.parse_args()
    if not options.interactive:
        options.__dict__["chkpoint"] = options.input.replace(".json", "_chk.json")

    return options

def non_interactive(options):
    try:
        # READING INPUT DATA + MAKE OUTPUT DIR
        print "Setting up the output folder and the initial configuration"
        json_data = topoIO.load_json(options.input)
        if not os.path.isdir(json_data["config"]["name"]):
            os.mkdir(json_data["config"]["name"])
        options.outdir   = json_data["config"]["name"]
        options.chkpoint = os.path.join(json_data["config"]["name"], options.chkpoint)
        utils.form_nomenclator(json_data)  # Set standard ID's for each STATUS

        # PROCESS PDB -> STATUS = 1 (if there are motifs to read)
        print "Reading the motifs (if any)"
        topoIO.read_pdbs(json_data)
        topoIO.print_json(json_data, options.chkpoint)
        if options.show: raw_input("READ PDB. Press for next...")

        # PROCESS MOTIFS -> STATUS = 2 (if there are motifs to read)
        print "Processing the motifs (if any)"
        utils.process_motifs(json_data, options)
        topoIO.print_json(json_data, options.chkpoint)
        if options.show: raw_input("PROCESS MOTIFS. Press for next...")

        # FORM COMBINATORIAL -> STATUS = 3
        print "Building and evaluating combinations"
        FormFabric().build(json_data, options)
        topoIO.print_json(json_data, options.chkpoint)
        htmlfile = os.path.join(json_data["config"]["name"], "combinations.html")
        tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
        shutil.copyfile(os.path.join(tpl, "bio-pv.min.js"),
                        os.path.join(json_data["config"]["name"], "bio-pv.min.js"))
        utils.make_html(json_data, htmlfile)
        if options.show: raw_input("EVALUATE COMBINATIONS. Press for next...")

        # PREPARE ALL FORMS -> STATUS = 4
        print "Preparing and printing the final outputs"
        utils.prepare_forms(json_data, options)
        topoIO.print_json(json_data, options.chkpoint)
        if options.show: raw_input("AVAILABLE FORMS PREPARED FOR FFL. Press for next...")

    except KeyboardInterrupt:
        print "Stopped by the user"

def backend(options):
    print "there"

def interactive(options):
    def init_worker():
        def sig_int(signal_num, frame):
            print 'signal: %s' % signal_num
        signal.signal(signal.SIGINT, sig_int)

    pool = Pool(2, init_worker)
    try:
        front = pool.apply_async(front_serve, [options])
        back  = pool.apply_async(back_serve, [options])
        front.get(sys.maxsize)
        back.get(sys.maxsize)
    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        print "Stopped by the user"


if __name__ == '__main__':

    # GETTING USER OPTIONS
    options = get_options()

    if not options.interactive: non_interactive(options)
    else:                       interactive(options)
