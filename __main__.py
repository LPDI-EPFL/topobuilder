# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-28 12:35:27
# @Last Modified by:   bonet
# @Last Modified time: 2016-05-02 23:59:32

import argparse
import os
import getpass

import topoIO
import utils

from .form.FormFabric import FormFabric


def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Plot Rosetta's fragments as an image")

    parser.add_argument('-input', dest='input', action='store', help='JSON input')
    parser.add_argument('-user', dest='user', action='store', help='castor username',
                        default = getpass.getuser())
    parser.add_argument('-show',  dest='show',  action='store_true',
                        help='Show mode (requires keypress)', default=False)

    options = parser.parse_args()
    options.__dict__["chkpoint"] = options.input.replace(".json", "_chk.json")

    return options


if __name__ == '__main__':

    # GETTING USER OPTIONS
    options = get_options()

    try:
        # READING INPUT DATA + MAKE OUTPUT DIR
        print "Setting up the output folder and the initial configuration"
        json_data = topoIO.load_json(options.input)
        if not os.path.isdir(json_data["config"]["name"]):
            os.mkdir(json_data["config"]["name"])
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
        print "Building and evaluating convinations"
        FormFabric().build(json_data, options)
        topoIO.print_json(json_data, options.chkpoint)
        htmlfile = os.path.join(json_data["config"]["name"], "combinations.html")
        utils.make_html(json_data, htmlfile)
        if options.show: raw_input("EVALUATE COMBINATIONS. Press for next...")

        # PREPARE ALL FORMS -> STATUS = 4
        print "Preparing and printing the final outputs"
        utils.prepare_forms(json_data, options)
        topoIO.print_json(json_data, options.chkpoint)
        if options.show: raw_input("AVAILABLE FORMS PREPARED FOR FFL. Press for next...")

    except KeyboardInterrupt:
        print "Stopped by the user"
