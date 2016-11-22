# -*-
# @project: topobuilder
# @file:    backend.py
#
# @author: jaume.bonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.cat
#
# @date:   2016-11-08 16:30:49
#
# @last modified by:   jaume.bonet
# @last modified time: 2016-11-09 16:37:39
#
# -*-
import bottle
import os
import sys
import json
import importlib

global wserver
wserver = None

@bottle.hook('after_request')
def enable_cors():
    """
    You need to add some headers to each request.
    Don't use the wildcard '*' for Access-Control-Allow-Origin in production.
    """
    bottle.response.headers['Access-Control-Allow-Origin']  = '*'
    bottle.response.headers['Access-Control-Allow-Methods'] = 'PUT, GET, POST, DELETE, OPTIONS'
    bottle.response.headers['Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'

@bottle.route('/', method='GET')
def root():
    return {
        'status': 'ok',
        'services': {'topobuilder': 'TopobuilderServer'}
    }
def init_server(qserver):
    k = 'services'
    sys.path.append(os.path.dirname(__file__))
    m = importlib.import_module( ".".join([k, root()[k][qserver]]), __package__)
    return getattr(m, root()['services'][qserver])()

@bottle.route('/query', method='POST')
def post():
    service = init_server(bottle.request.forms.get('service'))
    return json.dumps(service.process(bottle.request))

def serve(options):
    PORT = 8080
    host = 'localhost'
    bottle.run(host=host, port=PORT)
