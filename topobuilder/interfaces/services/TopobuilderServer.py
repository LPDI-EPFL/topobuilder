# -*-
# @project: topobuilder
# @file:    TopobuilderServer.py
#
# @author: jaume.bonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.cat
#
# @date:   2016-11-09 11:05:50
#
# @last modified by:   jaume.bonet
# @last modified time: 2016-11-11 16:34:45
#
# -*-
import os
import sys
import json

import TopobuilderUtils

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

class ServiceServer(object):
    def __init__(self):
        self.jobid = None

    def get_jobfile(self):
        return os.path.join(self.jobid, self.jobid + '.json')

    def process(self, options):
        self.jobid = options.forms.get('jobid')
        if options.forms.get('query') != 'check':
            if not os.path.isdir(self.jobid): os.mkdir(self.jobid)
        return getattr(self, options.forms.get('query'))(options)

    def check(self, options):
        if not os.path.isdir(self.jobid):
            return {'exists': False}
        elif not os.path.isfile(self.get_jobfile()):
            return {'exists': False}
        else:
            with open(self.get_jobfile()) as fd:
                data = json.loads(''.join([l.strip() for l in fd.readlines()]))
            return {'exists': True, 'pin': data['config']['pin']}

    def save(self, options):
        try:
            data = byteify(json.loads(options.forms.get('content')))
            content = json.dumps(data, indent=2, separators=(',', ':'))
            print content
            with open(os.path.join(self.get_jobfile()), 'w') as fd:
                fd.write(content)
            return {'success': 'ok'}
        except:
            return {'success': 'ok', 'msg': sys.exc_info()[0]}

    def load(self, options):
        if os.path.isfile(self.get_jobfile()):
            with open(self.get_jobfile()) as fd:
                data = json.loads(''.join([l.strip() for l in fd.readlines()]))
            return {'success': 'ok', 'data': data}
        else: return {'success': 'ko'}

    def fetch(self, options):
        pass


class TopobuilderServer(ServiceServer):
    def __init__(self):
        super(TopobuilderServer, self).__init__()

    def get_motif(self, options):
        upload = options.files.get('file')
        filedump = os.path.join(self.jobid, upload.filename)
        if not os.path.isfile(filedump):
            upload.save(filedump)
        data = TopobuilderUtils.read_motif(filedump, byteify(json.loads(options.forms.get('segments'))))
        return TopobuilderUtils.process_motifs(data)
