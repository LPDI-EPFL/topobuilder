# -*-
# @project: topobuilder
# @file:    frontend.py
#
# @author: jaume.bonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.cat
#
# @date:   2016-11-08 15:23:46
#
# @last modified by:   jaume.bonet
# @last modified time: 2016-11-08 16:24:24
#
# -*-
import os
import SimpleHTTPServer
import SocketServer
import webbrowser

def serve(options):
    try:
        os.chdir(os.path.normpath(os.path.dirname(__file__)))
        PORT    = 8008
        Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
        httpd   = SocketServer.TCPServer(("", PORT), Handler)
        webbrowser.open_new_tab('http:localhost:{0}/frontend/index.html'.format(PORT))
        httpd.serve_forever()
    except KeyboardInterrupt:
        httpd.shutdown()
