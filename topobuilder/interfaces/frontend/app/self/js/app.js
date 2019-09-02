/*
* @Author: jaume.bonet
* @Date:   2016-10-31 14:29:15
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-11 16:20:37
*/

'use strict';

(function(){
    var app = angular.module('topobuilder', ['ui.bootstrap', 'ngDialog', 'ngSanitize', 'ngFileSaver']);

    app.value('CONFIG', {
        name: null,
        pin: Array(),
        date: Array(),
        hLength: 13,
        eLength: 7,
        heDist: 11,
        hhDist: 12,
        epDist: 5,
        etDist: 5,
        maxlink: 0,
        keeplinkratio: false,
        interface: {
            host: 'localhost',
            port: 8080,
            lockedCanvas: false,
            proportion: 10,
            helix: {
                strokeColor: '#000000',
                strokeWidth: 1,
                fillColor: '#ccccff'
            },
            beta: {
                strokeColor: '#000000',
                strokeWidth: 1,
                fillColor: '#ffcccc'
            },
            cross: {
                strokeColor: '#000000',
                strokeWidth: 1,
                fillColor: '#ccffcc'
            },
            placeholders: {
                strokeColor: '#e6b0aa',
                strokeWidth: 1,
                fillColor: '#fdfefe'
            },
            motif: {
                strokeColor: '#000000',
                strokeWidth: 1,
                fillColor: '#ffffff'
            }
        }
    })
    .value('SSE', {
        layers: [],
        list: []
    })
    .value('MOTIFS', {
        desc: [],
        objects: []
    })
    .factory('configFUNC', function(CONFIG) {
        return {
            connection: function() {
                return 'http://'+CONFIG.interface.host+':'+CONFIG.interface.port;
            },
            projectHasStart: function() {
                return !(_.contains( [null, ""], CONFIG.name ) || CONFIG.pin.length == 0);
            },
            postForm: function( query ) {
                var formData = new FormData();
                formData.append('service', 'topobuilder');
                formData.append('jobid', CONFIG.name);
                formData.append('query', query);
                return formData
            }
        };
    });

    app.controller('MainCtrl', function($scope, CONFIG) {
        $scope.viewMotifLoader = false;
        $scope.noProjectID = function() {
            return CONFIG.name == null || CONFIG.name.length == 0;
        };
    });
    // app.directive('emitLastRepeaterElement', function() {
    //     return function(scope) {
    //         if (scope.$last){
    //             
    //         }
    //     };
    // });

    app.controller('EditorCtrl', function($scope, ngDialog, FileSaver, Blob, CONFIG, MOTIFS, SSE){
        $scope.canvas = BigEditorCanvas;
        $scope.canvas.makeGridCanvas( CONFIG.interface.proportion, true, true );
        $scope.add    = undefined;
        $scope.config = CONFIG;
        $scope.motifs = MOTIFS;
        $scope.active = null;
        $scope.phldr  = Array();

        $scope.add2Canvas = function( newsse ) {
            SSE.list.push(newsse);
            $scope.canvas.add(newsse);
        }

        // $scope.$on('LastRepeaterElement', function() {
        //     for (var index = 0; index < $scope.motifs.desc.length; index++) {
        //         console.log($scope.motifs.desc[index].id)
        //         var mc = new fabric.Canvas( $scope.motifs.desc[index].id, { selection: false });
        //         console.log($scope.motifs.objects[index])
        //         var mtf = fabric.util.object.clone($scope.motifs.objects[index]);
        //         // var mtf = $scope.motifs.objects[index];
        //         // mtf.initialize()
        //         console.log(mtf)
        //         mc.add(mtf);
        //         mtf.center();
        //         mtf.scale(0.5);
        //         mc.renderAll();
        //     }
        // });

        $scope.loadMiniCanvases = function() {
            for (var index = 0; index < $scope.motifs.desc.length; index++) {
                console.log($scope.motifs.desc[index].id)
                var mc = new fabric.Canvas( 'mc' + index, { selection: false });
                console.log($scope.motifs.objects[index])
                var mtf = fabric.util.object.clone($scope.motifs.objects[index]);
                // var mtf = $scope.motifs.objects[index];
                // mtf.initialize()
                console.log(mtf)
                mc.add(mtf);
                mtf.center();
                mtf.scale(0.5);
                mtf.set('evented', false);
                mtf.set('selectable', false);
                mc.renderAll();
                // mc.on('mouse:down', function( index ){
                //     var mtf = fabric.util.object.clone($scope.motifs.objects[index]);
                //     $scope.canvas.add(mtf);
                //     $scope.canvas.renderAll();
                // });
            }
        };
        // $scope.$on('UpdatedMotifList', $scope.loadMiniCanvases());
        $scope.$watch('motifs', $scope.loadMiniCanvases());
        $scope.motif2Canvas = function (index) {
            var mtf = fabric.util.object.clone($scope.motifs.objects[index]);
            $scope.canvas.add(mtf);
            mtf.set('evented', true);
            mtf.set('selectable', true);
            mtf.center();
            mtf.setCoords();
            $scope.canvas.renderAll();
        }
        $scope.addMode = function( mode ) {
            if (!$scope.canvas.getEventedItemsStack().length) {
                var sse = sseMaker(mode, $scope.canvas.getCenter(), CONFIG);
                $scope.add2Canvas(sse);
            }
            else {
                if ($scope.active == null) {
                    swal({ title:"Wait",
                           text: "Structures are added in relation to previous ones.<br>Select a structure in the editor canvas to add a new one.",
                           type: "warning",
                           html: true});
                }
                else {
                    var tmp = ssePlaceHolders($scope.active, mode, CONFIG);
                    _.each(tmp, function(p){
                        var overlap = _.find(_.map(SSE.list, function(s){ return {left:s.left, top:s.top}; }), {left:p.left, top:p.top});
                        if (overlap == undefined) {
                            $scope.canvas.add(p);
                            $scope.phldr.push(p);
                        }
                    });

                }
            }
        }

        // Canvas Mouse & Actions
        $scope.toggleLockCanvas = function() {
            if (CONFIG.interface.lockedCanvas) CONFIG.interface.lockedCanvas = $scope.canvas.unlockCanvas();
            else                               CONFIG.interface.lockedCanvas = $scope.canvas.lockCanvas();
        };
        $scope.centerCanvasContent = function() {
            $scope.canvas.centerEventedItems();
            if ($scope.active) {
                $scope.active.updateTopoValues( $scope.canvas.getCenter(), CONFIG.interface.proportion );
            }
        };
        $scope.getSVG = function() {
            var data = new Blob([$scope.canvas.toSVG()], { type: 'text/plain;charset=utf-8' });
            if (CONFIG.name == null) FileSaver.saveAs(data, 'form_editor.svg');
            else                     FileSaver.saveAs(data, CONFIG.name + '_editor.svg');
        };
        $scope.getPNG = function() {
            var png = $scope.canvas.toDataURL({format: 'png', multiplier: 2 }).replace("data:image/png;base64,","");
            function fixBinary (bin) {
                var length = bin.length;
                var buf = new ArrayBuffer(length);
                var arr = new Uint8Array(buf);
                for (var i = 0; i < length; i++) {
                  arr[i] = bin.charCodeAt(i);
                }
                return buf;
              }
            var data = new Blob([fixBinary(atob(png))], { type: 'image/png' });
            if ($scope.config.name == null) FileSaver.saveAs(data, 'form_editor.png');
            else                            FileSaver.saveAs(data, $scope.config.name + '_editor.png');
        };
        $scope.canvas.on('mouse:down', function(options) {
            _.each($scope.phldr, function(p){
                if (p == $scope.active) {
                    $scope.add2Canvas(placeholder2sse(p, CONFIG));
                }
                $scope.canvas.remove(p);
            });
            $scope.phldr = Array();
        });
        $scope.canvas.on('object:selected', function(options) {
            $scope.active = $scope.canvas.getActiveObject();
            $scope.active.on('selected', function(){
                $scope.active.updateTopoValues( $scope.canvas.getCenter(), CONFIG.interface.proportion );
            });
            $scope.active.on('mouseup', function(){
                $scope.active.updateTopoValues( $scope.canvas.getCenter(), CONFIG.interface.proportion );
                $scope.$apply();
            });
            $scope.$apply();
        });
        $scope.canvas.on('selection:cleared', function(options) {
            $scope.active = null;
            $scope.$apply();
        });
    });

    app.controller('ProjectCtrl', function($scope, CONFIG, SSE) {

        $scope.canvas2Layers = function() {
            return _.map(SSE.list, function(p){return p.tpb});
        }

        $scope.project2JSON = function() {
            var project = {config: _.clone(CONFIG), layers: $scope.canvas2Layers(), motifs: Array()};
            console.log(project);
            return project;
        }
    });

    app.directive('jbEditorManager', function(){
        return {
            scope: false,
            templateUrl:'app/self/static/editormanager.html',
        };
    });

})();