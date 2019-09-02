/*
* @Author: jaume.bonet
* @Date:   2016-11-09 07:27:50
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-11 14:52:04
*/

'use strict';

angular.module('topobuilder')
.directive('jbMotifLoader', function(MOTIFS, CONFIG){
    return {
        scope: {viewflag: "="},
        templateUrl:'app/self/static/motifloader.html',
        link: function($scope) {
            $scope.quit = function(){
                $scope.viewflag = false;
            }
        }
    }
})
.controller('MotifCtrl', function($scope, $http, MOTIFS, CONFIG) {
    $scope.motifs = MOTIFS;
    $scope.config  = CONFIG;
    $scope.newmotif = { id: null,
                        segments: Array()};
    $scope.motiffile = null;
    $scope.canvas = null;
    $scope.motif2D = null;

    $scope.fileEmpty = function(){
        return $scope.motiffile == null;
    };
    $scope.segmentIsFull = function( segment ) {
        var a = _.map(segment, function(val, key){ return val == null});
        for (var i = a.length - 1; i >= 0; i--) {
            if (a[i]) return false;
        }
        return true;
    }
    $scope.addSegmentField = function(){
        if ($scope.newmotif.segments.length == 0 || $scope.segmentIsFull($scope.newmotif.segments[$scope.newmotif.segments.length - 1])) {
            $scope.newmotif.segments.push({ini:null, end:null, chain:null, id:null, type:null, split:""});
        }
    };
    $scope.completeMotif = function() {
        if ($scope.motiffile == null) return false;
        var ideval = ($scope.newmotif.id != null && $scope.newmotif.id != "");
        if (!ideval) return false;
        var segeval = $scope.newmotif.segments.length > 0;
        if (!segeval) return false;
        return $scope.segmentIsFull($scope.newmotif.segments[$scope.newmotif.segments.length - 1]);
    }
    $scope.makeMotifView = function() {
        $scope.canvas = new fabric.Canvas('minigridcanvas', { selection: false });
        var midproconf = _.clone($scope.config);
        midproconf.interface.proportion = midproconf.interface.proportion / 1;
        $scope.canvas.makeGridCanvas( midproconf.interface.proportion, false, true );
        $scope.motif2D = makeMotifGroup($scope.newmotif.segments, midproconf).initTopoValues('M');
        $scope.canvas.add($scope.motif2D);
        $scope.motif2D.setLeft($scope.canvas.getCenter().left);
        $scope.motif2D.setTop($scope.canvas.getCenter().top);
        $scope.canvas.renderAll();
    }
    $scope.processMotif = function() {
        var formData = new FormData();
        formData.append('service', 'topobuilder');
        formData.append('file', $scope.motiffile);
        formData.append('query', 'get_motif');
        formData.append('segments', JSON.stringify($scope.newmotif));
        formData.append('jobid', $scope.config.name);
        $http({
            method: 'POST',
            url: 'http://'+$scope.config.interface.host+':'+$scope.config.interface.port+'/query',
            headers: { 'Content-Type': undefined },
            data: formData,
            transformRequest: angular.identity
        })
        .success(function (data) {
            console.log(data);
            $scope.newmotif = data;
            $scope.makeMotifView();
        })
        .error(function (data, status) {
            console.log(data);
            console.log(status);
        });
    }
    $scope.cancelMotif = function() {
        $scope.newmotif = { id: null,
                            segments: Array()};
        $scope.motiffile = null;
        if ($scope.canvas != null) {
            $scope.canvas.clear();
        }
        $scope.canvas = null;
        $scope.motif2D = null;
    }
    $scope.addMotif = function() {
        $scope.motifs.desc.push($scope.newmotif);
        $scope.motifs.objects.push($scope.motif2D);
        $scope.$emit('UpdatedMotifList');
        $scope.cancelMotif();
    }
});