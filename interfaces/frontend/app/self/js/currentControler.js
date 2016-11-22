/*
* @Author: jaume.bonet
* @Date:   2016-11-03 14:28:20
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-04 15:11:11
*/

'use strict';

angular.module('topobuilder')
.controller('currentCtrl', function($scope, CONFIG) {
    $scope.changedValue = function(){
        $scope.current.updateValuesTopo($scope.canvas.getCenter(), CONFIG.interface.proportion);
        $scope.current.setCoords();
        $scope.canvas.renderAll();
    };
    $scope.chagedAngle = function() {
        _.map($scope.current.tpb.tilt, function(num, key){
            if ($scope.current.tpb.tilt[key] > 360) {
                $scope.current.tpb.tilt[key] = $scope.current.tpb.tilt[key] - 360;
            }
            else if ($scope.current.tpb.tilt[key] < 0) {
                $scope.current.tpb.tilt[key] = 360 - Math.abs($scope.current.tpb.tilt[key]);
            }
        });
        $scope.$apply();
    }
    $scope.changedLength = function(){
        if ($scope.current.tpb.length == 0) {
            $scope.current.tpb.length = 1;
        }
        else if ($scope.current.tpb.length > 100) {
            $scope.current.tpb.length = 100;
        }
        else {
            $scope.current.tpb.length = Math.floor($scope.current.tpb.length);
        }
        $scope.$apply();
    };
})
.directive('jbCurrentManager', function(){
    return {
        scope: {current: '=current', canvas: '=canvas'},
        templateUrl:'app/self/static/currentmanager.html',
        controller: 'currentCtrl'
    };
});